############################################
# imports
############################################

import os
import yaml
import warnings
import random
import pandas as pd

from pathlib import Path
from copy import copy, deepcopy
from joblib import cpu_count, Parallel, delayed
from itertools import chain
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..utils._data_parser import (
    check_fasta_format,
    check_tsv_format,
    parse_fasta_header,
)

from ..utils._sequence_design import generate_random_sequence

############################################
# Oligo Database Class
############################################


class OligoDatabase:
    """This class generates all possible oligos that can be designed for a given list of regions (e.g. genes),
    based on the transcriptome or the gene CDS annotation or the whole genome provided as fasta file.

    The header of each sequence must start with '>' and contain the following information:
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
    where the region_id is compulsory and the other fileds are opional.

    Input Format (per sequence):
    >region_id::additional information::chromosome:start-end(strand)
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    Moreover, the database can be saved and loaded to/from a tsv file.

    :param min_oligos_per_region: Minimum number of oligos per region, if lower, region is removed from database, defaults to 0.
    :type min_oligos_per_region: int, optional
    :param write_regions_with_insufficient_oligos: write removed regions to file ``regions_with_insufficient_oligos.txt``, defaults to True
    :type write_regions_with_insufficient_oligos: bool, optional
    :param n_jobs: Number of parallel processes, if set to None all available CPUs are use, defaults to None.
    :type n_jobs: int, optional
    :param dir_output: Output directory, defaults to 'output'.
    :type dir_output: str, optional
    """

    def __init__(
        self,
        min_oligos_per_region: int = 0,
        write_regions_with_insufficient_oligos: bool = True,
        n_jobs: int = None,
        dir_output: str = "output",
    ):
        """Constructor"""
        self.min_oligos_per_region = min_oligos_per_region
        self.write_regions_with_insufficient_oligos = write_regions_with_insufficient_oligos

        self.metadata = {}

        if n_jobs is None:
            n_jobs = cpu_count() - 2
        self.n_jobs = n_jobs

        self.dir_output = os.path.abspath(os.path.join(dir_output, "oligo_database"))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # Initialize databse object
        self.database = {}
        self.oligosets = {}  # will be used later in the gereration of non overlpping sets

        # Initialize the file for regions with insufficient oligos
        if self.write_regions_with_insufficient_oligos:
            self.file_removed_regions = os.path.join(self.dir_output, "regions_with_insufficient_oligos.txt")
            with open(self.file_removed_regions, "a") as handle:
                handle.write(f"Region\tPipeline step\n")

    def load_metadata(self, metadata):
        """Loading metadata information for the oligo database.

        :param metadata: Dictionary or path to yaml file containing metadata of oligo database
        :type dict_metadata: dict or str
        """
        if self.metadata:
            warnings.warn("Metadata not empty! Overwriting metadata with new metadata from file!")

        if type(metadata) is str and os.path.exists(metadata):
            with open(metadata) as handle:
                self.metadata = yaml.safe_load(handle)
        elif type(metadata) is dict:
            self.metadata = metadata
        else:
            raise ValueError("Metadat has icorrect format!")

    def merge_databases(self, database1, database2):
        def _get_sequence_as_key(database):
            database_modified = {region: {} for region in database.keys()}
            for region, values in database.items():
                for oligo_id, oligo_info in values.items():
                    oligo_sequence = oligo_info["sequence"]
                    oligo_info.pop("sequence")
                    database_modified[region][oligo_sequence] = oligo_info
            return database_modified

        def _add_database_content(database_tmp, database_in):
            for region, values in database_in.items():
                for oligo_sequence, oligo_info in values.items():
                    if oligo_sequence in database_tmp[region]:
                        oligo_info_merged = self._collapse_info_for_duplicated_sequences(
                            database_tmp[region][oligo_sequence], oligo_info
                        )
                        database_tmp[region][oligo_sequence] = oligo_info_merged
                    else:
                        database_tmp[region][oligo_sequence] = oligo_info
            return database_tmp

        database_tmp = {region: {} for region in chain(database1.keys(), database2.keys())}
        database_tmp = _add_database_content(database_tmp, _get_sequence_as_key(database1))
        database_tmp = _add_database_content(database_tmp, _get_sequence_as_key(database2))

        database_merged = {region: {} for region in database_tmp.keys()}
        for region, value in database_tmp.items():
            i = 1
            for oligo_sequence, oligo_info in value.items():
                oligo_id = f"{region}::{i}"
                i += 1
                oligo_seq_info = {"sequence": oligo_sequence} | oligo_info
                database_merged[region][oligo_id] = oligo_seq_info

        return database_merged

    def load_sequences_from_database(
        self, file_database: str, region_ids: list[str] = None, database_overwrite: bool = False
    ):
        """Loading a previously generated oligos database and saves it in the ``database`` attribute as a dictionary.
        The order of columns in the database file is:

        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+

        Additional feat. includes additional info from fasta file and additional info computed by the filtering class.

        add sequences from anotehr existing database.
        merge duplicated sequences from the same region, considering star/end coordinates.
        extend the database, not overwrite it.

        :param file_database: File containing oligo database.
        :type file_database: str
        """

        def replace_none_with_null(obj):
            if isinstance(obj, dict):
                return {key: replace_none_with_null(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [replace_none_with_null(item) for item in obj]
            elif obj == "None":
                return None
            else:
                return obj

        if database_overwrite:
            warnings.warn("Overwriting database!")

        if os.path.exists(file_database):
            if not check_tsv_format(file_database):
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        file_tsv_content = pd.read_table(file_database, sep="\t")

        file_tsv_content = file_tsv_content.apply(
            lambda col: col.apply(lambda x: str(x).split("__MATCHSEQ__") if pd.notna(x) else [None])
        )

        database_tmp1 = file_tsv_content.to_dict(orient="records")
        database_tmp2 = defaultdict(dict)
        for entry in database_tmp1:
            entry = replace_none_with_null(entry)
            region_id, oligo_id = entry.pop("region_id")[0], entry.pop("oligo_id")[0]
            entry["sequence"] = entry["sequence"][0]
            database_tmp2[region_id][oligo_id] = entry

        database_tmp2 = dict(database_tmp2)

        if not database_overwrite and self.database:
            database_tmp2 = self.merge_databases(deepcopy(self.database), database_tmp2)

        if region_ids:
            database_tmp2 = self._filter_dabase_for_region(self, database_tmp2, region_ids)

        self.database = database_tmp2

    def load_sequences_from_fasta(
        self, file_fasta: str, region_ids: list[str] = None, database_overwrite: bool = False
    ):
        """add sequences to the database that are stored in a fasta file.
        load those sequences as they are (no processing into smaller windows) but process header.
        merge duplicated sequences from the same region, considering if start/end are the same.
        extend the database, not overwrite it.
        check if sequnce should be target or oligo, i.e. sequence type."""

        if database_overwrite:
            warnings.warn("Overwriting database!")

        if os.path.exists(file_fasta):
            if not check_fasta_format(file_fasta):
                raise ValueError("Fasta file has incorrect format!")
        else:
            raise ValueError("Fasta file does not exist!")

        # read sequences of fasta file
        with open(file_fasta, "r") as handle:
            fasta_sequences = list(SeqIO.parse(handle, "fasta"))

        region_sequences = {}
        for entry in fasta_sequences:
            region, additional_info, coordinates = parse_fasta_header(entry.id)
            oligo_info = coordinates | additional_info
            if region in region_sequences:
                if entry.seq in region_sequences[region]:
                    oligo_info_merged = self._collapse_info_for_duplicated_sequences(
                        region_sequences[region][entry.seq], oligo_info
                    )
                    region_sequences[region][str(entry.seq)] = oligo_info_merged
                else:
                    region_sequences[region][str(entry.seq)] = oligo_info
            else:
                region_sequences[region] = {str(entry.seq): oligo_info}

        database_tmp = {region: {} for region in region_sequences.keys()}
        for region, value in region_sequences.items():
            i = 1
            for oligo_sequence, oligo_info in value.items():
                oligo_id = f"{region}::{i}"
                i += 1
                oligo_seq_info = {"sequence": oligo_sequence} | oligo_info
                database_tmp[region][oligo_id] = oligo_seq_info

        if not database_overwrite and self.databse:
            database_tmp = self.merge_databases(deepcopy(self.database), database_tmp)

        if region_ids:
            database_tmp = self._filter_dabase_for_region(self, database_tmp, region_ids)

        self.database = database_tmp

    def create_sequences_sliding_window(self):
        """create new fasta files with tiled sequences from input sequences via sliding window.
        Load sequences from input fasta file, and run a sliding window over them to create
        tiled sequences with desired length.
        save the generated sequences in fasta file with header information taken from input fasta
        file but adjust coordinates."""
        pass

    def create_sequences_random(
        self,
        filename: str,
        sequence_length: int,
        num_sequences: int,
        base_alphabet_with_probability: list[float] = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    ):
        """create sequences from random sequences with desired length and sequence composition.
        save as fasta file, with appropriate header, e.g. region id and type=random_sequence."""

        def get_sequence_random(sequence_length, base_alphabet_with_probability):
            bases = list(base_alphabet_with_probability.keys())
            sequence = "".join(
                random.choices(
                    bases, weights=[base_alphabet_with_probability[n] for n in bases], k=sequence_length
                )
            )
            return sequence

        file_fasta = os.path.join(self.dir_output, f"{filename}.fna")

        sequences_list = []
        while len(sequences_list) < num_sequences:
            num_missing_sequences = num_sequences - len(sequences_list)
            missing_sequences = list(
                set(
                    [
                        get_sequence_random(sequence_length, base_alphabet_with_probability)
                        for i in range(num_missing_sequences)
                    ]
                )
            )
            sequences_list = list(set(sequences_list + missing_sequences))

        with open(file_fasta, "w") as handle_fasta:
            for i in range(num_sequences):
                header = f"randomsequence{i+1}::regiontype=random_sequence"
                seq = get_sequence_random(sequence_length, base_alphabet_with_probability)
                handle_fasta.write(f">{header}\n{seq}\n")
        return file_fasta

    def write_database(
        self,
        filename: str = "oligo_database",
    ):
        """
        Saves the oligo database as a tsv file.
        The order of the columns is:

        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+

        Additional feat. includes additional info from fasta file and additional info computed by the filtering class.

        :param filename: Database filename prefix, defaults to "oligo_database"
        :type filename: str
        :return: Path to database file (tsv file).
        :rtype: str
        """
        file_metadata = os.path.join(self.dir_output, filename + ".yaml")

        with open(file_metadata, "w") as handle:
            yaml.safe_dump(self.metadata, handle, sort_keys=True, default_flow_style=False)

        file_database = os.path.join(self.dir_output, filename + ".tsv")
        file_tsv_content = []

        for region_id, oligo_dict in self.database.items():
            for oligo_id, oligo_attributes in oligo_dict.items():
                entry = {"region_id": region_id, "oligo_id": oligo_id}
                entry.update(oligo_attributes)
                # merge all entries where the sequences match into one entry
                entry = {
                    key: "__MATCHSEQ__".join(map(str, values)) if isinstance(values, list) else values
                    for key, values in entry.items()
                }
                file_tsv_content.append(entry)

        file_tsv_content = pd.DataFrame(data=file_tsv_content)
        file_tsv_content.to_csv(file_database, sep="\t", index=False)

        return file_database, file_metadata

    def _collapse_info_for_duplicated_sequences(self, oligo_info1, oligo_info2):
        # Combine the dictionaries
        oligo_info = defaultdict(list)

        for d in (oligo_info1, oligo_info2):
            for key, values in d.items():
                oligo_info[key].extend(values)

        # Convert defaultdict back to a regular dictionary
        oligo_info = dict(oligo_info)

        return oligo_info

    def _filter_dabase_for_region(self, database, region_ids):
        # if a list of region ids is provided filter dict by this list
        # in addition check if all regions provided in region list exist
        # in the regions retreived from the fasta file
        keys = list(database.keys())
        for key in keys:
            if key not in region_ids:
                database.pop(key)
        for region_id in region_ids:
            if region_id not in keys:
                warnings.warn(f"Region {region_id} not available in reference file.")
                if self.write_regions_with_insufficient_oligos:
                    with open(self.file_removed_regions, "a") as hanlde:
                        hanlde.write(f"{region_id}\t{'Not in Annotation'}\n")
        return database

    ############################################
    # OLD
    ############################################
    def create_database_from_sequences(
        self,
        file_fasta: str,
        oligo_length_min: int,
        oligo_length_max: int,
        region_ids: list[str] = None,
    ):
        """
        Creates the database containing all the oligo sequence extracted form the given sequences (fasta file)
        and belonging the the specified region. If no list of region IDs is specified then all regions available
        in the fasta file will be used. The database created is not written automatically to the disk,
        the ``save_oligo_database`` method has to be called separately.

        :param file_fasta: Path to the fasta file, if None it is only possible to read a database.
        :type file_fasta: str
        :param oligo_length_min: Minimal length of oligo nucleotide.
        :type oligo_length_min: int
        :param oligo_length_max: Maximal length of oligo nucleotide.
        :type oligo_length_max: int
        :param region_ids: List of regions for which the oligos should be generated, defaults to None
        :type region_ids: list of str, optional
        """
        # check if files exist and are in correct format
        if os.path.exists(file_fasta):
            if not check_fasta_format(file_fasta):
                raise ValueError("Fasta file has incorrect format!")
        else:
            raise ValueError("Fasta file does not exist!")

        self.oligo_length_min = oligo_length_min
        self.oligo_length_max = oligo_length_max

        # read sequences of fasta file
        with open(file_fasta, "r") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))

        # group sequences in a dictionary by regions id
        region_sequences = {}
        for entry in sequences:
            region, _, _ = parse_fasta_header(entry.id)
            if region in region_sequences:
                region_sequences[region].append(entry)
            else:
                region_sequences[region] = [entry]

        # if a list of region ids is provided filter dict by this list
        # in addition check if all regions provided in region list exist
        # in the regions retreived from the fasta file
        if region_ids:
            keys = copy(list(region_sequences.keys()))
            for key in keys:
                if key not in region_ids:
                    region_sequences.pop(key)
            for region_id in region_ids:
                if region_id not in keys:
                    warnings.warn(f"Region {region_id} not available in reference file.")
                    if self.write_regions_with_insufficient_oligos:
                        with open(self.file_removed_regions, "a") as hanlde:
                            hanlde.write(f"{region_id}\t{'Not in Annotation'}\n")

        # create the database dictionary
        region_ids = region_sequences.keys()
        database = {}
        results = Parallel(n_jobs=self.n_jobs)(
            delayed(self._get_oligos_info)(region_id, region_sequences[region_id]) for region_id in region_ids
        )
        for result in results:
            database.update(result)

        self.database = database

    def _get_oligos_info(
        self,
        region_id,
        region_sequences,
    ):
        """Generate oligo sequences and save all corresponding information, like chromosome,
        start, end, strand, length, additional information from fasat header.

        :param region_id: ID of region, used as first layer key in the dict.
        :type region_id: str
        :param region_sequences: All sequence records for the region.
        :type region_sequences: Bio.Seq.Seq
        :return: Oligo Database
        :rtype: dict
        """

        def _merge(dicts):
            merged_dict = {}
            for dict in dicts:
                for key in dict.keys():
                    if key in merged_dict:
                        merged_dict[key].append(dict[key])
                    else:
                        merged_dict[key] = dict[key]

            return merged_dict

        database_entries = {}  # assiciate each id to the additional features
        tmp = {}  # associate each sequence to its id

        # parse the exon fasta sequence file
        for region_sequence in region_sequences:
            seq = region_sequence.seq
            _, add_inf, coordinates = parse_fasta_header(region_sequence.id)

            # should all have the same chromosome and strand entry
            # if no information provided in fasta file, those entries are None
            chromosome = coordinates["chromosome"][0]
            strand = coordinates["strand"][0]

            list_of_coordinates = []
            # if the fasta header DOES NOT contain coordinates information
            if coordinates["start"][0] is None:
                list_of_coordinates = [None for i in range(len(seq))]
            # if the fasta header DOES contain coordinates information
            else:
                for i in range(len(coordinates["start"])):
                    list_of_coordinates.extend(
                        list(
                            range(
                                # coordinates in fasta file use 1-base indixing, which go from 1 (for base 1) to n (for base n)
                                coordinates["start"][i],
                                # range produces values until end-1 (e.g. range(10) goes until 9) -> add +1
                                coordinates["end"][i] + 1,
                            )
                        )
                    )
            # sort reverse on minus strand becaus the sequence is translated into
            # the reverse complement by fasta -strand option
            if strand == "-":
                list_of_coordinates.reverse()

            # generate oligos for length range
            for oligo_length in range(self.oligo_length_min, self.oligo_length_max + 1):
                if len(seq) > oligo_length:
                    number_oligos = len(seq) - (oligo_length - 1)
                    oligos = [seq[i : i + oligo_length] for i in range(number_oligos)]

                    for i in range(number_oligos):
                        oligo_id = f"{region_id}::{i}"
                        oligo = oligos[i]
                        oligo_start_end = [
                            list_of_coordinates[i],
                            list_of_coordinates[(i + oligo_length - 1)],
                        ]
                        # sort coordinates for oligos on minus strand
                        start = min(oligo_start_end)  # 1-base index
                        end = max(oligo_start_end)
                        oligo_attributes = {
                            "sequence": oligo,
                            "chromosome": chromosome,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "length": oligo_length,
                            "additional_information": add_inf,
                        }
                        #### this generated some problems, weird output see test_database_notebook
                        if oligo in tmp:
                            tmp[oligo][oligo_id] = oligo_attributes
                        else:
                            tmp[oligo] = {oligo_id: oligo_attributes}

        database_entries = {}
        # collapse sequence entries like chromosome, strand etc into one entry
        for oligo_ids in tmp.values():
            print(oligo_ids)
            # use first id as final id for sequence
            oligo_id = list(oligo_ids.keys())[0]
            attributes_collapsed = {}
            for key in ["sequence", "chromosome", "start", "end", "strand", "length"]:
                attributes_collapsed[key] = list(set([attributes[key] for attributes in oligo_ids.values()]))
                if len(attributes_collapsed[key]) == 1 and key not in ["start", "end"]:
                    attributes_collapsed[key] = attributes_collapsed[key][0]
            # can't be collapsed caus sometimes contains lists of lists
            if attributes_collapsed["chromosome"] is None:
                attributes_collapsed["additional_information"] = []
            else:
                #### this generated some problems, weird output see test_database_notebook
                dicts_add_inf = [attributes["additional_information"] for attributes in oligo_ids.values()]
                if len(dicts_add_inf) > 1:
                    print(dicts_add_inf)
                attributes_collapsed["additional_information"] = _merge(dicts_add_inf)
            database_entries[oligo_id] = attributes_collapsed
        database = {region_id: database_entries}

        return database

    def to_sequence_list(self) -> list:
        """Converts the database into a list of sequences

        Returns:
            list of str: list of all sequences contained in the database
        """
        sequences = []
        for region_id, oligo_dict in self.database.items():
            for oligo_id, oligo_attributes in oligo_dict.items():
                entry = {"region_id": region_id, "oligo_id": oligo_id}
                entry.update(oligo_attributes)
                sequences.append(str(entry["sequence"]))

        return sequences

    def write_fasta_from_database(
        self,
        filename: str = "oligo_database",
    ):
        """Write sequences stored in database to fasta file.

        :param filename: Database filename prefix, defaults to "oligo_database"
        :type filename: str, optional
        :return: Path to fasta file.
        :rtype: str
        """
        file_fasta = os.path.join(self.dir_output, f"{filename}.fna")
        database = deepcopy(self.database)
        output_fasta = []

        with open(file_fasta, "w") as handle_fasta:
            for region_id, oligo in database.items():
                for oligo_id, oligo_attributes in oligo.items():
                    oligo_id = oligo_id.split(";")
                    oligo_length = oligo_attributes["length"]
                    for i in range(len(oligo_id)):
                        seq_description = f"{oligo_length}bp oligo nucleotide"
                        # parameter: seq, seq_id, seq_name, seq_description
                        output_fasta.append(
                            SeqRecord(
                                oligo_attributes["sequence"],
                                oligo_id[i],
                                oligo_id[i],
                                seq_description,
                            )
                        )  # write the sequence in the fasta file
            SeqIO.write(output_fasta, handle_fasta, "fasta")

        return file_fasta

    def write_oligosets(self, folder: str = "oligosets"):
        """Writes the data structure ``self.oligosets`` in a series of files, each contains the oligosets for one gene and is called "{gene}_oligosets.tsv".
        The files will be stored in a subdirectory of ``self.dir_output``.

        :param dir_oligosets: Subdirectory to store oligosets, defaults to "oligosets"
        :type dir_oligosets: str, optional
        :return: Path to oligosets directory.
        :rtype: str
        """
        dir_oligosets = os.path.join(self.dir_output, folder)
        Path(dir_oligosets).mkdir(parents=True, exist_ok=True)

        for region_id in self.oligosets.keys():
            file_oligosets = os.path.join(dir_oligosets, f"{region_id}_oligosets.tsv")
            self.oligosets[region_id].to_csv(file_oligosets, sep="\t", index=False)

        return dir_oligosets

    def remove_regions_with_insufficient_oligos(
        self,
        pipeline_step: str,
    ):
        """Deletes from the ``oligo_DB`` the regions (e.g. genes) which have less than ``min_oligos_per_region`` oligos,
        and optionally writes them in a file with the name of the step of the pipeline at which they have been deleted.

        :param pipeline_step: Step in the pipeline that lead to the removal of the region.
        :type pipeline_step: str
        """
        regions = copy(list(self.database.keys()))
        for region in regions:
            if len(list(self.database[region].keys())) <= self.min_oligos_per_region:
                del self.database[region]
                if region in self.oligosets:
                    del self.oligosets[region]
                if self.write_regions_with_insufficient_oligos:
                    with open(self.file_removed_regions, "a") as hanlde:
                        hanlde.write(f"{region}\t{pipeline_step}\n")
