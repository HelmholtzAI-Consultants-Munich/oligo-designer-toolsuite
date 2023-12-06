############################################
# imports
############################################

import os
import yaml
import warnings
import random
import pandas as pd

from pathlib import Path

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..utils._data_parser import (
    check_fasta_format,
    check_tsv_format,
    parse_fasta_header,
)

from ..utils._database_processor import collapse_info_for_duplicated_sequences, merge_databases


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
        dir_output: str = "output",
    ):
        """Constructor"""
        self.min_oligos_per_region = min_oligos_per_region
        self.write_regions_with_insufficient_oligos = write_regions_with_insufficient_oligos

        self.metadata = {}

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

    def load_database(
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
        region_ids = self._check_if_list(region_ids)

        if database_overwrite:
            warnings.warn("Overwriting database!")

        if os.path.exists(file_database):
            if not check_tsv_format(file_database):
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        file_tsv_content = pd.read_table(file_database, sep="\t")

        file_tsv_content = file_tsv_content.apply(
            lambda col: col.apply(lambda x: x if pd.notna(x) else [None])
        )

        file_tsv_content = file_tsv_content.apply(
            lambda col: col.apply(
                lambda x: eval(x) if isinstance(x, str) and x.startswith("[") and x.endswith("]") else x
            )
        )

        database_tmp1 = file_tsv_content.to_dict(orient="records")
        database_tmp2 = defaultdict(dict)
        for entry in database_tmp1:
            region_id, oligo_id = entry.pop("region_id"), entry.pop("oligo_id")
            entry["sequence"] = entry["sequence"]
            database_tmp2[region_id][oligo_id] = entry

        database_tmp2 = dict(database_tmp2)

        if not database_overwrite and self.database:
            database_tmp2 = merge_databases(self.database, database_tmp2)

        if region_ids:
            database_tmp2 = self._filter_dabase_for_region(database_tmp2, region_ids)
            self._check_if_region_in_database(database_tmp2, region_ids)

        self.database = database_tmp2

    def load_sequences_from_fasta(
        self, file_fasta_in: str, region_ids: list[str] = None, database_overwrite: bool = False
    ):
        """add sequences to the database that are stored in a fasta file.
        load those sequences as they are (no processing into smaller windows) but process header.
        merge duplicated sequences from the same region, considering if start/end are the same.
        extend the database, not overwrite it.
        check if sequnce should be target or oligo, i.e. sequence type."""

        region_ids = self._check_if_list(region_ids)

        if database_overwrite:
            warnings.warn("Overwriting database!")

        fasta_sequences = self._read_sequences(file_fasta_in, region_ids)

        region_sequences = {}
        for entry in fasta_sequences:
            region, additional_info, coordinates = parse_fasta_header(entry.id)
            oligo_info = coordinates | additional_info
            if region in region_sequences:
                if entry.seq in region_sequences[region]:
                    oligo_info_merged = collapse_info_for_duplicated_sequences(
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

        if not database_overwrite and self.database:
            database_tmp = merge_databases(self.database, database_tmp)

        # add this step to log regions which are not available in database
        if region_ids:
            self._check_if_region_in_database(database_tmp, region_ids)

        self.database = database_tmp

    def create_sequences_sliding_window(
        self, filename_out: str, file_fasta_in: list[str], sequence_length: int, region_ids: list[str] = None
    ):
        """create new fasta files with tiled sequences from input sequences via sliding window.
        Load sequences from input fasta file, and run a sliding window over them to create
        tiled sequences with desired length.
        save the generated sequences in fasta file with header information taken from input fasta
        file but adjust coordinates."""

        def get_sliding_window_sequence(entry, sequence_length, handle_fasta):
            entry_sequence = entry.seq
            region, additional_info, coordinates = parse_fasta_header(
                header=entry.id, parse_additional_info=False
            )
            # chromosome and strand information is the same for an entry but parsed as
            # list in the parse_fasta_header, hence, we take only the first element for each
            # if no information provided in fasta file, those entries are None
            chromosome = coordinates["chromosome"][0]
            strand = coordinates["strand"][0]
            list_of_coordinates = []
            # if the fasta header DOES NOT contain coordinates information
            if chromosome is None:
                list_of_coordinates = [None for i in range(len(seq))]
            # if the fasta header DOES contain coordinates information
            else:
                # coordinates in fasta file use 1-base indixing, which go from 1 (for base 1) to n (for base n)
                # range produces values until end-1 (e.g. range(10) goes until 9) -> add +1
                # the start and end coordinates can be a list for regions spanning split sequences (e.g. exon junctions, UTRs)
                for i in range(len(coordinates["start"])):
                    list_of_coordinates.extend(
                        list(
                            range(
                                coordinates["start"][i],
                                coordinates["end"][i] + 1,
                            )
                        )
                    )

            # sort reverse on minus strand becaus the sequence is translated into the reverse complement by fasta -strand option
            if strand == "-":
                list_of_coordinates.reverse()

            # generate sequences with sliding window and write to fasta file (use lock to ensure that hat only one process can write to the file at any given time)
            if len(entry_sequence) > sequence_length:
                num_sequences = len(entry_sequence) - (sequence_length - 1)
                sequences = [entry_sequence[i : i + sequence_length] for i in range(num_sequences)]

                for i in range(num_sequences):
                    seq = sequences[i]
                    seq_start_end = [
                        list_of_coordinates[i],
                        list_of_coordinates[(i + sequence_length - 1)],
                    ]
                    # sort coordinates for oligos on minus strand
                    start_seq = min(seq_start_end)  # 1-base index
                    end_seq = max(seq_start_end)

                    header = f"{region}::{additional_info}::{chromosome}:{start_seq}-{end_seq}({strand})"
                    handle_fasta.write(f">{header}\n{seq}\n")

        region_ids = self._check_if_list(region_ids)
        file_fasta_in = self._check_if_list(file_fasta_in)

        file_fasta_out = os.path.join(self.dir_output, f"{filename_out}.fna")
        with open(file_fasta_out, "w") as handle_fasta:
            for file_fasta in file_fasta_in:
                fasta_sequences = self._read_sequences(file_fasta, region_ids)
                # did not parallize this function because workers would be writing simultaneously to the same file
                # which could lead to data corruption. I tried using the file_lock=Lock() function from multiprocessing
                # as input to get_sliding_window_sequence(file_lock) and then use with file_lock: before writing to the file
                # but this lead to an pickle error from joblib which I could not solve, so I left it like that for now
                for entry in fasta_sequences:
                    get_sliding_window_sequence(entry, sequence_length, handle_fasta)

        return file_fasta_out

    def create_sequences_random(
        self,
        filename_out: str,
        length_sequences: int,
        num_sequences: int,
        name_sequences: str = "randomsequence",
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

        # make sure to have n unique sequences
        sequences_list = []
        while len(sequences_list) < num_sequences:
            num_missing_sequences = num_sequences - len(sequences_list)
            missing_sequences = list(
                set(
                    [
                        get_sequence_random(length_sequences, base_alphabet_with_probability)
                        for i in range(num_missing_sequences)
                    ]
                )
            )
            sequences_list = list(set(sequences_list + missing_sequences))

        file_fasta_out = os.path.join(self.dir_output, f"{filename_out}.fna")

        with open(file_fasta_out, "w") as handle_fasta:
            for seq in sequences_list:
                handle_fasta.write(f">{name_sequences}::regiontype=random_sequence\n{seq}\n")
        return file_fasta_out

    def remove_regions_with_insufficient_oligos(self, pipeline_step: str):
        """Deletes from the ``oligo_DB`` the regions (e.g. genes) which have less than
        ``min_oligos_per_region`` oligos, and optionally writes them in a file with the
        name of the step of the pipeline at which they have been deleted.

        :param pipeline_step: Step in the pipeline that lead to the removal of the region.
        :type pipeline_step: str
        """
        regions_to_remove = [
            region for region, oligos in self.database.items() if len(oligos) <= self.min_oligos_per_region
        ]

        for region in regions_to_remove:
            del self.database[region]
            self.oligosets.pop(region, None)

        if self.write_regions_with_insufficient_oligos and regions_to_remove:
            with open(self.file_removed_regions, "a") as handle:
                handle.write("\n".join(f"{region}\t{pipeline_step}" for region in regions_to_remove) + "\n")

    def get_sequence_list(self):
        """Converts the database into a list of sequences

        Returns:
            list of str: list of all sequences contained in the database
        """
        sequences = [
            str(oligo_attributes["sequence"])
            for region_id, oligo_dict in self.database.items()
            for oligo_id, oligo_attributes in oligo_dict.items()
        ]

        return sequences

    def calculate_num_targeted_transcripts(self):
        for region_id, oligo_dict in self.database.items():
            for oligo_id, oligo_attributes in oligo_dict.items():
                if "transcript_id" in oligo_attributes:
                    transcript_ids = oligo_attributes["transcript_id"]
                    oligo_dict[oligo_id]["num_targeted_transcripts"] = len(
                        set(
                            item
                            for sublist in (
                                transcript_ids if isinstance(transcript_ids[0], list) else [transcript_ids]
                            )
                            for item in sublist
                        )
                    )
                else:
                    oligo_dict[oligo_id]["num_targeted_transcripts"] = 0

    def save_database(
        self,
        region_ids: list[str] = None,
        filename_out: str = "oligo_database",
    ):
        """
        Saves the oligo database as a tsv file.
        The order of the columns is:

        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+

        Additional feat. includes additional info from fasta file and additional info computed by the filtering class.

        :param filename_out: Database filename_out prefix, defaults to "oligo_database"
        :type filename_out: str
        :return: Path to database file (tsv file).
        :rtype: str
        """
        if region_ids:
            region_ids = self._check_if_list(region_ids)
        else:
            region_ids = self.database.keys()

        file_metadata = os.path.join(self.dir_output, filename_out + ".yaml")

        with open(file_metadata, "w") as handle:
            yaml.safe_dump(self.metadata, handle, sort_keys=True, default_flow_style=False)

        file_database = os.path.join(self.dir_output, filename_out + ".tsv")
        file_tsv_content = []

        for region_id, oligo_dict in self.database.items():
            if region_id in region_ids:
                for oligo_id, oligo_attributes in oligo_dict.items():
                    entry = {"region_id": region_id, "oligo_id": oligo_id}
                    entry.update(oligo_attributes)
                    file_tsv_content.append(entry)

        file_tsv_content = pd.DataFrame(data=file_tsv_content)
        file_tsv_content.to_csv(file_database, sep="\t", index=False)

        return file_database, file_metadata

    def write_fasta_from_database(self, filename: str = "oligo_database"):
        """Write sequences stored in database to fasta file.

        :param filename: Database filename prefix, defaults to "oligo_database"
        :type filename: str, optional
        :return: Path to fasta file.
        :rtype: str
        """
        file_fasta = os.path.join(self.dir_output, f"{filename}.fna")
        output_fasta = []

        with open(file_fasta, "w") as handle_fasta:
            for region_id, oligo in self.database.items():
                for oligo_id, oligo_attributes in oligo.items():
                    seq_record = SeqRecord(
                        Seq(oligo_attributes["sequence"]),
                        id=oligo_id,
                        name=oligo_id.split("::")[0],
                        description="oligonucleotide",
                    )
                    output_fasta.append(seq_record)

            SeqIO.write(output_fasta, handle_fasta, "fasta")

        return file_fasta

    def write_oligosets(self, foldername_out: str = "oligo_sets"):
        """Writes the data structure ``self.oligosets`` in a series of files, each contains the oligosets for one gene and is called "{gene}_oligosets.tsv".
        The files will be stored in a subdirectory of ``self.dir_output``.

        :param dir_oligosets: Subdirectory to store oligosets, defaults to "oligosets"
        :type dir_oligosets: str, optional
        :return: Path to oligosets directory.
        :rtype: str
        """
        dir_oligosets = os.path.join(self.dir_output, foldername_out)
        Path(dir_oligosets).mkdir(parents=True, exist_ok=True)

        for region_id in self.oligosets.keys():
            file_oligosets = os.path.join(dir_oligosets, f"{region_id}_oligosets.tsv")
            self.oligosets[region_id].to_csv(file_oligosets, sep="\t", index=False)

        return dir_oligosets

    def _read_sequences(self, file_fasta_in, region_ids):
        # check if files exist and are in correct format
        if os.path.exists(file_fasta_in):
            if not check_fasta_format(file_fasta_in):
                raise ValueError("Fasta file has incorrect format!")
        else:
            raise ValueError("Fasta file does not exist!")
        # remove undesired regions already at the beginning
        if region_ids:
            fasta_sequences = []
            with open(file_fasta_in, "r") as handle:
                for entry in SeqIO.parse(handle, "fasta"):
                    region, _, _ = parse_fasta_header(entry.id, parse_additional_info=False)
                    if region in region_ids:
                        fasta_sequences.append(entry)
        else:
            with open(file_fasta_in, "r") as handle:
                fasta_sequences = list(SeqIO.parse(handle, "fasta"))

        return fasta_sequences

    def _check_if_region_in_database(self, database, region_ids):
        # check if all regions provided in region list exist
        # in the regions retreived from the fasta file
        keys = list(database.keys())
        for region_id in region_ids:
            if region_id not in keys:
                warnings.warn(f"Region {region_id} not available in reference file.")
                if self.write_regions_with_insufficient_oligos:
                    with open(self.file_removed_regions, "a") as hanlde:
                        hanlde.write(f"{region_id}\t{'Not in Annotation'}\n")

    def _filter_dabase_for_region(self, database, region_ids):
        # if a list of region ids is provided filter dict by this list
        for key in database.keys():
            if key not in region_ids:
                database.pop(key)
        return database

    def _check_if_list(self, obj):
        if obj:
            obj = [obj] if not isinstance(obj, list) else obj
        return obj
