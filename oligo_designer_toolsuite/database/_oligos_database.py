############################################
# imports
############################################

import os
import warnings
import pandas as pd

from pathlib import Path
from copy import copy, deepcopy
from joblib import cpu_count, Parallel, delayed

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..utils._data_parser import (
    check_fasta_format,
    check_tsv_format,
    parse_fasta_header,
)

############################################
# Oligo Database Class
############################################


class OligoDatabase:
    """This class generates all possible oligos that can be designed for a given list of regions (e.g. genes),
    based on the transcriptome or the gene CDS annotation or the whole genome provided as fasta file. The header of 
    each sequence must start with '>' and contain the following information: 
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand).

    Input Format (per sequence):
    >region_id::additional information::chromosome:start-end(strand)
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC
    
    Moreover, the database can be saved and loaded to/from a tsv file.

    Source, species, annotation_release and genome_assembly are set to 'unknown' if thay are not given as input.

    :param file_fasta: Path to the fasta file.
    :type file_fasta: str
    :param oligo_length_min: Minimal length of oligo nucleotide.
    :type oligo_length_min: int
    :param oligo_length_max: Maximal length of oligo nucleotide.
    :type oligo_length_max: int
    :param min_oligos_per_region: Minimum number of oligos per region, if lower, region is removed from database, defaults to 0.
    :type min_oligos_per_region: int, optional
    :param source: Source of annotations, e.g. NCBI, defaults to None.
    :type source: str, optional
    :param species: Species of annotation, e.g. Homo_sapiens, defaults to None.
    :type species: str, optional
    :param annotation_release: Release number of annotation, e.g. 110, defaults to None.
    :type annotation_release: str, optional
    :param genome_assembly: Genome assembly of annotation, e.g. GRCh38, defaults to None.
    :type genome_assembly: str, optional
    :param n_jobs: Number of parallel processes, if set to None all available CPUs are use, defaults to None.
    :type n_jobs: int, optional
    :param dir_output: Output directory, defaults to 'output'.
    :type dir_output: str, optional
    """

    def __init__(
        self,
        file_fasta: str,
        oligo_length_min: int,
        oligo_length_max: int,
        min_oligos_per_region: int = 0,
        source: str = None,
        species: str = None,
        annotation_release: str = None,
        genome_assembly: str = None,
        n_jobs: int = None,
        dir_output: str = "output",
    ):
        """Constructor"""
        self.file_fasta = file_fasta
        if os.path.exists(self.file_fasta):
            if not check_fasta_format(self.file_fasta):
                raise ValueError("Fasta file has incorrect format!")
        else:
            raise ValueError("Fasta file does not exist!")

        self.oligo_length_min = oligo_length_min
        self.oligo_length_max = oligo_length_max
        self.min_oligos_per_region = min_oligos_per_region

        if source is None:
            source = "custom"
            warnings.warn(f"No source defined. Using default source {source}!")

        if species is None:
            species = "unknown"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "unknown"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        if genome_assembly is None:
            genome_assembly = "unknown"
            warnings.warn(
                f"No genome assembly defined. Using default genome assembly {genome_assembly}!"
            )

        self.source = source
        self.species = species
        self.annotation_release = annotation_release
        self.genome_assembly = genome_assembly

        if n_jobs is None:
            n_jobs = cpu_count()
        self.n_jobs = n_jobs

        self.dir_output = os.path.abspath(os.path.join(dir_output, "oligo_database"))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # Initialize databse object
        self.database = {}
        self.oligosets = (
            {}
        )  # will be used later in the gereration of non overlpping sets

        # Initialize the file for regions with insufficient oligos
        self.file_removed_regions = os.path.join(
            dir_output, "regions_with_insufficient_oligos.txt"
        )
        with open(self.file_removed_regions, "a") as handle:
            handle.write(f"Region\tPipeline step\n")

    def create_oligo_database(
        self,
        region_ids: list[str] = None,
    ):
        """
        Creates the database containing all the oligo sequence extracted form the given sequences (fasta file)
        and belonging the the specified region. If no list of region IDs is specified then all regions available
        in the fasta file will be used. The database created is not written automatically to the disk,
        the ``save_oligo_database`` method has to be called separately.

        :param region_ids: List of regions for which the oligos should be generated, defaults to None
        :type region_ids: list of str, optional
        """
        with open(self.file_fasta, "r") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))

        region_sequences = {}
        for entry in sequences:
            region, _, _ = parse_fasta_header(entry.id)
            if region in region_sequences:
                region_sequences[region].append(entry)
            else:
                region_sequences[region] = [entry]

        if region_ids:
            keys = copy(list(region_sequences.keys()))
            for key in keys:
                if key not in region_ids:
                    region_sequences.pop(key)

        region_ids = region_sequences.keys()

        database = {}
        results = Parallel(n_jobs=self.n_jobs)(
            delayed(self._get_oligos_info)(region_id, region_sequences[region_id])
            for region_id in region_ids
        )
        for result in results:
            database.update(result)

        self.database = database

    def load_oligo_database(self, file_database: str):
        """
        Loads a previously generated oligos database and saves it in the ``database`` attribute as a dictionary.
        The order of columns is :

        +-----------+----------+----------------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | oligo_sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------------+------------+-------+-----+--------+--------+------------------+

        Additional feat. includes additional info from fasta file and additional info computed by the filtering class.

        :param file_tsv: Path of tsv file that contains the oligo database.
        :type file_tsv: str
        """
        if os.path.exists(file_database):
            if not check_tsv_format(file_database):
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        file_tsv_content = pd.read_table(file_database, sep="\t")
        cols = file_tsv_content.columns

        file_tsv_content.sequence = file_tsv_content.apply(
            lambda row: Seq(str(row.sequence)), axis=1
        )
        file_tsv_content.chromosome = file_tsv_content.chromosome.astype("str")
        file_tsv_content.start = file_tsv_content.start.str.split(";")
        file_tsv_content.start = file_tsv_content.apply(
            lambda row: list(map(int, row.start)), axis=1
        )
        file_tsv_content.end = file_tsv_content.end.str.split(";")
        file_tsv_content.end = file_tsv_content.apply(
            lambda row: list(map(int, row.end)), axis=1
        )
        file_tsv_content.length = file_tsv_content.length.astype("int")
        file_tsv_content.additional_information_fasta = (
            file_tsv_content.additional_information_fasta.str.split("__MATCHSEQ__")
        )

        database = {}
        for region in file_tsv_content.region_id.unique():
            database[region] = {
                row[1]: {cols[i]: row[i] for i in range(2, len(cols))}
                for row in list(
                    zip(*[file_tsv_content[col] for col in file_tsv_content])
                )
            }

        self.database = database

    def save_oligo_database(
        self,
        filename: str = "oligo_database",
    ):
        """
        Saves the oligo database as a tsv file.
        The order of the columns is:

        +-----------+----------+----------------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | oligo_sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------------+------------+-------+-----+--------+--------+------------------+

        Additional feat. includes additional info from fasta file and additional info computed by the filtering class.

        :param filename: Database filename prefix, defaults to "oligo_database"
        :type filename: str
        :return: Path to database file (tsv file).
        :rtype: str
        """
        file_database = os.path.join(self.dir_output, f"{filename}.tsv")
        file_tsv_content = []

        database = deepcopy(self.database)

        for region_id, oligo_dict in database.items():
            for oligo_id, oligo_attributes in oligo_dict.items():
                entry = {"region_id": region_id, "oligo_id": oligo_id}
                entry.update(oligo_attributes)
                entry["sequence"] = str(entry["sequence"])
                entry["start"] = ";".join(list(map(str, entry["start"])))
                entry["end"] = ";".join(list(map(str, entry["end"])))
                entry["additional_information_fasta"] = "__MATCHSEQ__".join(
                    entry["additional_information_fasta"]
                )
                file_tsv_content.append(pd.DataFrame(entry, index=[0]))

        file_tsv_content = pd.concat(file_tsv_content)
        file_tsv_content.to_csv(file_database, sep="\t", index=False)

        return file_database

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

    def write_oligosets(self, dir_oligosets: str = "oligosets"):
        """Writes the data structure ``self.oligosets`` in a series of files, each contains the oligosets for one gene and is called "{gene}_oligosets.tsv".
        The files will be stored in a subdirectory of ``self.dir_output``.

        :param dir_oligosets: Subdirectory to store oligosets, defaults to "oligosets"
        :type dir_oligosets: str, optional
        :return: Path to oligosets directory.
        :rtype: str
        """
        dir_oligosets = os.path.join(self.dir_output, dir_oligosets)
        Path(dir_oligosets).mkdir(parents=True, exist_ok=True)

        for region_id in self.oligosets.keys():
            file_oligosets = os.path.join(dir_oligosets, f"{region_id}_oligosets.tsv")
            self.oligosets[region_id].to_csv(file_oligosets, sep="\t", index=False)

        return dir_oligosets

    def remove_regions_with_insufficient_oligos(
        self,
        pipeline_step: str,
        write: bool = True,
    ):
        """Deletes from the ``oligo_DB`` the regions (e.g. genes) which have less than ``min_oligos_per_region`` oligos,
        and optionally writes them in a file with the name of the step of the pipeline at which they have been deleted.

        :param pipeline_step: Step in the pipeline that lead to the removal of the region.
        :type pipeline_step: str
        :param write: write removed regions to file ``regions_with_insufficient_oligos.txt``, defaults to True
        :type write: bool, optional
        """
        regions = copy(list(self.database.keys()))
        for region in regions:
            if len(list(self.database[region].keys())) <= self.min_oligos_per_region:
                del self.database[region]
                if region in self.oligosets:
                    del self.oligosets[region]
                if write:
                    with open(self.file_removed_regions, "a") as hanlde:
                        hanlde.write(f"{region}\t{pipeline_step}\n")

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
        database_entries = {}  # assiciate each id to the additional features
        oligo_sequence_ids = {}  # associate each sequence to its id

        # parse the exon fasta sequence file
        for region_sequence in region_sequences:
            seq = region_sequence.seq
            (
                region,
                additional_information,
                coordinates,
            ) = parse_fasta_header(region_sequence.id)

            # should all have the same chromosome and strand entry
            oligo_chrom = coordinates["chromosome"][0]
            oligo_strand = coordinates["strand"][0]

            list_of_coordinates = []
            # coordinates in fasta file use 0-based indixing, which go from 0 (for base 1) to n (for base n)
            # to calculate the start and end coordinates for oligos spanning split sequences (e.g. exon junctions)
            # we need to remove the additional 0-based index, i.e. ranging from start + 1 to end
            # e.g. for the first 10 nucleotides, the bed annotation is 0-10, but we need numbers from 1-10
            # once we have the correct start and end corrdinates we turn them again into 0-based coordinates
            # by subtracting 1. This needs to be done afterwards because on the minus strand the start has a
            # higher number than the end, which needs to be sorted and turned into 0-based index
            for i in range(len(coordinates["start"])):
                list_of_coordinates.extend(
                    list(
                        range(
                            coordinates["start"][i] + 1,
                            coordinates["end"][i] + 1,
                        )
                    )
                )
            # sort reverse on minus strand caus the sequence is translated into
            # the reverse complement by fasta -strand option
            if oligo_strand == "-":
                list_of_coordinates.reverse()

            for oligo_length in range(self.oligo_length_min, self.oligo_length_max + 1):
                if len(seq) > oligo_length:
                    number_oligos = len(seq) - (oligo_length - 1)
                    oligos = [seq[i : i + oligo_length] for i in range(number_oligos)]

                    for i in range(number_oligos):
                        oligo = oligos[i]
                        oligo_start_end = [
                            list_of_coordinates[i],
                            list_of_coordinates[(i + oligo_length - 1)],
                        ]
                        oligo_start_end.sort()
                        oligo_start = oligo_start_end[0] - 1  # turn into 0-based index
                        oligo_end = oligo_start_end[1]
                        oligo_id = f"{region_id}_{oligo_chrom}:{oligo_start}-{oligo_end}({oligo_strand})"
                        oligo_attributes = {
                            "sequence": oligo,
                            "chromosome": oligo_chrom,
                            "start": oligo_start,
                            "end": oligo_end,
                            "strand": oligo_strand,
                            "length": oligo_length,
                            "additional_information_fasta": [additional_information],
                        }

                        if oligo in oligo_sequence_ids:
                            if oligo_id in oligo_sequence_ids[oligo].keys():
                                entry_additional_information = oligo_sequence_ids[
                                    oligo
                                ][oligo_id]["additional_information_fasta"]
                                entry_additional_information.append(
                                    additional_information
                                )
                                # remove exon junctions entry if exon entry exists
                                oligo_sequence_ids[oligo][oligo_id][
                                    "additional_information_fasta"
                                ] = [
                                    info
                                    for info in entry_additional_information
                                    if "__JUNC__" not in info
                                ]
                            else:
                                oligo_sequence_ids[oligo][oligo_id] = oligo_attributes
                        else:
                            oligo_sequence_ids[oligo] = {oligo_id: oligo_attributes}

        database_entries = {}
        for oligo_ids in oligo_sequence_ids.values():
            oligo_id = ";".join([oligo_id for oligo_id in oligo_ids.keys()])
            attributes_collapsed = {}
            for key in [
                "sequence",
                "chromosome",
                "start",
                "end",
                "strand",
                "length",
            ]:
                attributes_collapsed[key] = list(
                    set([attributes[key] for attributes in oligo_ids.values()])
                )
                if len(attributes_collapsed[key]) == 1 and key not in ["start", "end"]:
                    attributes_collapsed[key] = attributes_collapsed[key][0]
            # can't be collapsed caus sometimes contains lists of lists
            attributes_collapsed["additional_information_fasta"] = sum(
                [
                    attributes["additional_information_fasta"]
                    for attributes in oligo_ids.values()
                ],
                [],
            )
            database_entries[oligo_id] = attributes_collapsed

        database = {region_id: database_entries}

        return database
