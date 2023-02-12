############################################
# imports
############################################

import os
import warnings
from Bio import SeqIO
from pathlib import Path

from ..utils._data_parser import check_fasta_format, parse_fasta_header

############################################
# Reference Database Class
############################################


class ReferenceDatabase:
    """This class stores a refernce database, which is created from
    a sequence fasta file. The header of each sequence must start with '>' and contain the
    following information: region_id, additional_information (optional) and coordinates
    (chrom, start, end, strand).

    Input Format (per sequence):
    >region_id::additional information::chromosome:start-end(strand)
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    :param file_fasta: Path to the fasta file.
    :type file_fasta: str
    :param source: Source of annotations, e.g. NCBI, defaults to None.
    :type source: str, optional
    :param species: Species of annotation, e.g. Homo_sapiens, defaults to None.
    :type species: str, optional
    :param annotation_release: Release number of annotation, e.g. 110, defaults to None.
    :type annotation_release: str, optional
    :param genome_assembly: Genome assembly of annotation, e.g. GRCh38, defaults to None.
    :type genome_assembly: str, optional
    :param dir_output: Output directory, defaults to 'output'.
    :type dir_output: str, optional
    """

    def __init__(
        self,
        file_fasta: str,
        source: str = None,
        species: str = None,
        annotation_release: str = None,
        genome_assembly: str = None,
        dir_output: str = "output",
    ):
        """Constructor"""
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

        self.file_fasta = file_fasta
        if os.path.exists(self.file_fasta):
            if not check_fasta_format(self.file_fasta):
                raise ValueError("Fasta file has incorrect format!")
        else:
            raise ValueError("Fasta file does not exist!")

        self.database = []

        self.dir_output = os.path.abspath(
            os.path.join(dir_output, "reference_database")
        )
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

    def load_fasta_into_database(self):
        """Load sequences from fasta file and stored in databse."""
        with open(self.file_fasta, "r") as handle:
            self.database = list(SeqIO.parse(handle, "fasta"))

    def write_fasta_from_database(
        self,
        filename: str,
    ):
        """Write sequences stored in database to fasta file.

        :param filename: Fasta file name prefix.
        :type filename: str
        :return: Path to fasta file.
        :rtype: str
        """
        file_database = os.path.join(self.dir_output, f"{filename}.fna")
        if self.database:
            with open(file_database, "w") as handle_fasta:
                SeqIO.write(self.database, handle_fasta, "fasta")
        else:
            raise ValueError("Database is empty! Nothing written to fasta file.")

        return file_database

    def filter_database(
        self,
        region_ids: list[str] = None,
    ):
        """Filter database based on region id and overwrite database with filtered database.

        :param region_ids: List of region ids.
        :type region_ids: list
        """
        if self.database:
            database_filtered = []
            for entry in self.database:
                region, _, _ = parse_fasta_header(entry.id)
                if region in region_ids:
                    database_filtered.append(entry)
            self.database = database_filtered
        else:
            raise ValueError("Can not filter. Database is empty!")
