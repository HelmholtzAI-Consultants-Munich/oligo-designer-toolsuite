############################################
# imports
############################################

from typing import Literal

############################################
# types
############################################

_TYPES_REF = Literal["fasta", "vcf"]
_TYPES_FILE = Literal["gff", "gtf", "fasta"]
_TYPES_FILE_SEQ = Literal["dna", "ncrna"]

############################################
# constants
############################################

SEPARATOR_OLIGO_ID = "::"
SEPARATOR_FASTA_HEADER_FIELDS = "::"
SEPARATOR_FASTA_HEADER_FIELDS_LIST = ";"

SUPPORTED_TAXA_SOURCES: dict[str, set[str]] = {
    "archaea": {"latest_assembly_versions", "reference"},
    "bacteria": {"latest_assembly_versions", "reference"},
    "fungi": {"latest_assembly_versions", "reference"},
    "invertebrate": {"annotation_releases", "latest_assembly_versions", "reference"},
    "metagenomes": {"latest_assembly_versions"},
    "plant": {"annotation_releases", "latest_assembly_versions", "reference"},
    "plants": {"annotation_releases", "latest_assembly_versions", "reference"},
    "protozoa": {"latest_assembly_versions", "reference"},
    "unknown": {"latest_assembly_versions"},
    "vertebrate_mammalian": {"annotation_releases", "latest_assembly_versions", "reference"},
    "vertebrate_other": {"annotation_releases", "latest_assembly_versions", "reference"},
    "viral": {"latest_assembly_versions"},
}
