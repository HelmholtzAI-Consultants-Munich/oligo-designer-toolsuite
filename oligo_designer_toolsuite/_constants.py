############################################
# imports
############################################

from typing import Literal

############################################
# types
############################################

_TYPES_SEQ = Literal["target", "target_short", "oligo", "oligo_short", "sequence_encoding_probe"]
_TYPES_REF = Literal["fasta", "vcf"]
_TYPES_FILE = Literal["gff", "gtf", "fasta"]
_TYPES_FILE_SEQ = Literal["dna", "ncrna"]

############################################
# constants
############################################

SEPARATOR_OLIGO_ID = "::"
SEPARATOR_FASTA_HEADER_FIELDS = "::"
SEPARATOR_FASTA_HEADER_FIELDS_LIST = ";"
