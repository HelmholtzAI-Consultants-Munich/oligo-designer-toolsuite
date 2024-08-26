############################################
# imports
############################################

from typing import Literal

############################################
# types
############################################

_TYPES_SEQ = Literal["target", "oligo", "sequence_encoding_probe"]
_TYPES_FILE = Literal["gff", "gtf", "fasta"]
_TYPES_FILE_SEQ = Literal["dna", "ncrna"]

############################################
# constants
############################################

SEPARATOR_OLIGO_ID = "::"
SEPARATOR_FASTA_HEADER_FIELDS = "::"
SEPARATOR_FASTA_HEADER_FIELDS_LIST = ";"
SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS = ","
