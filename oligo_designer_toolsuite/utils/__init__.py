"""
This module provides utilities for processing genomic databases, parsing sequences,
and performing various sequence-related tasks.

Classes:
- FastaParser: Parses FASTA files and provides methods for extracting sequence information.
- GffParser: Parses GFF/GTF files and provides methods for extracting annotation information.
"""

from ._sequence_parser import FastaParser, GffParser


__all__ = [
    "FastaParser",
    "GffParser",
]
