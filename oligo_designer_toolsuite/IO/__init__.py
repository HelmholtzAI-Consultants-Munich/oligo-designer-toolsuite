"""
This module manages the databases required for the oligos generations. In particular, it deals with the input and output of data.

:class CustomDB: Class that uses previously dowloaded .gtf and .fasta files to generate a Database containing all the oligos and magaes it.
:class NcbiDB: Class that dowloads a specific vesrion of the .gtf and .fasta files from the NCBI server to generate a Database containing all the oligos and magaes it.
:class EnsemblDB: Class that dowloads a specific vesrion of the .gtf and .fasta files from the Ensembl server to generate a Database containing all the oligos and magaes it.
"""
