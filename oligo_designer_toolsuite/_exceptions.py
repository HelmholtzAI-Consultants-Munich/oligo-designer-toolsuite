############################################
# imports
############################################

############################################
# Custom Exception Classes
############################################


class OligoDesignerError(Exception):
    """
    Base exception class for all Oligo Designer Toolsuite errors.

    All custom exceptions in this package should inherit from this class.
    This allows users to catch all package-specific errors with a single exception type.
    """


class ConfigurationError(OligoDesignerError):
    """
    Raised when there is an error in configuration files or parameters.

    This exception is used for issues with YAML configuration files, missing required
    parameters, invalid parameter values, or configuration file format errors.
    """


class DatabaseError(OligoDesignerError):
    """
    Raised when there is an error related to database operations.

    This exception is used for issues with OligoDatabase or ReferenceDatabase operations,
    such as empty databases, missing regions/oligos, or database format errors.
    """


class FileFormatError(OligoDesignerError):
    """
    Raised when a file format is incorrect or unsupported.

    This exception is used when files (FASTA, VCF, GFF, GTF, etc.) are malformed,
    missing required fields, or in an unsupported format.
    """


class FeatureNotImplementedError(OligoDesignerError):
    """
    Raised when a feature or functionality is not yet implemented.

    This exception is used when a method or feature is planned but not yet
    implemented in the codebase.
    """


class ExternalToolError(OligoDesignerError):
    """
    Raised when an external command line tool such as BLAST, Bowtie, bcftools or bedtools fails.
    """


class NetworkError(OligoDesignerError):
    """
    Raised when a remote server such as the NCBI or Ensembl FTP server cannot be reached.
    """


class EmptyResultError(OligoDesignerError, SystemExit):
    """
    Raised when no oligos are left to continue with.

    Inherits from SystemExit so an uncaught error still ends a command line run with
    exit code 1 and no traceback, like the ``sys.exit(1)`` it replaces.
    """

    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.code = 1
