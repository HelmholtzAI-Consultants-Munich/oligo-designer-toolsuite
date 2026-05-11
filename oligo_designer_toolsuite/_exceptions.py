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
