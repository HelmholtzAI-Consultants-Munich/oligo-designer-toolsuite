############################################
# imports
############################################

from abc import ABC, abstractmethod

############################################
# Property Filter Class
############################################


class PropertyFilterBase(ABC):
    """
    An abstract base class for property filters in oligo design pipelines.

    This class serves as a template for creating custom filters that can be applied to oligo sequences.
    """

    def __init__(self) -> None:
        """Constructor for the PropertyFilterBase class."""

    @abstractmethod
    def apply(self, sequence: str):
        """
        Abstract method to apply the filter to a given sequence.
        If the sequence fulfillts the constraints the function returns ``True`` and a dictionary that stores additional computed features.
        Note: a warning is thrown if this class is not reimplemented in the custom filter class.

        :param sequence: The oligo sequence to be filtered.
        :type sequence: str
        """
