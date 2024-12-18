############################################
# imports
############################################

from abc import ABC, abstractmethod

############################################
# Property Filter Class
############################################


class PropertyFilterBase(ABC):
    """
    An abstract base class for creating sequence property filters.

    The `PropertyFilterBase` class serves as a template for developing filters that evaluate sequences based on specific criteria.
    Subclasses must implement the `apply` method to define the filtering logic.

    :param sequence: The sequence to be evaluated by the filter.
    :type sequence: str
    """

    def __init__(self) -> None:
        """Constructor for the PropertyFilterBase class."""

    @abstractmethod
    def apply(self, sequence: str) -> bool:
        """
        Evaluate whether the sequence meets the filter's criteria.

        This abstract method must be implemented by subclasses to define the specific filtering logic for a given sequence.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence meets the filter's criteria, `False` otherwise.
        :rtype: bool
        """
