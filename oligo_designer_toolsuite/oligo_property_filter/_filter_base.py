############################################
# imports
############################################

from abc import ABC, abstractmethod

############################################
# Property Filter Class
############################################


class PropertyFilterBase(ABC):

    def __init__(self) -> None:
        """Constructor for the PropertyFilterBase class."""

    @abstractmethod
    def apply(self, sequence: str) -> bool:
        """ """
