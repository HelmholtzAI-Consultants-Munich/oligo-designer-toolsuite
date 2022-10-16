from abc import ABC, abstractmethod

from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt


class PreFilterBase(ABC):
    """Base class that gives the structure."""

    def __init__(self) -> None:
        pass

    @abstractmethod
    def apply(self, sequence):
        """Applies teh filters and returns if the sequence fulfillts the constraints and the additional features computed in a dictionary.
        If this method is not reimplemented in teh filters gives a warning and returns True.
        The aditional computed features must be float type.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True
        :rtype: bool
        """


class MaskedSequences(PreFilterBase):
    """Filters the sequences containing a masked nucleotide."""

    def __init__(self) -> None:
        super().__init__()

    def apply(self, sequence):
        """Applies the filter and returns True if there isn't any masked nucleotide.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrined is fulfilled
        :rtype: bool
        """

        if "N" in sequence:
            return False, {}
        return True, {}


class GCContent(PreFilterBase):
    """Filters the sequences by the GC content."""

    def __init__(self, GC_content_min, GC_content_max) -> None:
        super().__init__()
        self.GC_content_min = GC_content_min
        self.GC_content_max = GC_content_max

    def apply(self, sequence):
        """Applies the filter and returns True if the GC content is between the min and max values and the GC content.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrined is fulfilled and the GC content
        :rtype: bool and dict
        """

        GC_content = round(GC(sequence), 2)
        sequence_features = {"GC_content": GC_content}
        if self.GC_content_min < GC_content < self.GC_content_max:
            return True, sequence_features
        return False, {}  # if false the additional features are not been saved anyway


class MeltingTemperature(PreFilterBase):
    """Filters the sequences by the melting temperature."""

    def __init__(self, Tm_min, Tm_max, Tm_parameters, Tm_correction_parameters) -> None:
        """Initializes the class.

        :param Tm_min: minimum melting temperature
        :type Tm_min: float
        :param Tm_max: maximum melting temperature
        :type Tm_max: float
        :param Tm_parameters: parameters to compute the melting temperature
        :type Tm_parameters: dict
        :param Tm_correction_parameters: parameters to correct the melting temperature
        :type Tm_correction_parameters: dict
        """

        super().__init__()
        self.Tm_min = Tm_min
        self.Tm_max = Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_correction_parameters = Tm_correction_parameters

    def __get_Tm(self, sequence):
        """Computes the melting temperature of the sequence.

        :param sequence: sequence for which the melting temperature is computed
        :type sequence: str
        :return: melting temperature
        :rtype: float
        """

        Tm = mt.Tm_NN(sequence, **self.Tm_parameters)
        Tm_corrected = round(mt.chem_correction(Tm, **self.Tm_correction_parameters), 2)
        return Tm_corrected

    def apply(self, sequence):
        """Applies the filter and returns True if the melting temperature is between the min and max values and the melting temperature.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrined is fulfilled and the melting temperature
        :rtype: bool and dict
        """

        Tm = self.__get_Tm(sequence)
        sequence_features = {"melt_temp": Tm}
        if self.Tm_min < Tm < self.Tm_max:
            return True, sequence_features
        return False, {}


"""

def _filter_by_kmer():
    pass


def _mask_prohibited_sequences():
    pass
"""
