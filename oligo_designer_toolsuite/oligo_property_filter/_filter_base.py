from abc import ABC, abstractmethod

from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt


class PreFilterBase(ABC):
    """Base class that gives the structure."""

    def __init__(self) -> None:
        pass

    @abstractmethod
    def apply(self, sequence):
        """
        Applies the filters to a given sequence and if it fulfillts the constraints returns ``True`` and the additional features computed stored in a dictionary.
        If this method is not reimplemented in the filters classes it will give a warning.
        The aditional computed features must be float type.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: constraints outcome, dictionary with the additional features
        :rtype: bool, dict
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
    """Filters the sequences by the GC content.

    :param GC_content_min: minumum GC content value that the oligos need to have
    :type GC_content_min: float
    :param GC_content_max: maximum GC content value that the oligos need to have
    :type GC_content_max: float
    """

    def __init__(self, GC_content_min, GC_content_max) -> None:
        """Constructor"""
        super().__init__()
        self.GC_content_min = GC_content_min
        self.GC_content_max = GC_content_max

    def apply(self, sequence):
        """Applies the filter and returns True if the GC content is between the min and max values.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrined is fulfilled, the GC content
        :rtype: bool and dict
        """

        GC_content = round(GC(sequence), 2)
        sequence_features = {"GC_content": GC_content}
        if self.GC_content_min < GC_content < self.GC_content_max:
            return True, sequence_features
        return False, {}  # if false the additional features are not been saved anyway


class MeltingTemperature(PreFilterBase):
    """Filters the sequences by the melting temperature.

    :param Tm_min: minimum melting temperature
    :type Tm_min: float
    :param Tm_max: maximum melting temperature
    :type Tm_max: float
    :param Tm_parameters: parameters to compute the melting temperature, for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_correction_parameters: parameters to correct the melting temperature,for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_correction_parameters: dict
    """

    def __init__(self, Tm_min, Tm_max, Tm_parameters, Tm_correction_parameters) -> None:
        """Initializes the class."""

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
