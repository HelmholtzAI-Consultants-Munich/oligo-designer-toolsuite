from abc import ABC, abstractmethod

from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import Seq


class PropertyFilterBase(ABC):
    """Base class that gives the structure for Oligo Property Filters."""

    def __init__(self) -> None:
        pass

    @abstractmethod
    def apply(self, sequence: str):
        """
        Applies the filters to a given sequence and if the sequence fulfillts the constraints returns ``True`` and dictionary that stores additional computed features.
        If this method is not reimplemented in the filters classes it will give a warning.
        The aditional computed features must be float type.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: constraints outcome, dictionary with the additional features
        :rtype: bool, dict
        """


class MaskedSequences(PropertyFilterBase):
    """Checks if the sequences contains a masked nucleotide."""

    def __init__(self, 
            mask: str = 'N'
        ) -> None:
        """Constructor"""
        super().__init__()
        self.mask = mask

    def apply(self, sequence: Seq):
        """Applies the filter and returns True if there isn't any masked nucleotide.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constraint is fulfilled
        :rtype: bool
        """
        if self.mask in sequence:
            return False, {}
        return True, {}


class GCContent(PropertyFilterBase):
    """Checks if the GC content of a sequence lies within a user defined interval [GC_content_min, GC_content_max].

    :param GC_content_min: minumum GC content value that the oligos need to have
    :type GC_content_min: float
    :param GC_content_max: maximum GC content value that the oligos need to have
    :type GC_content_max: float
    """

    def __init__(self, 
            GC_content_min: float, 
            GC_content_max: float
        ) -> None:
        """Constructor"""
        super().__init__()
        self.GC_content_min = GC_content_min
        self.GC_content_max = GC_content_max
        assert GC_content_max >= GC_content_min, "Max value is lower that min value!"

    def apply(self, sequence: Seq):
        """Applies the filter and returns True if the GC content fulfills requirements.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return:True if the constraint is fulfilled, GC content
        :rtype: bool, dict
        """
        GC_content = round(gc_fraction(sequence)*100, 4)
        sequence_features = {"GC_content": GC_content}
        if self.GC_content_min < GC_content < self.GC_content_max:
            return True, sequence_features
        return False, {}  # if false the additional features are not been saved anyway


class MeltingTemperatureNN(PropertyFilterBase):
    """Checks if the melting temperature of a sequence lies within a user defined interval [Tm_min, Tm_max].
    The melting tenperature is computed using nearest neighbor thermodynamics.
    The parameters for melting temperature computation can be changed from default by providing ```Tm_parameters``` 
    dict with parameters specifications.
    The melting temperature can be corrected for salt ions by providing a ```Tm_salt_correction_parameters```dict 
    with parameters specifications.
    The melting temperature can be corrected for DMSO and formamide by providing a ```Tm_chem_correction_parameters```dict 
    with parameters specifications.

    :param Tm_min: minimum melting temperature
    :type Tm_min: float
    :param Tm_max: maximum melting temperature
    :type Tm_max: float
    :param Tm_parameters: parameters to compute the melting temperature, 
        set to ```{}``` (empty dict) if you wish to use Bio.SeqUtils.MeltingTemp default parameters
        for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_salt_correction_parameters: parameters to correct the melting temperature for salt ions, defaults to None, i.e. no correction
        set to ```{}``` (empty dict) if you wish to use Bio.SeqUtils.MeltingTemp default parameters
        for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
    :type Tm_salt_correction_parameters: dict
    :param Tm_chem_correction_parameters: parameters to correct the melting temperature for DMSO and formamide, defaults to None, i.e. no correction
        set to ```{}``` (empty dict) if you wish to use Bio.SeqUtils.MeltingTemp default parameters
        for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
    :type Tm_chem_correction_parameters: dict
    """
    def __init__(self,
            Tm_min: float,
            Tm_max: float,
            Tm_parameters: dict,
            Tm_salt_correction_parameters: dict = None,
            Tm_chem_correction_parameters: dict = None,
        ) -> None:
        """Constructor"""
        super().__init__()
        self.Tm_min = Tm_min
        self.Tm_max = Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters

    def __get_Tm(self, sequence: Seq):
        """Computes the melting temperature of the sequence and applies salt or chemical corrections if parameters provided.

        :param sequence: sequence for which the melting temperature is computed
        :type sequence: str
        :return: melting temperature
        :rtype: float
        """
        Tm = mt.Tm_NN(sequence, **self.Tm_parameters)
        if self.Tm_salt_correction_parameters is not None:
            Tm = mt.salt_correction(Tm, **self.Tm_salt_correction_parameters)
        if self.Tm_chem_correction_parameters is not None:
            Tm = mt.chem_correction(Tm, **self.Tm_chem_correction_parameters)
        return round(Tm, 4)

    def apply(self, sequence: Seq):
        """Applies the filter and returns True if the melting temperature fulfills requirements.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constraint is fulfilled, melting temperature
        :rtype: bool, dict
        """
        Tm = self.__get_Tm(sequence)
        sequence_features = {"melting_temperature": Tm}
        if self.Tm_min < Tm < self.Tm_max:
            return True, sequence_features
        return False, {Tm}


class ProhibitedSequences(PropertyFilterBase):
    """Filters the sequences containing a prohibited sequence."""

    def __init__(
        self,
    ) -> None:
        """Initializes the class."""

        super().__init__()

    def apply(self, sequence: Seq):
        """Applies the filter and returns True if the oligo does not contain prohibited sequences.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrined is fulfilled and the melting temperature
        :rtype: bool and dict
        """

class GCClamp(PropertyFilterBase):
    """ Filters the sequences by the presence of a GC Clamp: one of the n 3' terminal bases must be G or C
    :param n: number of terminal bases to check for a G or a C
    :type  n: int


    """
    def __init__(self, n:int) -> None:
        super().__init__()
        self.n = n

    def apply(self, sequence: Seq):
        """Applies the filter and returns True if there is a GC clamp

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrain is fulfilled
        :rtype: bool
        """
        for i in range (self.n):
            if (sequence[-i-1]=='G' or sequence[-i-1]=='C'):
                return True, {}
        return False, {}

class ConsecutiveRepeats(PropertyFilterBase):
    """Filters the sequences by consecutive repeats.

    :param Repeat_num_max: maximum number of consecutive repeats that the oligos can not have
    :type Repeat_num_max: int
    """

    def __init__(self, repeat_num_max: int) -> None:
        """Constructor"""
        super().__init__()

        self.repeat_num_max = repeat_num_max

    def apply(self, sequence: Seq):
        """Applies the filter and returns True if the number of consecutive repeats are lower the max values.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrined is fulfilled, the max consecutive repeats
        :rtype: bool and dict
        """
        repeat_num_max = self.repeat_num_max

        if repeat_num_max <= 0:
            self.repeat_num_max = 1
            repeat_num_max = self.repeat_num_max

        elif repeat_num_max >= len(sequence):
            repeat_num_max = len(sequence) - 1
            # max_num_repeat is changed to length of sequence - 1

        AA_consecutive_repeat_min = min([len(set(kmer)) for kmer in
                                         (sequence[i:i + repeat_num_max + 1] for i in
                                          range(len(sequence) - repeat_num_max))])

        sequence_features = {"Consecutive_repeat_min": AA_consecutive_repeat_min}

        if (AA_consecutive_repeat_min != 1) | (repeat_num_max == 1):
            return True, sequence_features
        return False, {}  # if false the additional features are not been saved anyway
