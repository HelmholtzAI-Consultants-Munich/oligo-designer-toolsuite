############################################
# imports
############################################
from typing import Union, List

from Bio.SeqUtils import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt

from seqfold import dg

from . import PropertyFilterBase


###TODO move this function to utils once database refactor is merged
def check_sequence(seq: str, valid_characters={"A", "C", "T", "G"}) -> bool:
    return all(char in valid_characters for char in seq)


############################################
# Oligo Property Filter Classes
############################################


class SoftMaskedSequenceFilter(PropertyFilterBase):
    """A filter to check if a DNA sequence is soft-masked (contains lowercase letters)."""

    def __init__(self) -> None:
        """Constructor for the SoftMaskedSequenceFilter class."""
        super().__init__()

    def apply(self, sequence: Seq) -> (bool, dict):
        """
        Applies the soft mask filter to a DNA sequence and returns True if the sequence does not contain lower-case letters.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and an empty dictionary.
        :rtype: (bool, dict)
        """
        if any(c.islower() for c in sequence):
            return False, {}
        return True, {}


class HardMaskedSequenceFilter(PropertyFilterBase):
    """A filter to check if a DNA sequence contains a specific mask character (default is "N").

    :param mask: The mask character to check for in the sequence. Default is "N".
    :type mask: str
    """

    def __init__(self, mask: str = "N") -> None:
        """Constructor for the HardMaskedSequenceFilter class."""
        super().__init__()
        self.mask = mask

    def apply(self, sequence: Seq):
        """
        Applies the hard mask filter to a DNA sequence and returns True if the sequence does not contain the mask letter.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and an empty dictionary.
        :rtype: (bool, dict)
        """
        if self.mask in sequence:
            return False, {}
        return True, {}


class ProhibitedSequenceFilter(PropertyFilterBase):
    """A filter that checks for the absence of a specified prohibited sequence in a DNA sequence.
    The filter will reject any sequence that contains the prohibited sequence.

    :param prohibited_sequence: The DNA sequence that is to be prohibited. This sequence should
                                be a valid DNA sequence consisting of characters A, T, C, and G only.
                                The user can privide either a single sequence or a list of sequences.
    :type prohibited_sequence: str, list[str]
    """

    def __init__(
        self,
        prohibited_sequence: Union[str, List[str]],
    ) -> None:
        """Constructor for the ProhibitedSequenceFilter class."""
        super().__init__()
        if not check_sequence(prohibited_sequence):
            raise ValueError("Prohibited sequence ({prohibited_sequences}) is not a DNA sequence.")
        if not isinstance(prohibited_sequence, list):
            prohibited_sequence = [prohibited_sequence]
        self.prohibited_sequence = [s.upper() for s in prohibited_sequence]

    def apply(self, sequence: Seq):
        """
        Applies the filter to a given DNA sequence to check if it contains the prohibited sequence.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and an empty dictionary.
        :rtype: (bool, dict)
        """
        for s in self.prohibited_sequence:
            if s in sequence.upper():
                return False, {}
        return True, {}


class HomopolymericRunsFilter(PropertyFilterBase):
    """A filter hat checks for the absence of a specified homopolymeric run in a DNA sequence.
    A homopolymeric run is defined as a sequence where the same nucleotide base repeats consecutively.

    :param base: The nucleotide base (A, T, C, or G) that is checked for consecutive repeats. The parameter
                 can be a single base or a list of bases. If a list is provided, the filter will check for consecutive
                 repeats for each base in the list.
    :type base: str, List[str]
    :param n: The minimum number of consecutive repeats of the base that defines a homopolymeric run. If ``base`` is a list,
              this parameter can be a single integer or a list of integers of equal lenght. If a single integer is provided, it will be assigned 
              to all the bases in the list. Alternatively, if a list of integers is provided each element will be assigned to the bases by matching indices.
    :type n: int, List[int]
    """

    def __init__(
        self,
        base: Union[str, List[str]],
        n: Union[int, List[int]],
    ) -> None:
        """Constructor for the HomopolymericRunsFilter class."""
        super().__init__()
        if not check_sequence(base):
            raise ValueError("Prohibited sequence ({base}) is not a DNA sequence.")
        # Check that the variables types are comaptible
        if not isinstance(base, list):
            if isinstance(n, list):
                raise TypeError("The variable n cannot be type list when base is type string.")
            else:
                base = [base]
                n = [n]
        elif isinstance(base, list):
            if isinstance(n, list) and len(base) != len(n):
                raise ValueError(f"The lists base and n must have the same length, but they have {len(base)} and {len(n)} repectively.")
            elif not isinstance(n, list):
                # n is the same for all the elements of base
                n = [n for _ in range(len(base))]
            elif not isinstance(n, list) and not isinstance(n, int):
                raise TypeError("The variable n is expected to be an integer or a list of integers.")
        else:
            raise TypeError("The variable base is expected to be a string or a list of strings.")
        
        # base and n are now lists of the same length
        self.base = [nucleotide.upper() for nucleotide in base]
        self.n = n
        self.homopolymeric_runs = [nucleotide * repepeats for nucleotide, repepeats in zip(self.base, self.n)]

    def apply(self, sequence: Seq):
        """Applies the filter to a given DNA sequence to check if it contains a homopolymeric run.

        :param sequence: The DNA sequence to be checked for homopolymeric runs.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and an empty dictionary.
        :rtype: (bool, dict)
        """
    
        for homopolymeric_run in self.homopolymeric_runs:
            if homopolymeric_run in sequence.upper():
                return False, {}
        return True, {}


class FivePrimeSequenceFilter(PropertyFilterBase):
    """A filter to check the presence or absence of a specified sequence at the 5'-end of a DNA sequence.

    :param five_prime_sequence: The sequence to check at the 5'-end of the DNA sequence.
    :type five_prime_sequence: str
    :param remove: If True, sequences starting with the specified sequence are filtered out. If False, only sequences starting with the specified sequence are retained.
    :type remove: bool
    """

    def __init__(self, five_prime_sequence: str, remove: bool = True) -> None:
        """Constructor for the FivePrimeSequenceFilter class."""
        super().__init__()
        self.five_prime_sequence = five_prime_sequence.upper()
        self.remove = remove

    def apply(self, sequence: str):
        """Applies the 5'-end sequence filter to a DNA sequence and eitehr keeps or removes the matching sequence, dependend on the parameter "remove".

        :param sequence: The DNA sequence to be checked.
        :type sequence: str
        :return: A tuple indicating if the sequence passes the filter and an empty dictionary.
        :rtype: (bool, dict)
        """
        if self.remove:
            if sequence.upper().startswith(self.five_prime_sequence):
                return False, {}
            return True, {}
        else:
            if sequence.upper().startswith(self.five_prime_sequence):
                return True, {}
            return False, {}


class ThreePrimeSequenceFilter(PropertyFilterBase):
    """A filter to check the presence or absence of a specified sequence at the 3'-end of a DNA sequence.

    :param three_prime_sequence: The sequence to check at the 3'-end of the DNA sequence.
    :type three_prime_sequence: str
    :param remove: If True, sequences ending with the specified sequence are filtered out. If False, only sequences ending with the specified sequence are retained.
    :type remove: bool
    """

    def __init__(self, three_prime_sequence: str, remove: bool = True) -> None:
        """Constructor for the ThreePrimeSequenceFilter class."""
        super().__init__()
        self.three_prime_sequence = three_prime_sequence.upper()
        self.remove = remove

    def apply(self, sequence: str):
        """Applies the 3'-end sequence filter to a DNA sequence and eitehr keeps or removes the matching sequence, dependend on the parameter "remove".

        :param sequence: The DNA sequence to be checked.
        :type sequence: str
        :return: A tuple indicating if the sequence passes the filter and an empty dictionary.
        :rtype: (bool, dict)
        """
        if self.remove:
            if sequence.upper().endswith(self.three_prime_sequence):
                return False, {}
            return True, {}
        else:
            if sequence.upper().endswith(self.three_prime_sequence):
                return True, {}
            return False, {}


class GCContentFilter(PropertyFilterBase):
    """A filter to check if the GC content of a DNA sequence falls within a specified range [GC_content_min, GC_content_max].

    :param GC_content_min: The minimum acceptable GC content as a percentage.
    :type GC_content_min: float
    :param GC_content_max: The maximum acceptable GC content as a percentage.
    :type GC_content_max: float
    """

    def __init__(self, GC_content_min: float, GC_content_max: float) -> None:
        """Constructor for the GCContentFilter class."""
        super().__init__()
        if GC_content_max <= GC_content_min:
            raise ValueError("GC_content_max is lower that GC_content_min!")
        self.GC_content_min = GC_content_min
        self.GC_content_max = GC_content_max

    def apply(self, sequence: Seq):
        """Applies the GC content filter to a DNA sequence.

        :param sequence: The DNA sequence to be checked for GC content.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and a dictionary containing
                 the actual GC content if the condition is met.
        :rtype: (bool, dict)
        """
        GC_content = round(gc_fraction(sequence) * 100, 4)
        if self.GC_content_min < GC_content < self.GC_content_max:
            return True, {"GC_content": GC_content}
        return False, {}


class GCClampFilter(PropertyFilterBase):
    """A filter to check if the n 3' terminal bases end of a DNA sequence contain at least one G or C bases.

    :param n: The number of bases from the 3' end of the sequence to check for the presence of G or C.
    :type n: int
    """

    def __init__(self, n: int) -> None:
        """Constructor for the GCClampFilter class."""
        super().__init__()
        self.n = n

    def apply(self, sequence: Seq):
        """Applies the GC clamp filter to the 3' end of a DNA sequence and returns True if there is a GC clamp.
        A GC clamp is the presence of a guanine (G) or cytosine (C) base in the last n bases (the 3' end) of an oligo.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and an empty dictionary.
        :rtype: (bool, dict)
        """
        for i in range(1, self.n + 1):
            if sequence.upper()[-i] == "G" or sequence.upper()[-i] == "C":
                return True, {}
        return False, {}


class MeltingTemperatureNNFilter(PropertyFilterBase):
    """A filter to determine if the melting temperature (Tm) of a DNA sequence falls within a specified range [Tm_min, Tm_max].
    It uses nearest-neighbor thermodynamic models to calculate the Tm, with optional salt and chemical corrections.

    The parameters for Tm calculation can be adjusted from default by providing the ``Tm_parameters`` dict with parameters specifications.
    The Tm can be corrected for salt ions by providing the ``Tm_salt_correction_parameters`` dict with parameters specifications.
    The Tm can be corrected for DMSO and formamide by providing a ``Tm_chem_correction_parameters`` dict with parameters specifications.

    :param Tm_min: The minimum acceptable melting temperature.
    :type Tm_min: float
    :param Tm_max: The maximum acceptable melting temperature.
    :type Tm_max: float
    :param Tm_parameters: Parameters for the nearest-neighbor thermodynamic model to calculate Tm.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_salt_correction_parameters: Optional parameters for salt correction of Tm calculations.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Optional parameters for chemical correction.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
    :type Tm_chem_correction_parameters: dict, optional
    """

    def __init__(
        self,
        Tm_min: float,
        Tm_max: float,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ) -> None:
        """Constructor for the MeltingTemperatureNNFilter class."""
        super().__init__()
        if Tm_max <= Tm_min:
            raise ValueError("Tm_max is lower that Tm_min!")
        self.Tm_min = Tm_min
        self.Tm_max = Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters

    ###TODO: move this function to utils as it is also used in the padlock arm filter once database refactor is merged
    def _get_Tm(self, sequence: Seq):
        """Internal method to calculate the melting temperature of a sequence.

        :param sequence: The DNA sequence for which Tm is calculated.
        :type sequence: Seq
        :return: The calculated melting temperature.
        :rtype: float
        """
        Tm = mt.Tm_NN(sequence, **self.Tm_parameters)
        if self.Tm_salt_correction_parameters is not None:
            Tm += mt.salt_correction(**self.Tm_salt_correction_parameters, seq=sequence)
        if self.Tm_chem_correction_parameters is not None:
            Tm = mt.chem_correction(Tm, **self.Tm_chem_correction_parameters)
        return round(Tm, 4)

    def apply(self, sequence: Seq):
        """Applies the melting temperature filter to a DNA sequence.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and a dictionary containing
                 the actual melting temperature if the condition is met.
        :rtype: (bool, dict)
        """
        Tm = self._get_Tm(sequence)
        if self.Tm_min < Tm < self.Tm_max:
            return True, {"melting_temperature": Tm}
        return False, {}


class SecondaryStructureFilter(PropertyFilterBase):
    """A filter to evaluate the stability of secondary structures formed by a DNA sequence, based on free energy (∆G) at a given temperature.
    Secondary structures can contain for instance hairpins, stacks, bulges or interior loops.
    The minimum free energy of the folded sequence should be weaker (more positive) than the given threshold.

    :param T: The temperature at which the free energy is calculated.
    :type T: float
    :param thr_DG: The threshold free energy value for determining structure stability.
    :type thr_DG: float
    """

    def __init__(self, T: float, thr_DG: float) -> None:
        """Constructor for the SecondaryStructureFilter class."""
        super().__init__()
        self.T = T
        self.thr_DG = thr_DG

    def apply(self, sequence: Seq):
        """Applies the secondary structure stability filter to a DNA sequence.

        :param sequence: The DNA sequence to be checked for secondary structure stability.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and a dictionary containing
                 the actual ∆G value if the condition is met.
        :rtype: (bool, dict)
        """
        DG = dg(sequence, temp=self.T)
        if DG > self.thr_DG:
            return True, {"secondary_structure_DG": DG}
        return False, {}
