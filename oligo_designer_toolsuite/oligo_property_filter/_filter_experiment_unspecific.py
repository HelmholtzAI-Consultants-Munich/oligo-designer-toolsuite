############################################
# imports
############################################
from typing import List, Union

from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import Seq, gc_fraction
from seqfold import dg

from . import PropertyFilterBase


###TODO move this function to utils once database refactor is merged
def check_sequence(seq: str, valid_characters={"A", "C", "T", "G"}) -> bool:
    return all(char.upper() in valid_characters for char in seq)


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
                                The parameter can be either a single sequence or a list of sequences.
    :type prohibited_sequence: str, list[str]
    """

    def __init__(
        self,
        prohibited_sequence: Union[str, List[str]],
    ) -> None:
        """Constructor for the ProhibitedSequenceFilter class."""
        super().__init__()
        if not isinstance(prohibited_sequence, list):
            prohibited_sequence = [prohibited_sequence]
        self.prohibited_sequence = [s.upper() for s in prohibited_sequence]
        # Check that the prohibited sequences are valid DNA sequences.
        for s in self.prohibited_sequence:
            if not check_sequence(s):
                raise ValueError(
                    "Prohibited sequence ({prohibited_sequences}) is not a DNA sequence."
                )

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

    :param base_n: A dictionary where keys are nucleotide bases (A, T, C, G) and values are the minimum
                   number of consecutive repeats that define a homopolymeric run for that base.
    :type base_n: dict
    """

    def __init__(
        self,
        base_n: dict,
    ) -> None:
        """Constructor for the HomopolymericRunsFilter class."""
        super().__init__()
        # check that the nucleotides provided are valid
        for b in base_n.keys():
            if not check_sequence(b):
                raise ValueError("Prohibited sequence ({base}) is not a DNA sequence.")
        # create all homopolymeric runs
        self.homopolymeric_runs = [base.upper() * n for base, n in base_n.items()]

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
    """A filter to check if the last `n_bases` of the 3' terminal end of a DNA sequence contain at least `n_GC` G or C bases.

    :param n_bases: The number of bases from the 3' end of the sequence to check for the presence of G or C.
    :type n_bases: int
    :param n_GC: The minimum number of G or C bases at the 3' end.
    :type n_GC: int
    """

    def __init__(self, n_bases: int, n_GC: int) -> None:
        """Constructor for the GCClampFilter class."""
        super().__init__()
        self.n_bases = n_bases
        self.n_GC = n_GC

    def apply(self, sequence: Seq):
        """Applies the GC clamp filter to the 3' end of a DNA sequence and returns True if there is a GC clamp.
        A GC clamp is the presence of a guanine (G) or cytosine (C) base in the last n bases (the 3' end) of an oligo.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and an empty dictionary.
        :rtype: (bool, dict)
        """
        GC_counetr = 0
        for i in range(1, self.n_bases + 1):
            if sequence.upper()[-i] == "G" or sequence.upper()[-i] == "C":
                GC_counetr += 1
            if GC_counetr >= self.n_GC:
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


class HomodimerFilter(PropertyFilterBase):
    """A filter to evaluate the potential formation of homodimers. A homodimer is formed when two strands of DNA
    bind to each other due to complementary sequences. This filter calculates the longest self-complementary
    sequence within a given DNA sequence and compares it against a maximum allowable length to prevent homodimer
    formation.

    :param max_len_selfcomp: The maximum length of self-complementary sequence allowed to avoid homodimer formation.
    :type max_len_selfcomp: int
    """

    def __init__(self, max_len_selfcomp: int) -> None:
        """Constructor for the HomodimerFilter class."""
        super().__init__()
        self.max_len_selfcomp = max_len_selfcomp

    def _calculate_len_selfcomp(self, sequence: Seq):
        """Calculates the length of the longest self-complementary sequence in the given DNA sequence.

        :param sequence: The DNA sequence to analyze for self-complementary sequences.
        :type sequence: Seq
        :return: The length of the longest self-complementary sequence found in the DNA sequence.
        :rtype: int
        """
        # we want to check if the reverse of our sequence is complementary to itself, e.g.
        # 5' - TAA CAA TAT ATA TTG TTA - 3' and it's reverse
        # 3' - ATT GTT ATA TAT AAC AAT - 5' are complementary to each other
        # but since we are comparing strings, we take the reverse complement,
        # which should be the exact same sequence in this case
        sequence_revcomp = sequence.reverse_complement()

        # initialize counters
        len_selfcomp_sub = 0
        len_selfcomp = 0

        # iterate through sequences
        for i in range(len(sequence)):
            if sequence[i] != sequence_revcomp[i]:
                len_selfcomp_sub = 0
            else:
                len_selfcomp_sub += 1
            len_selfcomp = max(len_selfcomp, len_selfcomp_sub)
        return len_selfcomp

    def apply(self, sequence: Seq):
        """Applies the homodimer filter to a DNA sequence to evaluate its potential for homodimer formation.

        :param sequence: The DNA sequence to be checked for homodimer formation potential.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and a dictionary containing
                 the maximum self-complementary length if the condition is met.
        :rtype: (bool, dict)
        """
        len_selfcomp = self._calculate_len_selfcomp(sequence)
        if len_selfcomp <= self.max_len_selfcomp:
            return True, {"len_selfcomp": len_selfcomp}
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
