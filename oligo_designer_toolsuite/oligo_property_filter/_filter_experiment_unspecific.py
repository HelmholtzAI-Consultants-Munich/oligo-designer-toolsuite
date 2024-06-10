############################################
# imports
############################################

from typing import List, Union

from Bio.SeqUtils import Seq

from oligo_designer_toolsuite.database import OligoAttributes
from oligo_designer_toolsuite.oligo_property_filter import PropertyFilterBase

from ..utils._checkers import check_if_dna_sequence

############################################
# Oligo Property Filter Classes
############################################


class SoftMaskedSequenceFilter(PropertyFilterBase):
    """A filter to check if a DNA sequence is soft-masked (contains lowercase letters)."""

    def __init__(self):
        """Constructor for the SoftMaskedSequenceFilter class."""
        super().__init__()

    def apply(self, sequence: Seq):
        """
        Applies the soft mask filter to a DNA sequence and returns True if the sequence does not contain lower-case letters.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        if any(c.islower() for c in sequence):
            return False
        return True


class HardMaskedSequenceFilter(PropertyFilterBase):
    """A filter to check if a DNA sequence contains a specific mask character (default is "N").

    :param mask: The mask character to check for in the sequence. Default is "N".
    :type mask: str
    """

    def __init__(self, mask: str = "N"):
        """Constructor for the HardMaskedSequenceFilter class."""
        super().__init__()
        self.mask = mask

    def apply(self, sequence: Seq):
        """
        Applies the hard mask filter to a DNA sequence and returns True if the sequence does not contain the mask letter.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        if self.mask in sequence:
            return False
        return True


class ProhibitedSequenceFilter(PropertyFilterBase):
    """A filter that checks for the absence of a specified prohibited sequence in a DNA sequence.
    The filter will reject any sequence that contains the prohibited sequence.

    :param prohibited_sequence: The DNA sequence that is to be prohibited. This sequence should
                                be a valid DNA sequence consisting of characters A, T, C, and G only.
                                The parameter can be either a single sequence or a list of sequences.
    :type prohibited_sequence: str, list[str]
    """

    def __init__(self, prohibited_sequence: Union[str, List[str]]):
        """Constructor for the ProhibitedSequenceFilter class."""
        super().__init__()
        if not isinstance(prohibited_sequence, list):
            prohibited_sequence = [prohibited_sequence]
        self.prohibited_sequence = [s.upper() for s in prohibited_sequence]
        # Check that the prohibited sequences are valid DNA sequences.
        for s in self.prohibited_sequence:
            if not check_if_dna_sequence(s):
                raise ValueError("Prohibited sequence ({prohibited_sequences}) is not a DNA sequence.")

    def apply(self, sequence: Seq):
        """
        Applies the filter to a given DNA sequence to check if it contains the prohibited sequence.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        for s in self.prohibited_sequence:
            if s in sequence.upper():
                return False
        return True


class HomopolymericRunsFilter(PropertyFilterBase):
    """A filter hat checks for the absence of a specified homopolymeric run in a DNA sequence.
    A homopolymeric run is defined as a sequence where the same nucleotide base repeats consecutively.

    :param base_n: A dictionary where keys are nucleotide bases (A, T, C, G) and values are the minimum
                   number of consecutive repeats that define a homopolymeric run for that base.
    :type base_n: dict
    """

    def __init__(self, base_n: dict):
        """Constructor for the HomopolymericRunsFilter class."""
        super().__init__()
        # check that the nucleotides provided are valid
        for b in base_n.keys():
            if not check_if_dna_sequence(b):
                raise ValueError("Prohibited sequence ({base}) is not a DNA sequence.")
        # create all homopolymeric runs
        self.homopolymeric_runs = [base.upper() * n for base, n in base_n.items()]

    def apply(self, sequence: Seq):
        """Applies the filter to a given DNA sequence to check if it contains a homopolymeric run.

        :param sequence: The DNA sequence to be checked for homopolymeric runs.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        for homopolymeric_run in self.homopolymeric_runs:
            if homopolymeric_run in sequence.upper():
                return False
        return True


class FivePrimeSequenceFilter(PropertyFilterBase):
    """A filter to check the presence or absence of a specified sequence at the 5'-end of a DNA sequence.

    :param five_prime_sequence: The sequence to check at the 5'-end of the DNA sequence.
    :type five_prime_sequence: str
    :param remove: If True, sequences starting with the specified sequence are filtered out. If False, only sequences starting with the specified sequence are retained.
    :type remove: bool
    """

    def __init__(self, five_prime_sequence: str, remove: bool = True):
        """Constructor for the FivePrimeSequenceFilter class."""
        super().__init__()
        self.five_prime_sequence = five_prime_sequence.upper()
        self.remove = remove

    def apply(self, sequence: str):
        """Applies the 5'-end sequence filter to a DNA sequence and eitehr keeps or removes the matching sequence, dependend on the parameter "remove".

        :param sequence: The DNA sequence to be checked.
        :type sequence: str
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        if self.remove:
            if sequence.upper().startswith(self.five_prime_sequence):
                return False
            return True
        else:
            if sequence.upper().startswith(self.five_prime_sequence):
                return True
            return False


class ThreePrimeSequenceFilter(PropertyFilterBase):
    """A filter to check the presence or absence of a specified sequence at the 3'-end of a DNA sequence.

    :param three_prime_sequence: The sequence to check at the 3'-end of the DNA sequence.
    :type three_prime_sequence: str
    :param remove: If True, sequences ending with the specified sequence are filtered out. If False, only sequences ending with the specified sequence are retained.
    :type remove: bool
    """

    def __init__(self, three_prime_sequence: str, remove: bool = True):
        """Constructor for the ThreePrimeSequenceFilter class."""
        super().__init__()
        self.three_prime_sequence = three_prime_sequence.upper()
        self.remove = remove

    def apply(self, sequence: str):
        """Applies the 3'-end sequence filter to a DNA sequence and eitehr keeps or removes the matching sequence, dependend on the parameter "remove".

        :param sequence: The DNA sequence to be checked.
        :type sequence: str
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        if self.remove:
            if sequence.upper().endswith(self.three_prime_sequence):
                return False
            return True
        else:
            if sequence.upper().endswith(self.three_prime_sequence):
                return True
            return False


class GCContentFilter(PropertyFilterBase):
    """A filter to check if the GC content of a DNA sequence falls within a specified range [GC_content_min, GC_content_max].

    :param GC_content_min: The minimum acceptable GC content as a percentage.
    :type GC_content_min: float
    :param GC_content_max: The maximum acceptable GC content as a percentage.
    :type GC_content_max: float
    """

    def __init__(self, GC_content_min: float, GC_content_max: float):
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
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        GC_content = OligoAttributes._calc_GC_content(sequence)
        if self.GC_content_min < GC_content < self.GC_content_max:
            return True
        return False


class GCClampFilter(PropertyFilterBase):
    """A filter to check if the last `n_bases` of the 3' terminal end of a DNA sequence contain at least `n_GC` G or C bases.

    :param n_bases: The number of bases from the 3' end of the sequence to check for the presence of G or C.
    :type n_bases: int
    :param n_GC: The minimum number of G or C bases at the 3' end.
    :type n_GC: int
    """

    def __init__(self, n_bases: int, n_GC: int):
        """Constructor for the GCClampFilter class."""
        super().__init__()
        self.n_bases = n_bases
        self.n_GC = n_GC

    def apply(self, sequence: Seq):
        """Applies the GC clamp filter to the 3' end of a DNA sequence and returns True if there is a GC clamp.
        A GC clamp is the presence of a guanine (G) or cytosine (C) base in the last n bases (the 3' end) of an oligo.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        GC_counetr = 0
        for i in range(1, self.n_bases + 1):
            if sequence.upper()[-i] == "G" or sequence.upper()[-i] == "C":
                GC_counetr += 1
            if GC_counetr >= self.n_GC:
                return True
        return False


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
    ):
        """Constructor for the MeltingTemperatureNNFilter class."""
        super().__init__()
        if Tm_max <= Tm_min:
            raise ValueError("Tm_max is lower that Tm_min!")
        self.Tm_min = Tm_min
        self.Tm_max = Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters

    def apply(self, sequence: Seq):
        """Applies the melting temperature filter to a DNA sequence.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        Tm = OligoAttributes._calc_TmNN(
            sequence,
            self.Tm_parameters,
            self.Tm_salt_correction_parameters,
            self.Tm_chem_correction_parameters,
        )
        if self.Tm_min < Tm < self.Tm_max:
            return True
        return False


class SelfComplementFilter(PropertyFilterBase):
    """A filter to evaluate the potential formation of self-complementary sequences.

    This filter calculates the longest self-complementary sequence within a given DNA sequence
    and compares it against a maximum allowable length.

    :param max_len_selfcomp: The maximum length of self-complementary sequence allowed.
    :type max_len_selfcomp: int
    """

    def __init__(self, max_len_selfcomp: int):
        """Constructor for the SelfComplementFilter class."""
        super().__init__()
        self.max_len_selfcomp = max_len_selfcomp

    def apply(self, sequence: Seq):
        """Applies the self complement filter to a DNA sequence.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        len_selfcomp = OligoAttributes._calc_length_selfcomplement(sequence)
        if len_selfcomp <= self.max_len_selfcomp:
            return True
        return False


class ComplementFilter(PropertyFilterBase):
    """A filter to evaluate the potential formation of complementary sequences between two DNA sequences.

    This filter calculates the longest complementary sequence between two DNA sequences and compares it against
    a maximum allowable length.

    :param max_len_complement: The maximum length of complementary sequence allowed.
    :type max_len_complement: int
    """

    def __init__(self, max_len_complement: int):
        """Constructor for the ComplementFilter class."""
        super().__init__()
        self.max_len_complement = max_len_complement

    def apply(self, sequence1: Seq, sequence2: Seq):
        """Applies the complement filter to a pair of DNA sequences.

        :param sequence1: The first DNA sequence to be checked.
        :type sequence1: Seq
        :param sequence2: The second DNA sequence to be checked.
        :type sequence2: Seq
        :return: True if the sequences meet the criteria, False otherwise.
        :rtype: bool
        """
        len_complement = OligoAttributes._calc_length_complement(sequence1, sequence2)
        if len_complement <= self.max_len_complement:
            return True
        return False


class SecondaryStructureFilter(PropertyFilterBase):
    """A filter to evaluate the stability of secondary structures formed by a DNA sequence, based on free energy (∆G) at a given temperature.
    Secondary structures can contain for instance hairpins, stacks, bulges or interior loops.
    The minimum free energy of the folded sequence should be weaker (more positive) than the given threshold.

    :param T: The temperature at which the free energy is calculated.
    :type T: float
    :param thr_DG: The threshold free energy value for determining structure stability.
    :type thr_DG: float
    """

    def __init__(self, T: float, thr_DG: float):
        """Constructor for the SecondaryStructureFilter class."""
        super().__init__()
        self.T = T
        self.thr_DG = thr_DG

    def apply(self, sequence: Seq):
        """Applies the secondary structure stability filter to a DNA sequence.

        :param sequence: The DNA sequence to be checked for secondary structure stability.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
        DG_secondary_structure = OligoAttributes._calc_DG_secondary_structure(sequence, self.T)
        if DG_secondary_structure > self.thr_DG:
            return True
        return False
