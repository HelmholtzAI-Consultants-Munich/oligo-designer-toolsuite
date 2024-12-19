############################################
# imports
############################################

from typing import List, Union

from Bio.SeqUtils import Seq

from oligo_designer_toolsuite.database import OligoAttributes
from oligo_designer_toolsuite.oligo_property_filter import PropertyFilterBase

from ..utils._checkers_and_helpers import check_if_dna_sequence

############################################
# Oligo Property Filter Classes
############################################


class SoftMaskedSequenceFilter(PropertyFilterBase):
    """
    A filter class to check for the presence of soft-masked sequences.

    The `SoftMaskedSequenceFilter` class is designed to identify sequences that contain soft-masked regions,
    typically represented by lowercase letters. This class can be used to filter out sequences that have
    been soft-masked, ensuring that only sequences with uppercase characters are passed.
    """

    def __init__(self) -> None:
        """Constructor for the SoftMaskedSequenceFilter class."""
        super().__init__()

    def apply(self, sequence: Seq) -> bool:
        """
        Evaluate the sequence to determine if it includes any lowercase letters.
        If any lowercase characters are present, the method returns `False`, indicating the sequence is soft-masked.
        Otherwise, it returns `True`.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence has no lowercase characters, `False` otherwise.
        :rtype: bool
        """
        if any(c.islower() for c in sequence):
            return False
        return True


class HardMaskedSequenceFilter(PropertyFilterBase):
    """
    A filter class to check for the presence of hard-masked sequences.

    The `HardMaskedSequenceFilter` class is designed to identify sequences that contain specific hard-masked regions,
    represented by a designated character (defaults to "N"). This class can be used to filter out sequences that contain
    such masked regions, ensuring only sequences without the specified mask are passed.

    :param mask: The character used to represent hard-masked regions in sequences, defaults to "N".
    :type mask: str
    """

    def __init__(self, mask: str = "N") -> None:
        """Constructor for the HardMaskedSequenceFilter class."""
        super().__init__()
        self.mask = mask

    def apply(self, sequence: Seq) -> bool:
        """
        Evaluate the sequence to determine if it includes the hard-mask character defined during initialization.
        If the character is present, the method returns `False`, indicating the sequence is hard-masked.
        Otherwise, it returns `True`.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence does not contain the mask character, `False` otherwise.
        :rtype: bool
        """
        if self.mask in sequence:
            return False
        return True


class ProhibitedSequenceFilter(PropertyFilterBase):
    """
    A filter class for identifying and excluding sequences that contain specific prohibited subsequences.

    The `ProhibitedSequenceFilter` class is designed to filter out sequences that contain any of the specified prohibited subsequences.
    This is particularly useful in contexts where certain sequences are known to cause issues or are otherwise undesirable.
    The prohibited sequences are case-insensitive and validated to ensure they are valid DNA sequences.

    :param prohibited_sequence: A single prohibited sequence or a list of prohibited sequences to filter out.
    :type prohibited_sequence: Union[str, List[str]]
    """

    def __init__(self, prohibited_sequence: Union[str, List[str]]) -> None:
        """Constructor for the ProhibitedSequenceFilter class."""
        super().__init__()
        if not isinstance(prohibited_sequence, list):
            prohibited_sequence = [prohibited_sequence]
        self.prohibited_sequence = [s.upper() for s in prohibited_sequence]
        # Check that the prohibited sequences are valid DNA sequences.
        for s in self.prohibited_sequence:
            if not check_if_dna_sequence(s):
                raise ValueError("Prohibited sequence ({prohibited_sequences}) is not a DNA sequence.")

    def apply(self, sequence: Seq) -> bool:
        """
        Evaluate the sequence to determine if it contains any of the prohibited subsequences specified during initialization.
        If a match is found, the method returns `False`, indicating the sequence contains a prohibited subsequence.
        Otherwise, it returns `True`.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence does not contain any prohibited subsequences, `False` otherwise.
        :rtype: bool
        """
        for s in self.prohibited_sequence:
            if s in sequence.upper():
                return False
        return True


class HomopolymericRunsFilter(PropertyFilterBase):
    """
    A filter class for excluding sequences containing homopolymeric runs of specified nucleotides.

    The `HomopolymericRunsFilter` class is used to filter out sequences that contain homopolymeric runs—stretches
    of the same nucleotide repeated multiple times—based on criteria defined by the user.
    The filter is useful in various genomic applications where homopolymeric runs may cause issues.

    :param base_n: A dictionary where keys are nucleotides (e.g., 'A', 'T') and values are the minimum number of repeats to define a homopolymeric run.
    :type base_n: dict
    """

    def __init__(self, base_n: dict) -> None:
        """Constructor for the HomopolymericRunsFilter class."""
        super().__init__()
        # check that the nucleotides provided are valid
        for b in base_n.keys():
            if not check_if_dna_sequence(b):
                raise ValueError("Prohibited sequence ({base}) is not a DNA sequence.")
        # create all homopolymeric runs
        self.homopolymeric_runs = [base.upper() * n for base, n in base_n.items()]

    def apply(self, sequence: Seq) -> bool:
        """
        Evaluate a DNA sequence to determine if it contains any homopolymeric runs,
        which are stretches of the same nucleotide repeated a specified number of times.
        If such a run is found, the method returns `False`, indicating that the sequence should be excluded.
        Otherwise, it returns `True`.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence does not contain any prohibited homopolymeric runs, `False` otherwise.
        :rtype: bool
        """
        for homopolymeric_run in self.homopolymeric_runs:
            if homopolymeric_run in sequence.upper():
                return False
        return True


class FivePrimeSequenceFilter(PropertyFilterBase):
    """
    A filter class for identifying or removing sequences based on a specified 5' (five prime) sequence.

    The `FivePrimeSequenceFilter` class is designed to filter sequences according to
    whether they start with a specific nucleotide sequence at the 5' end.
    The user can choose to remove sequences that match this criterion or to retain only those that do.

    :param five_prime_sequence: The nucleotide sequence to check at the 5' end of each sequence.
    :type five_prime_sequence: str
    :param remove: Whether to remove sequences that start with the specified 5' sequence (`True` to remove, `False` to retain).
    :type remove: bool
    """

    def __init__(self, five_prime_sequence: str, remove: bool = True) -> None:
        """Constructor for the FivePrimeSequenceFilter class."""
        super().__init__()
        self.five_prime_sequence = five_prime_sequence.upper()
        self.remove = remove

    def apply(self, sequence: str) -> bool:
        """
        Evaluate whether the input sequence begins with the specified 5' sequence.
        Depending on the `remove` parameter set during initialization, it will return `True` or `False`
        to indicate whether the sequence passes the filter.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :return: `True` if the sequence should be kept according to the filter criteria, `False` otherwise.
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
    """
    A filter class for identifying or removing sequences based on a specified 3' (three prime) sequence.

    The `ThreePrimeSequenceFilter` class is designed to filter sequences according to
    whether they end with a specific nucleotide sequence at the 3' end.
    The user can choose to remove sequences that match this criterion or to retain only those that do.

    :param three_prime_sequence: The nucleotide sequence to check at the 3' end of each sequence.
    :type three_prime_sequence: str
    :param remove: Whether to remove sequences that end with the specified 3' sequence (`True` to remove, `False` to retain).
    :type remove: bool
    """

    def __init__(self, three_prime_sequence: str, remove: bool = True) -> None:
        """Constructor for the ThreePrimeSequenceFilter class."""
        super().__init__()
        self.three_prime_sequence = three_prime_sequence.upper()
        self.remove = remove

    def apply(self, sequence: str) -> bool:
        """
        Evaluate whether the input sequence ends with the specified 3' sequence.
        Depending on the `remove` parameter set during initialization, it will return `True` or `False`
        to indicate whether the sequence passes the filter.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :return: `True` if the sequence should be kept according to the filter criteria, `False` otherwise.
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
    """
    A filter class for evaluating sequences based on their GC content.

    The `GCContentFilter` class is designed to filter sequences by ensuring that their GC content falls within a specified range.
    This is useful in various applications where the stability and binding properties of sequences are important,
    as GC content can influence these factors.

    :param GC_content_min: The minimum acceptable GC content percentage for a sequence.
    :type GC_content_min: float
    :param GC_content_max: The maximum acceptable GC content percentage for a sequence.
    :type GC_content_max: float
    """

    def __init__(self, GC_content_min: float, GC_content_max: float) -> None:
        """Constructor for the GCContentFilter class."""
        super().__init__()
        if GC_content_max <= GC_content_min:
            raise ValueError("GC_content_max is lower that GC_content_min!")
        self.GC_content_min = GC_content_min
        self.GC_content_max = GC_content_max

    def apply(self, sequence: Seq) -> bool:
        """
        Ccalculate the GC content of the provided sequence and check if it falls between
        the minimum and maximum GC content thresholds set during initialization.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the GC content of the sequence is within the specified range, `False` otherwise.
        :rtype: bool
        """
        GC_content = OligoAttributes._calc_GC_content(sequence)
        if self.GC_content_min < GC_content < self.GC_content_max:
            return True
        return False


class GCClampFilter(PropertyFilterBase):
    """
    A filter class that checks for the presence of a GC clamp at the 3' end of a sequence.

    The `GCClampFilter` class is designed to ensure that a specified number of bases at the 3' end of a sequence
    contain a minimum number of G or C nucleotides, forming a GC clamp.
    GC clamps are known to enhance the stability of DNA binding and are often desired in oligonucleotide design.

    :param n_bases: The number of bases to consider from the 3' end of the sequence.
    :type n_bases: int
    :param n_GC: The minimum number of G or C nucleotides required within the specified number of bases.
    :type n_GC: int
    """

    def __init__(self, n_bases: int, n_GC: int) -> None:
        """Constructor for the GCClampFilter class."""
        super().__init__()
        self.n_bases = n_bases
        self.n_GC = n_GC

    def apply(self, sequence: Seq) -> bool:
        """
        Eevaluate the last `n_bases` of the provided sequence to determine if it
        contains at least `n_GC` G or C nucleotides, indicating the presence of a GC clamp.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence has the required GC clamp, `False` otherwise.
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
    """
    A filter class that evaluates sequences based on their melting temperature (Tm) using the nearest-neighbor thermodynamic model.

    The `MeltingTemperatureNNFilter` class is used to determine if a DNA sequence's melting temperature falls within a specified range.
    This filter is particularly useful in oligonucleotide design where precise control over Tm is critical for effective hybridization.

    :param Tm_min: The minimum acceptable melting temperature.
    :type Tm_min: float
    :param Tm_max: The maximum acceptable melting temperature.
    :type Tm_max: float
    :param Tm_parameters: Parameters for calculating the melting temperature.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_salt_correction_parameters: Parameters for salt correction in Tm calculation (optional).
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Parameters for chemical correction in Tm calculation (optional).
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

    def apply(self, sequence: Seq) -> bool:
        """
        Calculate the melting temperature of the given sequence using the nearest-neighbor method and
        check if it lies between the specified minimum and maximum Tm values.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence's Tm is within the specified range, `False` otherwise.
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
    """
    A filter class for detecting and excluding self-complementary DNA sequences.

    The `SelfComplementFilter` class is designed to identify sequences that are self-complementary,
    which means they can form secondary structures such as hairpins.
    This filter ensures that only sequences with a self-complementarity length below a specified threshold are considered valid.

    :param max_len_selfcomplement: The maximum allowed length for self-complementary sequences.
    :type max_len_selfcomplement: int
    """

    def __init__(self, max_len_selfcomplement: int) -> None:
        """Constructor for the SelfComplementFilter class."""
        super().__init__()
        self.max_len_selfcomplement = max_len_selfcomplement

    def apply(self, sequence: Seq) -> bool:
        """
        Calculate the length of the self-complementary region of a given sequence and
        check if it is within the allowed maximum length.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence's self-complementary length is within the specified limit, `False` otherwise.
        :rtype: bool
        """
        len_selfcomp = OligoAttributes._calc_length_complement(sequence, sequence[::-1])
        if len_selfcomp <= self.max_len_selfcomplement:
            return True
        return False


class ComplementFilter(PropertyFilterBase):
    """
    A filter class for detecting and excluding sequences based on their complementarity to a comparison sequence.

    The `ComplementFilter` class is designed to identify if the length of the complementary overlap between a
    given sequence and a comparison sequence is within a specified maximum limit. It is useful for avoiding
    sequences that might interact with other sequences.

    :param comparison_sequence: The sequence against which complementarity is assessed.
    :type comparison_sequence: Seq
    :param max_len_complement: The maximum allowed length of complementarity.
    :type max_len_complement: int
    """

    def __init__(self, comparison_sequence: Seq, max_len_complement: int) -> None:
        """Constructor for the ComplementFilter class."""
        super().__init__()
        self.max_len_complement = max_len_complement
        self.comparison_sequence = comparison_sequence

    def apply(self, sequence: Seq) -> bool:
        """
        Calculate the length of the complementary overlap between the given sequence
        and the comparison sequence and check if it is within the allowed maximum length.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: True if the overlap length is less than or equal to the maximum allowed length, False otherwise.
        :rtype: bool
        """
        len_complement = OligoAttributes._calc_length_complement(sequence, self.comparison_sequence)
        if len_complement <= self.max_len_complement:
            return True
        return False


class SecondaryStructureFilter(PropertyFilterBase):
    """
    A filter class for excluding sequences based on their secondary structure free energy (ΔG).

    The `SecondaryStructureFilter` class is designed to assess the stability of a sequence's secondary structure
    by evaluating its free energy at a given temperature. Sequences with a ΔG above a specified threshold are
    considered stable and are accepted by this filter.

    :param T: The temperature at which the secondary structure is evaluated, in degrees Celsius.
    :type T: float
    :param thr_DG: The threshold for the free energy (ΔG) of the secondary structure. Sequences with a ΔG value above this threshold are accepted.
    :type thr_DG: float
    """

    def __init__(self, T: float, thr_DG: float) -> None:
        """Constructor for the SecondaryStructureFilter class."""
        super().__init__()
        self.T = T
        self.thr_DG = thr_DG

    def apply(self, sequence: Seq) -> bool:
        """
        Calculate the free energy ΔG of the secondary structure of a given sequence and
        check if it exceeds the defined threshold.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: True if the sequence's secondary structure ΔG is greater than the threshold, indicating that the sequence is acceptable. False otherwise.
        :rtype: bool
        """
        DG_secondary_structure = OligoAttributes._calc_DG_secondary_structure(sequence, self.T)
        if DG_secondary_structure > self.thr_DG:
            return True
        return False
