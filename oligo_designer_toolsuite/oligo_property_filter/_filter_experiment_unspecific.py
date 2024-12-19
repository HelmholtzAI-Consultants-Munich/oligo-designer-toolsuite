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
    def __init__(self) -> None:
        """Constructor for the SoftMaskedSequenceFilter class."""
        super().__init__()

    def apply(self, sequence: Seq) -> bool:
        if any(c.islower() for c in sequence):
            return False
        return True


class HardMaskedSequenceFilter(PropertyFilterBase):
    def __init__(self, mask: str = "N") -> None:
        """Constructor for the HardMaskedSequenceFilter class."""
        super().__init__()
        self.mask = mask

    def apply(self, sequence: Seq) -> bool:
        if self.mask in sequence:
            return False
        return True


class ProhibitedSequenceFilter(PropertyFilterBase):
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
        for s in self.prohibited_sequence:
            if s in sequence.upper():
                return False
        return True


class HomopolymericRunsFilter(PropertyFilterBase):
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
        for homopolymeric_run in self.homopolymeric_runs:
            if homopolymeric_run in sequence.upper():
                return False
        return True


class FivePrimeSequenceFilter(PropertyFilterBase):

    def __init__(self, five_prime_sequence: str, remove: bool = True) -> None:
        """Constructor for the FivePrimeSequenceFilter class."""
        super().__init__()
        self.five_prime_sequence = five_prime_sequence.upper()
        self.remove = remove

    def apply(self, sequence: str) -> bool:

        if self.remove:
            if sequence.upper().startswith(self.five_prime_sequence):
                return False
            return True
        else:
            if sequence.upper().startswith(self.five_prime_sequence):
                return True
            return False


class ThreePrimeSequenceFilter(PropertyFilterBase):

    def __init__(self, three_prime_sequence: str, remove: bool = True) -> None:
        """Constructor for the ThreePrimeSequenceFilter class."""
        super().__init__()
        self.three_prime_sequence = three_prime_sequence.upper()
        self.remove = remove

    def apply(self, sequence: str) -> bool:

        if self.remove:
            if sequence.upper().endswith(self.three_prime_sequence):
                return False
            return True
        else:
            if sequence.upper().endswith(self.three_prime_sequence):
                return True
            return False


class GCContentFilter(PropertyFilterBase):

    def __init__(self, GC_content_min: float, GC_content_max: float) -> None:
        """Constructor for the GCContentFilter class."""
        super().__init__()
        if GC_content_max <= GC_content_min:
            raise ValueError("GC_content_max is lower that GC_content_min!")
        self.GC_content_min = GC_content_min
        self.GC_content_max = GC_content_max

    def apply(self, sequence: Seq) -> bool:

        GC_content = OligoAttributes._calc_GC_content(sequence)
        if self.GC_content_min < GC_content < self.GC_content_max:
            return True
        return False


class GCClampFilter(PropertyFilterBase):

    def __init__(self, n_bases: int, n_GC: int) -> None:
        """Constructor for the GCClampFilter class."""
        super().__init__()
        self.n_bases = n_bases
        self.n_GC = n_GC

    def apply(self, sequence: Seq) -> bool:
        GC_counetr = 0
        for i in range(1, self.n_bases + 1):
            if sequence.upper()[-i] == "G" or sequence.upper()[-i] == "C":
                GC_counetr += 1
            if GC_counetr >= self.n_GC:
                return True
        return False


class MeltingTemperatureNNFilter(PropertyFilterBase):

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

    def __init__(self, max_len_selfcomplement: int) -> None:
        """Constructor for the SelfComplementFilter class."""
        super().__init__()
        self.max_len_selfcomplement = max_len_selfcomplement

    def apply(self, sequence: Seq) -> bool:

        len_selfcomp = OligoAttributes._calc_length_complement(sequence, sequence[::-1])
        if len_selfcomp <= self.max_len_selfcomplement:
            return True
        return False


class ComplementFilter(PropertyFilterBase):
    def __init__(self, comparison_sequence: Seq, max_len_complement: int) -> None:
        """Constructor for the ComplementFilter class."""
        super().__init__()
        self.max_len_complement = max_len_complement
        self.comparison_sequence = comparison_sequence

    def apply(self, sequence: Seq) -> bool:
        len_complement = OligoAttributes._calc_length_complement(sequence, self.comparison_sequence)
        if len_complement <= self.max_len_complement:
            return True
        return False


class SecondaryStructureFilter(PropertyFilterBase):
    def __init__(self, T: float, thr_DG: float) -> None:
        """Constructor for the SecondaryStructureFilter class."""
        super().__init__()
        self.T = T
        self.thr_DG = thr_DG

    def apply(self, sequence: Seq) -> bool:
        DG_secondary_structure = OligoAttributes._calc_DG_secondary_structure(sequence, self.T)
        if DG_secondary_structure > self.thr_DG:
            return True
        return False
