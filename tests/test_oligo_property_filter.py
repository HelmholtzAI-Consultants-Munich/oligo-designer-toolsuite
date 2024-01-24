############################################
# imports
############################################


from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.oligo_property_filter import (
    GCContent,
    MaskedSequences,
    MeltingTemperatureNN,
    GCClamp,
    ConsecutiveRepeats,
    SecondaryStructure,
    ThreePrimeSequence,
    FivePrimeSequence,
    RepeatMaskingFilter,
)

from oligo_designer_toolsuite.oligo_property_filter import (
    PadlockArms,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    PropertyFilter,
)

############################################
# Global Parameters
############################################

Tm_parameters = {
    "check": True,
    "strict": True,
    "c_seq": None,
    "shift": 0,
    "nn_table": getattr(mt, "DNA_NN3"),
    "tmm_table": getattr(mt, "DNA_TMM1"),
    "imm_table": getattr(mt, "DNA_IMM1"),
    "de_table": getattr(mt, "DNA_DE1"),
    "dnac1": 50,  # [nM]
    "dnac2": 0,
    "selfcomp": False,
    "dNTPs": 0,
    "saltcorr": 7,
    "Na": 1.25,  # [mM]
    "K": 75,  # [mM]
    "Tris": 20,  # [mM]
    "Mg": 10,  # [mM]
}

Tm_chem_correction_parameters = {
    "DMSO": 0,
    "DMSOfactor": 0.75,
    "fmdfactor": 0.65,
    "fmdmethod": 1,
    "GC": None,
    "fmd": 20,
}

############################################
# Tests
############################################


def test_masked_sequences_filter():
    """Test if masked sequences filter works, i.e. sequences containing a mask string should be removed."""
    masked_sequences_filter = MaskedSequences(mask="N")

    seq_remove = Seq("TGTCGGATCTCNTCAACAAGCTGGTCNTGA")
    res, _ = masked_sequences_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the mask condition has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = masked_sequences_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the mask conditions has not been accepted!"


def test_GC_content_filter():
    """Tets if GC content filter removes sequences that are not within min and max bounds."""
    GC_content_filter = GCContent(GC_content_min=40, GC_content_max=60)

    seq_remove = Seq("TCGGGCGGGAGATCCAGGTGGCGCGCAAAG")
    res, _ = GC_content_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the GC content condition has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = GC_content_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the GC content conditions has not been accepted!"


def test_Tm_filter():
    """Test if melting temperature filter removes sequences that are not within min and max bounds."""
    # Test if Tm filter works with default parameters
    Tm_filter = MeltingTemperatureNN(Tm_min=52, Tm_max=67, Tm_parameters={})

    seq_remove = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, _ = Tm_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the Tm condition with has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = Tm_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the Tm conditions has not been accepted!"

    # Test if Tm filter works with user-defined Tm parameters
    Tm_filter = MeltingTemperatureNN(
        Tm_min=52,
        Tm_max=67,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )

    seq_remove = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, _ = Tm_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the Tm condition with user-defined parameters has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = Tm_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the Tm conditions with user-defined parameters has not been accepted!"


def test_Tm_arms_filter():
    """Test if filter for padlock arm melting temperature removes sequences that don't fullfill
    user-defined requirements concerning melting temperature and arm length.
    """
    Tm_arms_filter = PadlockArms(
        min_arm_length=10,
        max_arm_Tm_dif=2,
        arm_Tm_min=40,
        arm_Tm_max=43,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )

    seq_remove = Seq("CTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, _ = Tm_arms_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the Tm condition has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = Tm_arms_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the Tm conditions has not been accepted!"


def test_consecutive_repeats_filter():  # add tests for property filters
    consecutive_repeats_filter = ConsecutiveRepeats(3)

    seq_remove = Seq("CTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, _ = consecutive_repeats_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the ConsecutiveRepeat condition has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = consecutive_repeats_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the ConsecutiveRepeat conditions has not been accepted!"


def test_GC_clamp_filter():  # add tests for property filters
    GC_clamp_filter = GCClamp(2)

    seq_remove = Seq("TGGGCCTTTCCAAGCCCCCATTTGAGCTA")
    res, _ = GC_clamp_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the GC_Clamp condition has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGC")
    res, _ = GC_clamp_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the GC_Clamp conditions has not been accepted!"


def test_secondary_structure_filter():  # add tests for property filters
    secondary_structure_filter = SecondaryStructure(37, -5)

    seq_remove = Seq("GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC")
    res, _ = secondary_structure_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the secondary_structure condition has been accepted!"

    seq_keep = Seq("ATGACAGTATGACGATGACGATGACGATGAC")
    res, delta_g = secondary_structure_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep} with {delta_g}) fulfilling the secondary_structure conditions has not been accepted!"


def test_three_prime_sequence_filter():
    """Test if the 3' sequence filter works, i.e., sequences with a certain 3' end should be removed."""
    three_prime_sequence_filter = ThreePrimeSequence(three_prime_sequence="GA")

    seq_remove = "TGTCGGATCTCTTCAACAAGCTGGTCATGA"
    res, _ = three_prime_sequence_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) with a matching 3' end has been accepted!"

    seq_keep = "TGTCGGATCTCNTCAACAAGCTGGTCNTGG"
    res, _ = three_prime_sequence_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) with a non-matching 3' end has not been accepted!"


def test_five_prime_sequence_filter():
    """Test if the 5' sequence filter works, i.e., sequences with a certain 5' end should be removed."""
    five_prime_sequence_filter = FivePrimeSequence(five_prime_sequence="TT")

    seq_remove = "TTTCGGATCCGAATNCAAGCTGGTCATGA"
    res, _ = five_prime_sequence_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) with a matching 5' end has been accepted!"

    seq_keep = "TCGGATCCGAATNCAAGCTGGTCATGA"
    res, _ = five_prime_sequence_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) with a non-matching 5' end has not been accepted!"


def test_repeat_masking_filter():
    """Test if the RepeatMaskingFilter correctly identifies and filters sequences with soft-masked regions."""
    repeat_masking_filter = RepeatMaskingFilter()

    # Sequence with lower-case letters indicating soft-masked regions
    seq_with_mask = "TGTCGGatctntCAACaagctggtcATGA"
    res, _ = repeat_masking_filter.apply(seq_with_mask)
    assert (
        not res
    ), f"Error: A sequence with soft-masked regions ({seq_with_mask}) has been accepted!"

    # Sequence without lower-case letters, should be kept
    seq_without_mask = "TGTCGGATCTNTCAACAAGCTGGTCATGA"
    res, _ = repeat_masking_filter.apply(seq_without_mask)
    assert (
        res
    ), f"Error: A sequence without soft-masked regions ({seq_without_mask}) has not been accepted!"


def test_property_filter():
    """Test if property filter correctly applies all given filters."""
    masked_sequences_filter = MaskedSequences(mask="N")
    GC_content_filter = GCContent(GC_content_min=40, GC_content_max=60)
    Tm_filter = MeltingTemperatureNN(
        Tm_min=52,
        Tm_max=67,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )
    Tm_arms_filter = PadlockArms(
        min_arm_length=10,
        max_arm_Tm_dif=2,
        arm_Tm_min=40,
        arm_Tm_max=43,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )
    consecutive_repeats_filter = ConsecutiveRepeats(3)
    GC_clamp_filter = GCClamp(2)
    secondary_structure_filter = SecondaryStructure(37, -5)

    filters = [
        masked_sequences_filter,
        GC_content_filter,
        Tm_filter,
        Tm_arms_filter,
        consecutive_repeats_filter,
        GC_clamp_filter,
        secondary_structure_filter,
    ]
    property_filter = PropertyFilter(filters=filters)

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = property_filter._filter_sequence(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling all conditions has not been accepted!"
