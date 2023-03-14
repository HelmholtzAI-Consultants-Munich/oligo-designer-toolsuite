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
    Secondary_struct
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

    filters = [
        masked_sequences_filter,
        GC_content_filter,
        Tm_filter,
        Tm_arms_filter,
    ]
    property_filter = PropertyFilter(filters=filters)

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = property_filter._filter_sequence(seq_keep)
    assert (
            res == True
    ), f"error: A sequence ({seq_keep}) fulfilling all conditions has not been accepted!"


def test_ConsecutiveRepeats_filter():  # add tests for property filters
    repeat_num_max = ConsecutiveRepeats(3)

    seq_remove = Seq("CTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, _ = repeat_num_max.apply(seq_remove)
    assert (
            res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the ConsecutiveRepeat condition has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = repeat_num_max.apply(seq_keep)
    assert (
            res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the ConsecutiveRepeat conditions has not been accepted!"

def test_GC_Clamp_filter():  # add tests for property filters
    GC_Clamp = GCClamp(2)

    seq_remove = Seq("TGGGCCTTTCCAAGCCCCCATTTGAGCTA")
    res, _ = GC_Clamp.apply(seq_remove)
    assert (
            res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the GC_Clamp condition has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGC")
    res, _ = GC_Clamp.apply(seq_keep)
    assert (
            res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the GC_Clamp conditions has not been accepted!"

def test_secondary_structure_filter():  # add tests for property filters
    secondary_structure = Secondary_struct(76.0, 0.0)

    seq_remove = Seq("GCUUUAGAGAUCGUUUCGAAUGAUAAUCGUUCGAAACGUUCUCCGAAGC")
    res, _ = GCClamp.apply(seq_remove)
    assert (
            res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the secondary_structure condition has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = GCClamp.apply(seq_keep)
    assert (
            res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the secondary_structure conditions has not been accepted!"



def test_merfish_property_filter():  # add tests for property filters
    repeat_num_max = ConsecutiveRepeats(3)
    GC_Clamp = GCClamp(2)
    secondary_structure = Secondary_struct(76.0, 0.0)
    GC_content_filter = GCContent(GC_content_min=40, GC_content_max=60)
    Tm_filter = MeltingTemperatureNN(
        Tm_min=52,
        Tm_max=67,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )

    filters = [
        repeat_num_max,
        GC_Clamp,
        secondary_structure,
        Tm_filter,
        GC_content_filter
    ]

    property_filter = PropertyFilter(filters=filters)

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = property_filter._filter_sequence(seq_keep)
    assert (
            res == True
    ), f"error: A sequence ({seq_keep}) fulfilling all conditions has not been accepted!"
