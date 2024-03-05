############################################
# imports
############################################


from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.oligo_property_filter import (
    FivePrimeSequenceFilter,
    GCClampFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomodimerFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PadlockArmsFilter,
    ProhibitedSequenceFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SoftMaskedSequenceFilter,
    ThreePrimeSequenceFilter,
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
    "saltcorr": 7,
    "Na": 50,  # [mM]
    "K": 75,  # [mM]
    "Tris": 20,  # [mM]
    "Mg": 10,  # [mM]
    "dNTPs": 0,
}

Tm_chem_correction_parameters = {
    "DMSO": 0,
    "DMSOfactor": 0.75,
    "fmdfactor": 0.65,
    "fmdmethod": 1,
    "GC": None,
    "fmd": 20,
}

Tm_salt_correction_parameters = {
    "method": 7,
    "Na": 50,  # [mM]
    "K": 75,  # [mM]
    "Tris": 20,  # [mM]
    "Mg": 10,  # [mM]
    "dNTPs": 0,
}

############################################
# Tests
############################################


def test_masked_sequence_filters():
    """Test if masked sequence filters work, e.g. sequences containing a mask string or lower case letters should be removed."""
    softmasked_sequence_filter = SoftMaskedSequenceFilter()

    seq_remove = Seq("TGTCGGATCTCcTCAACAAGCTGGTCtTGA")
    res, _ = softmasked_sequence_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [SoftMaskedSequenceFilter]"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, feature = softmasked_sequence_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [SoftMaskedSequenceFilter]"
    print(feature)

    hardmasked_sequence_filter1 = HardMaskedSequenceFilter(mask="N")

    seq_remove = Seq("TGTCGGATCTCNTCAACAAGCTGGTCNTGA")
    res, _ = hardmasked_sequence_filter1.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [HardMaskedSequenceFilter]"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, feature = hardmasked_sequence_filter1.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [HardMaskedSequenceFilter]"
    print(feature)

    hardmasked_sequence_filter2 = HardMaskedSequenceFilter(mask="Q")

    seq_remove = Seq("TGTCGGATCTCQTCAACAAGCTGGTCQTGA")
    res, _ = hardmasked_sequence_filter2.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [HardMaskedSequenceFilter]"

    seq_keep = Seq("TGTCGGATCTCTNNAACAAGCTGGTCATGA")
    res, feature = hardmasked_sequence_filter2.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [HardMaskedSequenceFilter]"
    print(feature)


def test_sequence_content_filters():
    """Test if sequence content filters work, e.g. sequences containing / not containing certain nucleotides should be removed."""
    prohibited_sequence_filter1 = ProhibitedSequenceFilter(prohibited_sequence="ACT")

    seq_remove = Seq("GGGGGGGGGGGGGGACT")
    res, _ = prohibited_sequence_filter1.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [ProhibitedSequenceFilter]"

    seq_keep = Seq("GGGGGGGGGGGGGGATC")
    res, feature = prohibited_sequence_filter1.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [ProhibitedSequenceFilter]"
    print(feature)

    prohibited_sequence_filter2 = ProhibitedSequenceFilter(
        prohibited_sequence=["ACT", "CCGC"]
    )

    seq_remove = Seq("GGGGGGGGGGGGGGACT")
    res, _ = prohibited_sequence_filter2.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [ProhibitedSequenceFilter]"

    seq_remove = Seq("GGGGGGGGGGGGGGCCGC")
    res, _ = prohibited_sequence_filter2.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [ProhibitedSequenceFilter]"

    seq_keep = Seq("GGGGGGGGGGGGGGATC")
    res, feature = prohibited_sequence_filter2.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [ProhibitedSequenceFilter]"
    print(feature)

    homopolymeric_run_filter = HomopolymericRunsFilter(base_n={"A": 4, "C": 5})

    seq_remove = Seq("GGGGGGGGGGGGGGAAAAA")
    res, _ = homopolymeric_run_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted!"

    seq_keep = Seq("GGGGGGGGGGGGGGAAA")
    res, feature = homopolymeric_run_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted!"
    print(feature)

    three_prime_filter = ThreePrimeSequenceFilter(
        three_prime_sequence="TT", remove=False
    )

    seq_remove = Seq("GGGGGGGGGGGGGGAAAAA")
    res, _ = three_prime_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [ThreePrimeSequenceFilter]"

    seq_keep = Seq("GGGGGGGGGGGGGGAAATT")
    res, feature = three_prime_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [ThreePrimeSequenceFilter]"
    print(feature)

    five_prime_filter = FivePrimeSequenceFilter(five_prime_sequence="TT", remove=True)

    seq_remove = Seq("TTGGGGGGGGGGGGGGAAAAA")
    res, _ = five_prime_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [FivePrimeSequenceFilter]"

    seq_keep = Seq("GGGGGGGGGGGGGGAAATT")
    res, feature = five_prime_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [FivePrimeSequenceFilter]"
    print(feature)


def test_GC_content_filters():
    """Test if GC content filters work, e.g. sequences having certain GC content or no GC clamp should be removed."""
    GC_content_filter = GCContentFilter(GC_content_min=40, GC_content_max=60)

    seq_remove = Seq("TCGGGCGGGAGATCCAGGTGGCGCGCAAAG")
    res, _ = GC_content_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [GCContentFilter]"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, feature = GC_content_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [GCContentFilter]"
    print(feature)

    GC_clamp_filter = GCClampFilter(n_bases=3, n_GC=1)

    seq_remove = Seq("TCGGGCGGGAGATCCAGGTGGCGCGCAAAAA")
    res, _ = GC_clamp_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [GCClampFilter]"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGAA")
    res, feature = GC_clamp_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [GCClampFilter]"
    print(feature)


def test_sequence_structure_filters():
    """Test if melting temperature filters work, e.g. sequences having certain Tm or secondary structure probability should be removed."""
    # Test if Tm filter works with default parameters
    Tm_filter1 = MeltingTemperatureNNFilter(Tm_min=52, Tm_max=67, Tm_parameters={})

    seq_remove = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, _ = Tm_filter1.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition with has been accepted! [MeltingTemperatureNNFilter]"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, feature = Tm_filter1.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [MeltingTemperatureNNFilter]"
    print(feature)

    # Test if Tm filter works with user-defined Tm parameters
    Tm_filter2 = MeltingTemperatureNNFilter(
        Tm_min=52,
        Tm_max=67,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
        Tm_salt_correction_parameters=Tm_salt_correction_parameters,
    )

    seq_remove = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, _ = Tm_filter2.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition with user-defined parameters has been accepted! [MeltingTemperatureNNFilter]"

    seq_keep = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTAAAAA")
    res, feature = Tm_filter2.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions with user-defined parameters has not been accepted! [MeltingTemperatureNNFilter]"
    print(feature)

    secondary_structure_filter = SecondaryStructureFilter(T=37, thr_DG=0)

    seq_remove = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, _ = secondary_structure_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition with has been accepted! [SecondaryStructureFilter]"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, feature = secondary_structure_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [SecondaryStructureFilter]"
    print(feature)

    homodimer_filer = HomodimerFilter(max_len_selfcomp=6)

    seq_remove = Seq("TAACAATATATATTGTTA")
    res, _ = homodimer_filer.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition with has been accepted!"

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, feature = homodimer_filer.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted!"
    print(feature)


def test_experiment_specific_filters():
    """Test if experiment specific filters work, e.g. sequences not forming proper padlock arms should be removed."""
    padlock_arms_filter = PadlockArmsFilter(
        arm_length_min=5,
        arm_Tm_dif_max=5,
        arm_Tm_min=40,
        arm_Tm_max=60,
        Tm_parameters=Tm_parameters,
        Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )

    seq_remove = Seq("TGTCGGATCTCTTCAACAAGCTGGTCAT")
    res, _ = padlock_arms_filter.apply(seq_remove)
    assert (
        res == False
    ), f"error: A sequence ({seq_remove}) not fulfilling the condition with has been accepted! [PadlockArmsFilter]"

    seq_keep = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
    res, feature = padlock_arms_filter.apply(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [PadlockArmsFilter]"
    print(feature)


def test_property_filter():
    """Test if property filter correctly applies all given filters."""
    masked_sequences_filter = HardMaskedSequenceFilter(mask="N")
    GC_content_filter = GCContentFilter(GC_content_min=40, GC_content_max=60)
    GC_clamp_filter = GCClampFilter(n_bases=2, n_GC=1)

    Tm_filter = MeltingTemperatureNNFilter(
        Tm_min=52,
        Tm_max=67,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )
    secondary_structure_filter = SecondaryStructureFilter(37, -5)

    Tm_arms_filter = PadlockArmsFilter(
        arm_length_min=5,
        arm_Tm_dif_max=5,
        arm_Tm_min=40,
        arm_Tm_max=60,
        Tm_parameters=Tm_parameters,
        Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )

    filters = [
        masked_sequences_filter,
        GC_content_filter,
        GC_clamp_filter,
        Tm_filter,
        secondary_structure_filter,
        Tm_arms_filter,
    ]
    property_filter = PropertyFilter(filters=filters)

    seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
    res, _ = property_filter._filter_sequence(seq_keep)
    assert (
        res == True
    ), f"error: A sequence ({seq_keep}) fulfilling all conditions has not been accepted!"
