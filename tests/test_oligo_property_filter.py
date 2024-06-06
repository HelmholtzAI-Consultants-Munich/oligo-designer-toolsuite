############################################
# imports
############################################

import os
import shutil
import unittest

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import OligoDatabase
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
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# setup
############################################

# Global Parameters
TM_PARAMETERS = {
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

TM_PARAMETERS_CHEM_CORRECTION = {
    "DMSO": 0,
    "DMSOfactor": 0.75,
    "fmdfactor": 0.65,
    "fmdmethod": 1,
    "GC": None,
    "fmd": 20,
}

TM_PARAMETERS_SALT_CORRECTION = {
    "method": 7,
    "Na": 50,  # [mM]
    "K": 75,  # [mM]
    "Tris": 20,  # [mM]
    "Mg": 10,  # [mM]
    "dNTPs": 0,
}

############################################
# tests
############################################


class TestMaskedSequenceFilters(unittest.TestCase):
    """Test if masked sequence filters work, e.g. sequences containing a mask string or lower case letters should be removed."""

    def setUp(self):
        self.softmasked_sequence_filter = SoftMaskedSequenceFilter()
        self.hardmasked_sequence_filter_N = HardMaskedSequenceFilter(mask="N")
        self.hardmasked_sequence_filter_Q = HardMaskedSequenceFilter(mask="Q")

    def test_softmasked_filter(self):
        seq_remove = Seq("TGTCGGATCTCcTCAACAAGCTGGTCtTGA")
        res = self.softmasked_sequence_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [SoftMaskedSequenceFilter]"

        seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
        res = self.softmasked_sequence_filter.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [SoftMaskedSequenceFilter]"

    def test_hardmasked_filter_N(self):
        seq_remove = Seq("TGTCGGATCTCNTCAACAAGCTGGTCNTGA")
        res = self.hardmasked_sequence_filter_N.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [HardMaskedSequenceFilter]"

        seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
        res = self.hardmasked_sequence_filter_N.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [HardMaskedSequenceFilter]"

    def test_hardmasked_filter_Q(self):
        seq_remove = Seq("TGTCGGATCTCQTCAACAAGCTGGTCQTGA")
        res = self.hardmasked_sequence_filter_Q.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [HardMaskedSequenceFilter]"

        seq_keep = Seq("TGTCGGATCTCTNNAACAAGCTGGTCATGA")
        res = self.hardmasked_sequence_filter_Q.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [HardMaskedSequenceFilter]"


class TestSequenceContentFilters(unittest.TestCase):
    """Test if sequence content filters work, e.g. sequences containing / not containing certain nucleotides should be removed."""

    def setUp(self):
        self.prohibited_sequence_filter_str = ProhibitedSequenceFilter(prohibited_sequence="ACT")
        self.prohibited_sequence_filter_list = ProhibitedSequenceFilter(prohibited_sequence=["ACT", "CCGC"])
        self.homopolymeric_run_filter = HomopolymericRunsFilter(base_n={"A": 4, "C": 5})
        self.three_prime_filter = ThreePrimeSequenceFilter(three_prime_sequence="TT", remove=False)
        self.five_prime_filter = FivePrimeSequenceFilter(five_prime_sequence="TT", remove=True)

    def test_prohibites_sequence_filter_str(self):
        seq_remove = Seq("GGGGGGGGGGGGGGACT")
        res = self.prohibited_sequence_filter_str.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [ProhibitedSequenceFilter]"

        seq_keep = Seq("GGGGGGGGGGGGGGATC")
        res = self.prohibited_sequence_filter_str.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [ProhibitedSequenceFilter]"

    def test_prohibites_sequence_filter_list(self):
        seq_remove = Seq("GGGGGGGGGGGGGGACT")
        res = self.prohibited_sequence_filter_list.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [ProhibitedSequenceFilter]"

        seq_remove = Seq("GGGGGGGGGGGGGGCCGC")
        res = self.prohibited_sequence_filter_list.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [ProhibitedSequenceFilter]"

        seq_keep = Seq("GGGGGGGGGGGGGGATC")
        res = self.prohibited_sequence_filter_list.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [ProhibitedSequenceFilter]"

    def test_homopolymeric_run_filter(self):
        seq_remove = Seq("GGGGGGGGGGGGGGAAAAA")
        res = self.homopolymeric_run_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted!"

        seq_keep = Seq("GGGGGGGGGGGGGGAAA")
        res = self.homopolymeric_run_filter.apply(seq_keep)
        assert res == True, f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted!"

    def test_three_prime_filter(self):
        seq_remove = Seq("GGGGGGGGGGGGGGAAAAA")
        res = self.three_prime_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [ThreePrimeSequenceFilter]"

        seq_keep = Seq("GGGGGGGGGGGGGGAAATT")
        res = self.three_prime_filter.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [ThreePrimeSequenceFilter]"

    def test_five_prime_filter(self):
        seq_remove = Seq("TTGGGGGGGGGGGGGGAAAAA")
        res = self.five_prime_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [FivePrimeSequenceFilter]"

        seq_keep = Seq("GGGGGGGGGGGGGGAAATT")
        res = self.five_prime_filter.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [FivePrimeSequenceFilter]"


class TestGCContentFilters(unittest.TestCase):
    """Test if GC content filters work, e.g. sequences having certain GC content or no GC clamp should be removed."""

    def setUp(self):
        self.GC_content_filter = GCContentFilter(GC_content_min=40, GC_content_max=60)
        self.GC_clamp_filter = GCClampFilter(n_bases=3, n_GC=1)

    def test_GC_content_filter(self):
        seq_remove = Seq("TCGGGCGGGAGATCCAGGTGGCGCGCAAAG")
        res = self.GC_content_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [GCContentFilter]"

        seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
        res = self.GC_content_filter.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [GCContentFilter]"

    def test_GC_clamp_filter(self):
        seq_remove = Seq("TCGGGCGGGAGATCCAGGTGGCGCGCAAAAA")
        res = self.GC_clamp_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition has been accepted! [GCClampFilter]"

        seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGAA")
        res = self.GC_clamp_filter.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [GCClampFilter]"


class TestSequenceStructureFilters(unittest.TestCase):
    """Test if melting temperature filters work, e.g. sequences having certain Tm or secondary structure probability should be removed."""

    def setUp(self):
        self.Tm_filter_default = MeltingTemperatureNNFilter(Tm_min=52, Tm_max=67, Tm_parameters={})
        self.Tm_filter_user = MeltingTemperatureNNFilter(
            Tm_min=52,
            Tm_max=67,
            Tm_parameters=TM_PARAMETERS,
            Tm_chem_correction_parameters=TM_PARAMETERS_CHEM_CORRECTION,
            Tm_salt_correction_parameters=TM_PARAMETERS_SALT_CORRECTION,
        )
        self.secondary_structure_filter = SecondaryStructureFilter(T=37, thr_DG=0)
        self.homodimer_filter = HomodimerFilter(max_len_selfcomp=6)

    def test_Tm_filter_default(self):
        seq_remove = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
        res = self.Tm_filter_default.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition with has been accepted! [MeltingTemperatureNNFilter]"

        seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
        res = self.Tm_filter_default.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [MeltingTemperatureNNFilter]"

    def test_Tm_filter_user(self):
        seq_remove = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
        res = self.Tm_filter_user.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition with user-defined parameters has been accepted! [MeltingTemperatureNNFilter]"

        seq_keep = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTAAAAA")
        res = self.Tm_filter_user.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions with user-defined parameters has not been accepted! [MeltingTemperatureNNFilter]"

    def tes_secondary_structure_filter(self):
        seq_remove = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
        res = self.secondary_structure_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition with has been accepted! [SecondaryStructureFilter]"

        seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
        res = self.secondary_structure_filter.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [SecondaryStructureFilter]"

    def test_homodimer_filter(self):
        seq_remove = Seq("TAACAATATATATTGTTA")
        res = self.homodimer_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition with has been accepted!"

        seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
        res = self.homodimer_filter.apply(seq_keep)
        assert res == True, f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted!"


class TestExperimentSpecificFilters(unittest.TestCase):
    """Test if melting temperature filters work, e.g. sequences having certain Tm or secondary structure probability should be removed."""

    def setUp(self):
        self.padlock_arms_filter = PadlockArmsFilter(
            arm_length_min=5,
            arm_Tm_dif_max=5,
            arm_Tm_min=40,
            arm_Tm_max=60,
            Tm_parameters=TM_PARAMETERS,
            Tm_salt_correction_parameters=TM_PARAMETERS_SALT_CORRECTION,
            Tm_chem_correction_parameters=TM_PARAMETERS_CHEM_CORRECTION,
        )

    def test_padlock_filter(self):
        seq_remove = Seq("TGTCGGATCTCTTCAACAAGCTGGTCAT")
        res = self.padlock_arms_filter.apply(seq_remove)
        assert (
            res == False
        ), f"error: A sequence ({seq_remove}) not fulfilling the condition with has been accepted! [PadlockArmsFilter]"

        seq_keep = Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
        res = self.padlock_arms_filter.apply(seq_keep)
        assert (
            res == True
        ), f"error: A sequence ({seq_keep}) fulfilling the conditions has not been accepted! [PadlockArmsFilter]"

    def test_detection_oligo_filter(self):
        # TODO: add test for this filter
        pass


class TestPropertyFilter(unittest.TestCase):
    """Test if melting temperature filters work, e.g. sequences having certain Tm or secondary structure probability should be removed."""

    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_property_filters")

        masked_sequences_filter = HardMaskedSequenceFilter(mask="N")
        GC_content_filter = GCContentFilter(GC_content_min=40, GC_content_max=60)
        GC_clamp_filter = GCClampFilter(n_bases=2, n_GC=1)

        Tm_filter = MeltingTemperatureNNFilter(
            Tm_min=52,
            Tm_max=67,
            Tm_parameters=TM_PARAMETERS,
            Tm_chem_correction_parameters=TM_PARAMETERS_CHEM_CORRECTION,
        )
        secondary_structure_filter = SecondaryStructureFilter(37, -5)

        padlock_arms_filter = PadlockArmsFilter(
            arm_length_min=5,
            arm_Tm_dif_max=5,
            arm_Tm_min=40,
            arm_Tm_max=60,
            Tm_parameters=TM_PARAMETERS,
            Tm_salt_correction_parameters=TM_PARAMETERS_SALT_CORRECTION,
            Tm_chem_correction_parameters=TM_PARAMETERS_CHEM_CORRECTION,
        )

        self.filters = [
            masked_sequences_filter,
            GC_content_filter,
            GC_clamp_filter,
            Tm_filter,
            secondary_structure_filter,
            padlock_arms_filter,
        ]

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_property_filter(self):
        os.makedirs(self.tmp_path, exist_ok=True)
        property_filter = PropertyFilter(filters=self.filters)

        seq_keep = Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
        res = property_filter._filter_sequence(seq_keep)
        assert res == True, f"error: A sequence ({seq_keep}) fulfilling all conditions has not been accepted!"

    def test_property_filter_on_database(self):
        property_filter = PropertyFilter(filters=self.filters)
        # check if apply function for property filter works
        oligo_sequence_generator = OligoSequenceGenerator(dir_output=self.tmp_path)

        file_fasta = oligo_sequence_generator.create_sequences_random(
            filename_out="random_sequences1",
            length_sequences=30,
            num_sequences=100,
            name_sequences="random_sequences1",
            base_alphabet_with_probability={"A": 0.1, "C": 0.3, "G": 0.4, "T": 0.2},
        )

        oligos = OligoDatabase(
            min_oligos_per_region=2, write_regions_with_insufficient_oligos=True, dir_output=self.tmp_path
        )
        oligos.load_database_from_fasta(
            files_fasta=file_fasta,
            sequence_type="oligo",
            region_ids=["random_sequences1"],
            database_overwrite=True,
        )

        property_filter.apply(sequence_type="oligo", oligo_database=oligos, n_jobs=2)
