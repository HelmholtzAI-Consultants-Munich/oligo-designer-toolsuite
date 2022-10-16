import sys
import unittest

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

sys.path.append("../oligo_designer_toolsuite")

from oligo_designer_toolsuite.oligo_pre_filter._filter_base import (
    GCContent,
    MaskedSequences,
    MeltingTemperature,
)
from oligo_designer_toolsuite.oligo_pre_filter._filter_padlock_probes import PadlockArms
from oligo_designer_toolsuite.oligo_pre_filter.pre_filter import PreFilter


class TestPreFilters(unittest.TestCase):
    """Test that the filtering classes return the corrected output when given in input a sequence"""

    @classmethod
    def setUpClass(self):
        """Define the filter classes and their parameters"""
        self.Tm_parameters = {
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

        self.Tm_correction_parameters = {
            "DMSO": 0,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
            "fmd": 20,
        }

        self.masked_sequences = MaskedSequences()
        self.GC_content = GCContent(GC_content_min=40, GC_content_max=60)
        self.melting_temperature = MeltingTemperature(
            Tm_min=52,
            Tm_max=67,
            Tm_parameters=self.Tm_parameters,
            Tm_correction_parameters=self.Tm_correction_parameters,
        )
        self.arms_tm = PadlockArms(
            min_arm_length=10,
            max_Tm_dif=2,
            Tm_min=40,
            Tm_max=43,
            Tm_parameters=self.Tm_parameters,
            Tm_correction_parameters=self.Tm_correction_parameters,
        )

        filters = [
            self.masked_sequences,
            self.GC_content,
            self.melting_temperature,
            self.arms_tm,
        ]
        self.pre_filter = PreFilter(filters=filters)

    def tets_positive_outcome(self):
        """Tests that a correct sequences passes all the filers."""
        fulfills, _ = self.pre_filter.filter_sequence(
            Seq("TGTCGGATCTCTTCAACAAGCTGGTCATGA")
        )
        self.assertEqual(
            fulfills,
            True,
            "A sequence fulfilling the conditions has not been accepted!",
        )

    def test_negative_mask(self):
        "Tests that the masked sequences filter returns false"
        fulfills, _ = self.masked_sequences.apply(Seq("TGTCGGATCTCTTCAACANGCTGGTCATGA"))
        self.assertEqual(
            fulfills,
            False,
            "A sequence not fulfilling the mask condition has been accepted!",
        )

    def test_negative_GC(self):
        "Tests that the GC content filter returns false"
        fulfills, _ = self.GC_content.apply(Seq("TCGGGCGGGAGATCCAGGTGGCGCGCAAAG"))
        self.assertEqual(
            fulfills,
            False,
            "A sequence not fulfilling the GC condition has been accepted!",
        )

    def test_negative_tm(self):
        "Tests that the melting temperature filter returns false"
        fulfills, _ = self.melting_temperature.apply(
            Seq("TGGCTTGGGCCTTTCCAAGCCCCCATTTGAGCT")
        )
        self.assertEqual(
            fulfills,
            False,
            "A sequence not fulfilling the mnelting temperature condition has been accepted!",
        )

    def test_negative_arm(self):
        "Tests that the arms melting tempretature content filter returns false"
        fulfills, _ = self.arms_tm.apply(Seq("CTTGGGCCTTTCCAAGCCCCCATTTGAGCT"))
        self.assertEqual(
            fulfills,
            False,
            "A sequence not fulfilling the arm condition has been accepted!",
        )


# run the tets
# unittest.main()
