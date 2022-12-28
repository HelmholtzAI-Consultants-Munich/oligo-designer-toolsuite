import filecmp
import shutil
import sys
import unittest

from Bio.SeqUtils import MeltingTemp as mt

sys.path.append("../oligo_designer_toolsuite")

from oligo_designer_toolsuite.IO import CustomOligoDB, CustomReferenceDB
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContent,
    MaskedSequences,
    MeltingTemperature,
    PadlockArms,
    PropertyFilter,
)
from oligo_designer_toolsuite.utils._data_parser import check_gtf_format


class TestDBGeneration(unittest.TestCase):
    """Tests that the oligos oligos_DB and reference oligos_DB generated are correct"""

    @classmethod
    def setUpClass(cls) -> None:
        """Define the classes and the filter parameters"""
        cls.genes = [
            "WASH7P",
            "DDX11L1",
            "TRNT",
            "NOC2L",
            "PLEKHN1",
            "AGRN",
            "UBE2J2",
            "DVL1",
            "MIB2",
            "LOC112268402_1",
        ]

        cls.Tm_parameters = {
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

        cls.Tm_correction_parameters = {
            "DMSO": 0,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
            "fmd": 20,
        }

        masked_sequences = MaskedSequences()
        GC_content = GCContent(GC_content_min=40, GC_content_max=60)
        melting_temperature = MeltingTemperature(
            Tm_min=52,
            Tm_max=67,
            Tm_parameters=cls.Tm_parameters,
            Tm_correction_parameters=cls.Tm_correction_parameters,
        )
        arms_tm = PadlockArms(
            min_arm_length=10,
            max_arm_Tm_dif=2,
            arm_Tm_min=38,
            arm_Tm_max=49,
            Tm_parameters=cls.Tm_parameters,
            Tm_correction_parameters=cls.Tm_correction_parameters,
        )

        filters = [masked_sequences, GC_content, melting_temperature, arms_tm]
        # cls.oligos_db = NcbiOligoDB(probe_length_min=30, probe_length_max=40, dir_output='tests/output')

        # If the anotation and fasta files are already saved on the machine, it is possible to direclty use them
        # instead of downloading them again.
        dir_annotation = "/home/francesco/Desktop/Work/NCBI"
        annotation = dir_annotation + "/GCF_000001405.40_GRCh38.p14_genomic.gtf"
        sequence = dir_annotation + "/GCF_000001405.40_GRCh38.p14_genomic.fna"
        cls.oligos_db = CustomOligoDB(
            probe_length_min=30,
            probe_length_max=40,
            file_annotation=annotation,
            file_sequence=sequence,
            species="unknown",
            genome_assembly="unknown",
            annotation_release="unknown",
            annotation_source="unknown",
            dir_output="tests/output",
        )
        # define the reference class
        cls.reference_db = CustomReferenceDB(
            file_annotation=cls.oligos_db.file_annotation,
            file_sequence=cls.oligos_db.file_sequence,
            species=cls.oligos_db.species,
            genome_assembly=cls.oligos_db.genome_assembly,
            annotation_release=cls.oligos_db.annotation_release,
            dir_output=cls.oligos_db.dir_output,
        )

        cls.property_filter = PropertyFilter(filters=filters)

    def test_transcriptome(self):
        """Test that the reference DB created is correct"""
        self.reference_db.create_reference_DB(
            dir_reference_DB="reference", block_size=self.oligos_db.probe_length_max + 5
        )
        file_computed = "tests/output/gene_trascript_computed.fna"
        # let's check only the first 200 000 lines
        with open(self.reference_db.file_reference_DB, mode="r") as complete:
            with open(file_computed, mode="w+") as truncated:
                for x in range(200000):
                    line = complete.readline()
                    truncated.write(line)
        file_correct = "tests/data/gene_trascript.fna"
        self.assertTrue(
            filecmp.cmp(file_computed, file_correct),
            "The gene transcript computed do not correspond to the correct one",
        )

    def test_list_sequences(self):
        """Test that the oligos oligos_DB created is correct"""

        def oligos_DB_to_list(oligos_DB):
            sequences = []
            for gene in oligos_DB.keys():
                for sequence in oligos_DB[gene].keys():
                    sequences.append(oligos_DB[gene][sequence]["probe_sequence"])
            sequences.sort()  # needed to compare
            return sequences

        def list_from_file(file):
            with open(file) as handle:
                lines = handle.readlines()
                sequences = [line.rstrip() for line in lines]
            return sequences

        self.oligos_db.create_oligos_DB(genes=self.genes)
        self.oligos_db = self.property_filter.apply(self.oligos_db)
        sequences_computed = oligos_DB_to_list(self.oligos_db.oligos_DB)
        sequences_correct = list_from_file("tests/data/sequences_10_genes.txt")
        sequences_correct.sort()
        self.assertListEqual(
            sequences_computed,
            sequences_correct,
            "The sequences computed do not correspond to the correct ones",
        )

    def test_read_write_oligos_DB_tsv(self):
        """Tetst that if write and read the oligos oligos_DB in tsv format, the oligos oligos_DB does not change."""
        self.oligos_db.create_oligos_DB(genes=[self.genes[0]])
        self.oligos_db.write_oligos_DB(format="tsv", dir_oligos_DB="prefiltered_probes")
        DB_correct = self.oligos_db.oligos_DB
        self.oligos_db.read_oligos_DB(
            format="tsv", file_oligos_DB_tsv=self.oligos_db.file_oligos_DB_tsv
        )  # overwrite the dict
        for probe_id in DB_correct[self.genes[0]].keys():
            self.assertDictEqual(
                DB_correct[self.genes[0]][probe_id],
                self.oligos_db.oligos_DB[self.genes[0]][probe_id],
                f"The oligos oligos_DB changes when it is written and read for {probe_id} in tsv fromat. \nThe correct dict is \n {DB_correct[self.genes[0]][probe_id]} \n and the computed one is \n {self.oligos_db.oligos_DB[self.genes[0]][probe_id]}",
            )

    def test_write_oligos_DB_gtf(self):
        """Test that the gtf file created is in the correct format"""
        self.oligos_db.create_oligos_DB(genes=[self.genes[0]])
        self.oligos_db.write_oligos_DB(format="gtf", dir_oligos_DB="prefiltered_probes")
        self.assertTrue(check_gtf_format(self.oligos_db.file_oligos_DB_gtf))

    def test_read_write_oligos_DB_gtf(self):
        """Test that if write and read the oligos oligos_DB in gtf format, the oligos oligos_DB does not change."""
        self.oligos_db.create_oligos_DB(genes=[self.genes[0]])
        self.oligos_db.write_oligos_DB(format="gtf", dir_oligos_DB="prefiltered_probes")
        DB_correct = self.oligos_db.oligos_DB
        self.oligos_db.read_oligos_DB(
            format="gtf",
            file_oligos_DB_gtf=self.oligos_db.file_oligos_DB_gtf,
            file_oligos_DB_fasta=self.oligos_db.file_oligos_DB_fasta,
        )  # overwrite the dict
        for probe_id in DB_correct[self.genes[0]].keys():
            self.assertDictEqual(
                DB_correct[self.genes[0]][probe_id],
                self.oligos_db.oligos_DB[self.genes[0]][probe_id],
                f"The oligos oligos_DB changes when it is written and read for {probe_id} in gtf format. \nThe correct dict is \n {DB_correct[self.genes[0]][probe_id]} \n and the computed one is \n {self.oligos_db.oligos_DB[self.genes[0]][probe_id]}",
            )

    @classmethod
    def tearDownClass(cls) -> None:
        """Delete all the files dowloaded."""
        shutil.rmtree("tests/output")
        del cls.oligos_db


# run the tests
# unittest.main()
