############################################
# imports
############################################

import os
import shutil
import unittest

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_calculator import (
    DetectOligoProperty,
    DGSecondaryStructureProperty,
    GCContentProperty,
    IsoformConsensusProperty,
    LengthComplementProperty,
    LengthProperty,
    LengthSelfComplementProperty,
    NumTargetedTranscriptsProperty,
    PadlockArmsProperty,
    PropertyCalculator,
    ReverseComplementSequenceProperty,
    SeedregionProperty,
    SeedregionSiteProperty,
    ShortenedSequenceProperty,
    TmNNProperty,
)

############################################
# setup
############################################

# Global Parameters
FILE_DATABASE_OLIGO_PROPERTIES = "tests/data/databases/database_oligo_properties.tsv"

############################################
# tests
############################################


class TestOligoProperties(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_properties")

        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="test_oligo_properties",
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_calculate_oligo_length(self) -> None:
        length_property = LengthProperty()
        calculator = PropertyCalculator(properties=[length_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        length1 = oligo_database.get_oligo_property_value(
            property="length_oligo", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        length2 = oligo_database.get_oligo_property_value(
            property="length_oligo", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        assert length1 == 20, "error: wrong oligo length"
        assert length2 == 29, "error: wrong oligo length"

    def test_calculate_reverse_complement_sequence(self) -> None:
        reverse_complement_property = ReverseComplementSequenceProperty(
            sequence_type_reverse_complement="oligo"
        )
        calculator = PropertyCalculator(properties=[reverse_complement_property])
        oligo_database = calculator.apply(
            oligo_database=self.oligo_database, sequence_type="target", n_jobs=1
        )

        target = oligo_database.get_oligo_property_value(
            property="target", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        oligo = oligo_database.get_oligo_property_value(
            property="oligo", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert target == "ATCGTCATCCATTGGGGCAT", "error: wrong target sequence"
        assert oligo == "ATGCCCCAATGGATGACGAT", "error: wrong oligo sequence"

    def test_calculate_shortened_sequence(self) -> None:
        # check if it works for oligos
        sequence_type = "oligo"

        shortened_sequence_property = ShortenedSequenceProperty(sequence_length=10, reverse=False)
        calculator = PropertyCalculator(properties=[shortened_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=self.oligo_database, sequence_type=sequence_type, n_jobs=1
        )

        sequence_short_oligo = oligo_database.get_oligo_property_value(
            property=f"{sequence_type}_short", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        shortened_sequence_property = ShortenedSequenceProperty(sequence_length=10, reverse=True)
        calculator = PropertyCalculator(properties=[shortened_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=self.oligo_database, sequence_type=sequence_type, n_jobs=1
        )

        sequence_short_oligo_reverse_read = oligo_database.get_oligo_property_value(
            property=f"{sequence_type}_short", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        # check if it works for targets
        sequence_type = "target"

        shortened_sequence_property = ShortenedSequenceProperty(sequence_length=5, reverse=False)
        calculator = PropertyCalculator(properties=[shortened_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=self.oligo_database, sequence_type=sequence_type, n_jobs=1
        )

        sequence_short_target = oligo_database.get_oligo_property_value(
            property=f"{sequence_type}_short", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        shortened_sequence_property = ShortenedSequenceProperty(sequence_length=5, reverse=True)
        calculator = PropertyCalculator(properties=[shortened_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=self.oligo_database, sequence_type=sequence_type, n_jobs=1
        )

        sequence_short_target_reverse_read = oligo_database.get_oligo_property_value(
            property=f"{sequence_type}_short", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )
        assert sequence_short_oligo == "GGCTAGGGAA", "error: wrong short sequence calculated"
        assert sequence_short_target == "CTCTA", "error: wrong short sequence calculated"
        assert sequence_short_oligo_reverse_read == "TCCAATAGAG", "error: wrong short sequence calculated"
        assert sequence_short_target_reverse_read == "TAGCC", "error: wrong short sequence calculated"

    def test_calculate_num_targeted_transcripts(self) -> None:
        num_targeted_transcripts_property = NumTargetedTranscriptsProperty()
        calculator = PropertyCalculator(properties=[num_targeted_transcripts_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        num_targeted_transcripts1 = oligo_database.get_oligo_property_value(
            property="num_targeted_transcripts", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        num_targeted_transcripts2 = oligo_database.get_oligo_property_value(
            property="num_targeted_transcripts", flatten=True, region_id="region_2", oligo_id="region_2::1"
        )
        num_targeted_transcripts3 = oligo_database.get_oligo_property_value(
            property="num_targeted_transcripts", flatten=True, region_id="region_3", oligo_id="region_3::4"
        )

        assert num_targeted_transcripts1 == 2, "error: wrong number targeted transcripts"
        assert num_targeted_transcripts2 == 1, "error: wrong number targeted transcripts"
        assert num_targeted_transcripts3 == 28, "error: wrong number targeted transcripts"

    def test_calculate_isoform_consensus(self) -> None:
        num_targeted_transcripts_property = NumTargetedTranscriptsProperty()
        calculator = PropertyCalculator(properties=[num_targeted_transcripts_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)
        isoform_consensus_property = IsoformConsensusProperty()
        calculator = PropertyCalculator(properties=[isoform_consensus_property])
        oligo_database = calculator.apply(oligo_database=oligo_database, sequence_type="oligo", n_jobs=1)

        isoform_consensus1 = oligo_database.get_oligo_property_value(
            property="isoform_consensus", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        isoform_consensus2 = oligo_database.get_oligo_property_value(
            property="isoform_consensus", flatten=True, region_id="region_2", oligo_id="region_2::1"
        )

        assert isoform_consensus1 == 100, "error: wrong isoform consensus, should be 100%"
        assert isoform_consensus2 == 50, "error: wrong isoform consensus, should be 50%"

    def test_calculate_seedregion(self) -> None:
        seedregion_property = SeedregionProperty(start=0.4, end=0.6)
        calculator = PropertyCalculator(properties=[seedregion_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        seedregion_start = oligo_database.get_oligo_property_value(
            property="seedregion_start", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        seedregion_end = oligo_database.get_oligo_property_value(
            property="seedregion_end", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert (seedregion_start == 8) and (seedregion_end == 12), "error: wrong seedregion calculated"

    def test_calculate_seedregion_site(self) -> None:
        seedregion_site_property = SeedregionSiteProperty(
            seedregion_size=5, seedregion_site_name="ligation_site"
        )
        calculator = PropertyCalculator(properties=[seedregion_site_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        seedregion_start = oligo_database.get_oligo_property_value(
            property="seedregion_start", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        seedregion_end = oligo_database.get_oligo_property_value(
            property="seedregion_end", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert (seedregion_start == 6) and (seedregion_end == 15), "error: wrong seedregion calculated"

    def test_calculate_GC_content(self) -> None:
        gc_content_property = GCContentProperty()
        calculator = PropertyCalculator(properties=[gc_content_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        GC_content = oligo_database.get_oligo_property_value(
            property="GC_content_oligo", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert GC_content == 50, "error: wrong GC content calculated"

    def test_calculate_TmNN(self) -> None:
        tm_nn_property = TmNNProperty(Tm_parameters={})
        calculator = PropertyCalculator(properties=[tm_nn_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        TmNN = oligo_database.get_oligo_property_value(
            property="TmNN_oligo", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert TmNN == 53.57, "error: wrong Tm calculated"

    def test_calculate_length_selfcomplement(self) -> None:
        length_selfcomplement_property = LengthSelfComplementProperty()
        calculator = PropertyCalculator(properties=[length_selfcomplement_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        length_selfcomplement = oligo_database.get_oligo_property_value(
            property="length_selfcomplement_oligo",
            flatten=True,
            region_id="region_3",
            oligo_id="region_3::3",
        )

        assert length_selfcomplement == 18, "error: wrong length of selfcomplement calculated"

    def test_calculate_length_complement(self) -> None:
        length_complement_property = LengthComplementProperty(comparison_sequence="AGTC")
        calculator = PropertyCalculator(properties=[length_complement_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        length_complement = oligo_database.get_oligo_property_value(
            property="length_complement_oligo_AGTC",
            flatten=True,
            region_id="region_1",
            oligo_id="region_1::3",
        )

        assert length_complement == 3, "error: wrong length of complement calculated"

    def test_calculate_secondary_structure_DG(self) -> None:
        dg_secondary_structure_property = DGSecondaryStructureProperty(T=37)
        calculator = PropertyCalculator(properties=[dg_secondary_structure_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        DG_secondary_structure = oligo_database.get_oligo_property_value(
            property="DG_secondary_structure_oligo",
            flatten=True,
            region_id="region_1",
            oligo_id="region_1::1",
        )

        assert DG_secondary_structure == 0.8, "error: wrong DG calculated"

    def test_calculate_padlock_arms(self) -> None:
        padlock_arms_property = PadlockArmsProperty(
            arm_length_min=3,
            arm_Tm_dif_max=15,
            arm_Tm_min=30,
            arm_Tm_max=80,
            Tm_parameters={},
        )
        calculator = PropertyCalculator(properties=[padlock_arms_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        ligation_site = oligo_database.get_oligo_property_value(
            property="ligation_site", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        assert ligation_site == 14, "error: wrong padlock arms calculated"

    def test_calculate_detect_oligo(self) -> None:
        detect_oligo_property = DetectOligoProperty(
            detect_oligo_length_min=8,
            detect_oligo_length_max=12,
            min_thymines=2,
        )
        calculator = PropertyCalculator(properties=[detect_oligo_property])
        oligo_database = calculator.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1)

        detect_oligo_even = oligo_database.get_oligo_property_value(
            property="detect_oligo_even", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )
        detect_oligo_long_left = oligo_database.get_oligo_property_value(
            property="detect_oligo_long_left", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )
        detect_oligo_long_right = oligo_database.get_oligo_property_value(
            property="detect_oligo_long_right", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        assert detect_oligo_even == "AGGGAATCGAAT", "error: wrong detection oligo even calculated"
        assert detect_oligo_long_left is None, "error: wrong detection oligo left calculated"
        assert detect_oligo_long_right is None, "error: wrong detection oligo right calculated"
