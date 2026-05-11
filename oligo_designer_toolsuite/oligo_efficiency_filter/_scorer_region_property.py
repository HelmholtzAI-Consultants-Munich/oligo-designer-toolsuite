############################################
# imports
############################################

from typing import Any

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import BaseScorer
from oligo_designer_toolsuite.oligo_property_calculator._property_functions import calc_isoform_consensus
from oligo_designer_toolsuite.utils import cast_to_list

############################################
# Sequence Property Scorer Classes
############################################


class OverlapTargetedExonsScorer(BaseScorer):
    """
    Scores oligos based on whether they overlap with a set of targeted exons.

    If an oligo overlaps at least one of the specified targeted exons, it receives a weighted score.
    This can be used to prioritize probes targeting specific exonic regions.

    :param targeted_exons: List of exon identifiers that should be targeted by oligos.
    :type targeted_exons: list
    :param score_weight: Weight to apply if an oligo overlaps with a targeted exon.
    :type score_weight: float
    :param property_name: Name of the property to use for scoring.
    :type property_name: str, optional
    """

    def __init__(self, targeted_exons: list[str], score_weight: float, property_name: str = "exon_number"):
        """Constructor for the OverlapTargetedExonsScorer class."""

        self.targeted_exons = sorted(targeted_exons)
        self.score_weight = score_weight
        self.property_name = property_name

    def apply(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_id: str,
        sequence_type: str,
        **_: Any,
    ) -> float:
        """
        Apply the targeted exon overlap scoring strategy to a given oligo.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the score is computed.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.  Note: This parameter is not used in this function.
        :type sequence_type: str
        :param _: Additional keyword arguments accepted for extensibility. These are ignored by this method.
        :type _: Any
        :return: Weighted score based on overlap with targeted exons.
        :rtype: float
        """
        exon_numbers = oligo_database.get_oligo_property_value(
            self.property_name, flatten=True, region_id=region_id, oligo_id=oligo_id
        )
        exon_numbers = cast_to_list(exon_numbers) if exon_numbers else None

        if exon_numbers is None:
            # if no exon numbers are available, consider the oligo to be in a targeted exon
            in_targeted_exons = True
        else:
            in_targeted_exons = any(
                any(exon in self.targeted_exons for exon in str(exon_number).split("__JUNC__"))
                for exon_number in exon_numbers
            )

        score = self.score_weight * (1 - in_targeted_exons)

        return score


class OverlapUTRScorer(BaseScorer):
    """
    Scores oligos based on whether they originate from untranslated regions (UTRs).

    Oligos that map to either 5' or 3' UTRs receive a weighted score. This is useful when
    targeting regulatory regions of transcripts.

    :param score_weight: Weight to apply if an oligo originates from a UTR.
    :type score_weight: float
    :param property_name: Name of the property to use for scoring the region type.
    :type property_name: str, optional
    """

    def __init__(self, score_weight: float, property_name: str = "regiontype") -> None:
        """Constructor for the OverlapUTRScorer class."""

        self.score_weight = score_weight
        self.property_name = property_name

    def apply(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_id: str,
        sequence_type: str,
        **_: Any,
    ) -> float:
        """
        Apply the UTR overlap scoring strategy to a given oligo.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the score is computed.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.  Note: This parameter is not used in this function.
        :type sequence_type: str
        :param _: Additional keyword arguments accepted for extensibility. These are ignored by this method.
        :type _: Any
        :return: Weighted score based on UTR overlap.
        :rtype: float
        """
        regiontype = oligo_database.get_oligo_property_value(
            property=self.property_name, region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        if regiontype:
            sequence_originates_from_UTR = "three_prime_UTR" in regiontype or "five_prime_UTR" in regiontype
        else:
            # if no regiontype is available, consider the oligo to originate from a UTR
            sequence_originates_from_UTR = True

        score = self.score_weight * (1 - sequence_originates_from_UTR)

        return score


class IsoformConsensusScorer(BaseScorer):
    """
    Scores oligos based on their presence across transcript isoforms (isoform consensus).

    The consensus (in percent) is normalized to a 0–1 range and inverted so that
    higher consensus yields lower scores.

    :param score_weight: Weight to apply to the consensus score.
    :type score_weight: float
    :param property_name_transcript_id: Name of the property to use for scoring the transcript ID.
    :type property_name_transcript_id: str, optional
    :param property_name_number_total_transcripts: Name of the property to use for scoring the number of total transcripts.
    :type property_name_number_total_transcripts: str, optional
    """

    def __init__(
        self,
        score_weight: float,
        property_name_transcript_id: str = "transcript_id",
        property_name_number_total_transcripts: str = "number_total_transcripts",
    ) -> None:
        """Constructor for the IsoformConsensusScorer class."""
        self.score_weight = score_weight
        self.property_name_transcript_id = property_name_transcript_id
        self.property_name_number_total_transcripts = property_name_number_total_transcripts

    def apply(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_id: str,
        sequence_type: str,
        **_: Any,
    ) -> float:
        """
        Apply the isoform consensus scoring strategy to a given oligo.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the score is computed.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.  Note: This parameter is not used in this function.
        :type sequence_type: str
        :param _: Additional keyword arguments accepted for extensibility. These are ignored by this method.
        :type _: Any
        :return: Weighted score based on isoform consensus.
        :rtype: float
        """
        transcript_id = oligo_database.get_oligo_property_value(
            property=self.property_name_transcript_id, region_id=region_id, oligo_id=oligo_id, flatten=True
        )

        number_transcripts = oligo_database.get_oligo_property_value(
            property=self.property_name_number_total_transcripts,
            region_id=region_id,
            oligo_id=oligo_id,
            flatten=True,
        )

        if transcript_id and number_transcripts:
            isoform_consensus = calc_isoform_consensus(
                transcript_id=cast_to_list(transcript_id),
                number_total_transcripts=cast_to_list(number_transcripts),
            )
            # isoform consensus is given in % (0–100); use 1 - consensus/100 so higher
            # consensus leads to a lower score (better oligo).
            isoform_consensus = 1 - (isoform_consensus / 100)
        else:
            # if information not available, don't consider isoform consensus in scoring
            isoform_consensus = 0

        score = self.score_weight * isoform_consensus

        return score
