############################################
# imports
############################################


from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_calculator import BaseProperty
from oligo_designer_toolsuite.utils import cast_to_list

from ._property_functions import calc_isoform_consensus, calc_num_targeted_transcripts

############################################
# Region Property Classes
############################################


class NumTargetedTranscriptsProperty(BaseProperty):
    """
    A property class for calculating the number of unique transcripts targeted by an oligonucleotide.
    """

    def __init__(self) -> None:
        """Constructor for the NumTargetedTranscriptsProperty class."""
        super().__init__()

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the number of unique transcripts targeted by the oligonucleotide.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.  Note: This parameter is not used for this property.
        :type sequence_type: str
        :return: A dictionary containing the calculated number of targeted transcripts property.
        :rtype: dict
        """
        transcript_id = oligo_database.get_oligo_property_value(
            property="transcript_id", region_id=region_id, oligo_id=oligo_id, flatten=True
        )

        if transcript_id:
            num_targeted_transcripts = calc_num_targeted_transcripts(
                transcript_id=cast_to_list(transcript_id)
            )
        else:
            num_targeted_transcripts = None

        return {
            "num_targeted_transcripts": num_targeted_transcripts,
        }


class IsoformConsensusProperty(BaseProperty):
    """
    A property class for calculating the isoform consensus for an oligonucleotide.
    """

    def __init__(self) -> None:
        """Constructor for the IsoformConsensusProperty class."""
        super().__init__()

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the isoform consensus for the oligonucleotide.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.  Note: This parameter is not used for this property.
        :type sequence_type: str
        :return: A dictionary containing the calculated isoform consensus property.
        :rtype: dict
        """
        number_total_transcripts = oligo_database.get_oligo_property_value(
            property="number_total_transcripts", region_id=region_id, oligo_id=oligo_id, flatten=True
        )

        transcript_id = oligo_database.get_oligo_property_value(
            property="transcript_id", region_id=region_id, oligo_id=oligo_id, flatten=True
        )

        if transcript_id and number_total_transcripts:
            isoform_consensus = calc_isoform_consensus(
                transcript_id=cast_to_list(transcript_id),
                number_total_transcripts=cast_to_list(number_total_transcripts),
            )
        else:
            isoform_consensus = None

        return {
            "isoform_consensus": isoform_consensus,
        }
