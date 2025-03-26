############################################
# imports
############################################

from abc import ABC, abstractmethod
from typing import Tuple

import pandas as pd

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoAttributes, OligoDatabase

############################################
# Oligo Scoring Classes
############################################


class OligoScoringBase(ABC):
    """
    Abstract base class for scoring oligonucleotides based on various criteria.

    This class provides a framework for developing specific scoring strategies for oligos.
    It defines the basic structure for scoring oligonucleotides by requiring the implementation
    of the `get_score` method in subclasses. The `apply` method is provided to score all oligos
    in a specified region and update the database with the calculated scores.
    """

    def apply(
        self, oligo_database: OligoDatabase, region_id: str, sequence_type: _TYPES_SEQ
    ) -> Tuple[OligoDatabase, pd.Series]:
        """
        Applies the scoring function to all oligos within a specified region of the OligoDatabase, updating the
        database with the calculated scores and returning the updated database along with a pandas Series of scores.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param sequence_type: The type of sequence to be used for score calculation.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: A tuple containing the updated OligoDatabase and a pandas Series with oligo scores.
        :rtype: Tuple[OligoDatabase, pd.Series]
        """
        oligos_ids = list(oligo_database.database[region_id].keys())
        oligos_scores = pd.Series(index=oligos_ids, dtype=float)
        for oligo_id in oligos_ids:
            score = round(self.get_score(oligo_database, region_id, oligo_id, sequence_type), 4)
            oligo_database.database[region_id][oligo_id]["oligo_score"] = score
            oligos_scores[oligo_id] = score
        return oligo_database, oligos_scores

    @abstractmethod
    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:
        """
        Abstract method to compute the score of a specific oligonucleotide sequence. This method should be
        implemented by subclasses to define the specific scoring logic.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param region_id: Oligo ID to process.
        :type region_id: str
        :param sequence_type: The type of sequence to be used for score calculation.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: The computed score of the nucleotide.
        :rtype: float
        """


class GCOligoScoring(OligoScoringBase):
    """
    A class for scoring oligonucleotides based on their GC content.

        :math:`score = |GC_{opt} - GC_{oligo}|`

    The `GCOligoScoring` class calculates the score for an oligo by evaluating the difference between
    its actual GC content and a predefined optimal GC content. The smaller the difference,
    the better the score. This scoring method helps in selecting oligos with GC content
    closer to the desired optimal value, which is critical for various applications such as
    PCR primer design and hybridization efficiency.

    :param GC_content_opt: The optimal GC content value used as the target for scoring.
    :type GC_content_opt: float
    """

    def __init__(
        self,
        GC_content_opt: float,
    ) -> None:
        """Constructor for the GCOligoScoring class."""
        self.GC_content_opt = GC_content_opt

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:
        """
        Calculates the GC content score for a specific nucleotide based on the optimal GC content.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param region_id: Oligo ID to process.
        :type region_id: str
        :param sequence_type: The type of sequence to be used for score calculation.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: The absolute difference between the optimal GC content and the GC content of the nucleotide sequence.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_attribute_value(
            attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        GC_content_oligo = OligoAttributes._calc_GC_content(sequence=sequence)
        return abs(self.GC_content_opt - GC_content_oligo)


class WeightedGCUtrScoring(OligoScoringBase):
    """
    A class used to score oligonucleotides based on their GC content and whether they originate from UTR regions.

        :math:`score = w_{GC}|GC_{opt} - GC_{oligo}| + w_{UTR}I_{UTR}`

    The `WeightedGCUtrScoring` class calculates the score of a nucleotide by evaluating its GC content
    in relation to a specified optimal value and by applying an additional weight if the sequence is from
    untranslated regions (UTRs). This allows for a more nuanced scoring based on both sequence composition
    and origin.

    :param GC_content_opt: The optimal GC content value used as the target for scoring.
    :type GC_content_opt: float
    :param GC_weight: The weight applied to the GC content difference in the scoring calculation.
    :type GC_weight: float
    :param UTR_weight: The weight applied if the sequence originates from a UTR.
    :type UTR_weight: float
    """

    def __init__(
        self,
        GC_content_opt: float,
        GC_weight: float,
        UTR_weight: float,
    ) -> None:
        """Constructor for the WeightedGCUtrScoring class."""
        self.GC_content_opt = GC_content_opt
        self.GC_weight = GC_weight
        self.UTR_weight = UTR_weight

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:
        """
        Calculates the weighted score for a given nucleotide based on its GC content and whether it originates
        from a UTR (Untranslated Region).

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param region_id: Oligo ID to process.
        :type region_id: str
        :param sequence_type: The type of sequence to be used for score calculation.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: The calculated score based on the weighted difference from optimal GC content and UTR consideration.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_attribute_value(
            attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        GC_content_oligo = OligoAttributes._calc_GC_content(sequence=sequence)

        regiontype = oligo_database.get_oligo_attribute_value(
            attribute="regiontype", region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        if regiontype:
            sequence_originates_from_UTR = "three_prime_UTR" in regiontype or "five_prime_UTR" in regiontype
        else:
            sequence_originates_from_UTR = False

        score = (
            self.GC_weight * abs(self.GC_content_opt - GC_content_oligo)
            + self.UTR_weight * sequence_originates_from_UTR
        )
        return score


class WeightedIsoformTmScoring(OligoScoringBase):

    def __init__(
        self,
        Tm_content_opt: float,
        Tm_weight: float,
        isoform_weight: float,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ) -> None:
        """Constructor for the WeightedIsoformTmScoring class."""
        self.Tm_content_opt = Tm_content_opt
        self.Tm_weight = Tm_weight
        self.isoform_weight = isoform_weight

        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:
        """
        Calculates the weighted score for a given nucleotide based on its GC content and whether it originates
        from a UTR (Untranslated Region).

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param region_id: Oligo ID to process.
        :type region_id: str
        :param sequence_type: The type of sequence to be used for score calculation.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: The calculated score based on the weighted difference from optimal GC content and UTR consideration.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_attribute_value(
            attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        transcript_id = oligo_database.get_oligo_attribute_value(
            attribute="transcript_id", region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        number_transcripts = oligo_database.get_oligo_attribute_value(
            attribute="number_transcripts", region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        if transcript_id and number_transcripts:
            # isoform consensus is given in % (0-100), hence we have to devide by 100
            # we use 1 - isoform consensus as score because here the lower the score the better
            isoform_consensus = 1 - (
                OligoAttributes._calc_isoform_consensus(transcript_id, number_transcripts) / 100
            )
        else:
            # if information on available, don't consider isoform consensus in scoring
            isoform_consensus = 0

        Tm_oligo = OligoAttributes._calc_TmNN(
            sequence=sequence,
            Tm_parameters=self.Tm_parameters,
            Tm_salt_correction_parameters=self.Tm_salt_correction_parameters,
            Tm_chem_correction_parameters=self.Tm_chem_correction_parameters,
        )

        score = self.Tm_weight * abs(self.Tm_content_opt - Tm_oligo) + self.isoform_weight * isoform_consensus
        return score


class WeightedTmGCOligoScoring(OligoScoringBase):
    """
    A class used to score oligonucleotides based on their melting temperature (Tm) and GC content.

        :math:`score = w_{Tm}\\left[I_{Tm_{oligo} \\ge Tm_{opt}}\\left(\\frac{|Tm_{oligo} - Tm_{opt}|}{Tm_{max} - Tm_{opt}}\\right) + I_{Tm_{oligo} < Tm_{opt}}\\left(\\frac{|Tm_{oligo} - Tm_{opt}|}{Tm_{opt} - Tm_{min}}\\right)\\right]`
        :math:`+ w_{GC}\\left[I_{GC_{oligo} \\ge GC_{opt}}\\frac{|GC_{oligo} - GC_{opt}|}{GC_{max} - GC_{opt}} + I_{GC_{oligo} < GC_{opt}}\\frac{|GC_{oligo} - GC_{opt}|}{GC_{opt} - GC_{min}}\\right]`

    The `WeightedTmGCOligoScoring` class evaluates nucleotides by calculating a weighted score that considers both
    the deviation of the oligo's melting temperature from an optimal value and the deviation of its GC content from a desired
    percentage. By combining these two factors with user-defined weights, this class provides a comprehensive score that
    reflects both the thermal stability and nucleotide composition of the oligo.

    :param Tm_min: The minimum acceptable melting temperature.
    :type Tm_min: float
    :param Tm_opt: The optimal melting temperature for scoring.
    :type Tm_opt: float
    :param Tm_max: The maximum acceptable melting temperature.
    :type Tm_max: float
    :param GC_content_min: The minimum acceptable GC content percentage.
    :type GC_content_min: float
    :param GC_content_opt: The optimal GC content percentage for scoring.
    :type GC_content_opt: float
    :param GC_content_max: The maximum acceptable GC content percentage.
    :type GC_content_max: float
    :param Tm_parameters: Parameters for calculating the melting temperature.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_weight: The weight assigned to the Tm scoring component.
    :type Tm_weight: float
    :param GC_weight: The weight assigned to the GC content scoring component.
    :type GC_weight: float
    :param Tm_salt_correction_parameters: Parameters for salt correction in Tm calculation (optional).
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Parameters for chemical correction in Tm calculation (optional).
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
    :type Tm_chem_correction_parameters: dict, optional
    """

    def __init__(
        self,
        Tm_min: float,
        Tm_opt: float,
        Tm_max: float,
        GC_content_min: float,
        GC_content_opt: float,
        GC_content_max: float,
        Tm_parameters: dict,
        Tm_weight: float,
        GC_weight: float,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ) -> None:
        """Constructor for the WeightedTmGCOligoScoring class."""
        self.Tm_min = Tm_min
        self.Tm_opt = Tm_opt
        self.Tm_max = Tm_max
        self.GC_content_min = GC_content_min
        self.GC_content_opt = GC_content_opt
        self.GC_content_max = GC_content_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters
        self.Tm_weight = Tm_weight
        self.GC_weight = GC_weight
        self._generate_scoring_functions()

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:
        """
        Computes the score of a nucleotide based on the weighted difference from optimal melting temperature (Tm)
        and GC content, using specified parameters for the calculation.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param region_id: Oligo ID to process.
        :type region_id: str
        :param sequence_type: The type of sequence to be used for score calculation.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: The calculated score based on the weighted difference from optimal GC content and Tm.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_attribute_value(
            attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        Tm_oligo = OligoAttributes._calc_TmNN(
            sequence=sequence,
            Tm_parameters=self.Tm_parameters,
            Tm_salt_correction_parameters=self.Tm_salt_correction_parameters,
            Tm_chem_correction_parameters=self.Tm_chem_correction_parameters,
        )
        GC_oligo = OligoAttributes._calc_GC_content(sequence=sequence)

        Tm_dif = Tm_oligo - self.Tm_opt
        GC_dif = GC_oligo - self.GC_content_opt
        score = self.Tm_weight * self.Tm_error(Tm_dif) + self.GC_weight * self.GC_error(GC_dif)
        return score

    def _generate_scoring_functions(self) -> None:
        """
        Generates scoring functions for Tm and GC content errors.

        This method creates error functions for melting temperature (Tm) and GC content by using the minimum, optimal, and maximum values provided.
        These error functions are later used to score oligonucleotides based on their deviation from the optimal values.

        :return: None
        """
        self.Tm_error = self._generate_error_function(self.Tm_min, self.Tm_opt, self.Tm_max)
        self.GC_error = self._generate_error_function(
            self.GC_content_min, self.GC_content_opt, self.GC_content_max
        )

    def _generate_error_function(self, min_val: float, opt_val: float, max_val: float):
        """
        Generates an error function for scoring based on deviation from optimal values.

        This method creates a scoring function that calculates the relative error for a given difference (`dif`) from an optimal value.
        The function scales the error based on the proximity of the difference to either the minimum or maximum thresholds,
        depending on whether the deviation is positive or negative.

        :param min_val: The minimum acceptable value.
        :type min_val: float
        :param opt_val: The optimal value.
        :type opt_val: float
        :param max_val: The maximum acceptable value.
        :type max_val: float
        :return: A function that calculates the relative error based on the input difference.
        :rtype: function
        """

        dif_max = abs(max_val - opt_val)
        dif_min = abs(opt_val - min_val)
        if dif_max == dif_min:
            return lambda dif: abs(dif) / dif_max
        else:
            return lambda dif: abs(dif) / (dif_max if dif > 0 else dif_min)


class WeightedIsoformTmGCOligoScoring(WeightedTmGCOligoScoring):
    """
    A class for scoring oligonucleotides based on melting temperature (Tm), GC content, and isoform consensus.

        :math:`score = w_{Tm}\\left[I_{Tm_{oligo} \\ge Tm_{opt}}\\left(\\frac{|Tm_{oligo} - Tm_{opt}|}{Tm_{max} - Tm_{opt}}\\right) + I_{Tm_{oligo} < Tm_{opt}}\\left(\\frac{|Tm_{oligo} - Tm_{opt}|}{Tm_{opt} - Tm_{min}}\\right)\\right]`
        :math:`+ w_{GC}\\left[I_{GC_{oligo} \\ge GC_{opt}}\\frac{|GC_{oligo} - GC_{opt}|}{GC_{max} - GC_{opt}} + I_{GC_{oligo} < GC_{opt}}\\frac{|GC_{oligo} - GC_{opt}|}{GC_{opt} - GC_{min}}\\right] + w_{IC} IC`

    This class extends `WeightedTmGCOligoScoring` by incorporating isoform targeting efficiency into the scoring criteria.
    It calculates a composite score using weights assigned to Tm, GC content, and isoform consensus, helping in the selection
    of oligos with optimal thermal stability, nucleotide composition, and isoform specificity.

    :param Tm_min: Minimum acceptable melting temperature.
    :type Tm_min: float
    :param Tm_opt: Optimal melting temperature.
    :type Tm_opt: float
    :param Tm_max: Maximum acceptable melting temperature.
    :type Tm_max: float
    :param GC_content_min: Minimum acceptable GC content.
    :type GC_content_min: float
    :param GC_content_opt: Optimal GC content.
    :type GC_content_opt: float
    :param GC_content_max: Maximum acceptable GC content.
    :type GC_content_max: float
    :param Tm_parameters: Parameters for calculating the melting temperature.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param isoform_weight: Weight given to the isoform consensus in the scoring.
    :type isoform_weight: float
    :param Tm_weight: Weight given to the melting temperature in the scoring.
    :type Tm_weight: float
    :param GC_weight: Weight given to the GC content in the scoring.
    :type GC_weight: float
    :param Tm_salt_correction_parameters: Parameters for salt correction in Tm calculation (optional).
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Parameters for chemical correction in Tm calculation (optional).
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
    :type Tm_chem_correction_parameters: dict, optional
    """

    def __init__(
        self,
        Tm_min: float,
        Tm_opt: float,
        Tm_max: float,
        GC_content_min: float,
        GC_content_opt: float,
        GC_content_max: float,
        Tm_parameters: dict,
        isoform_weight: float,
        Tm_weight: float,
        GC_weight: float,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ) -> None:
        """Constructor for the WeightedIsoformTmGCOligoScoring class."""
        self.Tm_min = Tm_min
        self.Tm_opt = Tm_opt
        self.Tm_max = Tm_max
        self.GC_content_min = GC_content_min
        self.GC_content_opt = GC_content_opt
        self.GC_content_max = GC_content_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters
        self.isoform_weight = isoform_weight
        self.Tm_weight = Tm_weight
        self.GC_weight = GC_weight
        self._generate_scoring_functions()

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:
        """
        Calculates the score of a nucleotide based on its melting temperature (Tm), GC content,
        and isoform consensus, using specified weights for each factor. The score is derived by combining
        the weighted errors in Tm and GC content with the isoform consensus, where a lower score is preferable.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param region_id: Oligo ID to process.
        :type region_id: str
        :param sequence_type: The type of sequence to be used for score calculation.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: The calculated score based on the Tm, GC content, and isoform consensus.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_attribute_value(
            attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        transcript_id = oligo_database.get_oligo_attribute_value(
            attribute="transcript_id", region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        number_transcripts = oligo_database.get_oligo_attribute_value(
            attribute="number_transcripts", region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        if transcript_id and number_transcripts:
            # isoform consensus is given in % (0-100), hence we have to devide by 100
            # we use 1 - isoform consensus as score because here the lower the score the better
            isoform_consensus = 1 - (
                OligoAttributes._calc_isoform_consensus(transcript_id, number_transcripts) / 100
            )
        else:
            # if information on available, don't consider isoform consensus in scoring
            isoform_consensus = 0

        Tm_oligo = OligoAttributes._calc_TmNN(
            sequence=sequence,
            Tm_parameters=self.Tm_parameters,
            Tm_salt_correction_parameters=self.Tm_salt_correction_parameters,
            Tm_chem_correction_parameters=self.Tm_chem_correction_parameters,
        )
        GC_oligo = OligoAttributes._calc_GC_content(sequence=sequence)

        Tm_dif = Tm_oligo - self.Tm_opt
        GC_dif = GC_oligo - self.GC_content_opt
        score = (
            self.Tm_weight * self.Tm_error(Tm_dif)
            + self.GC_weight * self.GC_error(GC_dif)
            + self.isoform_weight * isoform_consensus
        )
        return score


class WeightedIsoformTmGCOligoScoringTargetedExons(WeightedIsoformTmGCOligoScoring):
    r"""
    A class for scoring oligonucleotides based on melting temperature (Tm), GC content, isoform consensus, and targeted exons.

    This class extends `WeightedIsoformTmGCOligoScoring` by incorporating targeted exons into the scoring criteria.
    It adds an additional weight if the oligo targets an exon in `targeted_exons`.

    The score is calculated as:

    .. math::

        \text{score} =
        \begin{cases}
            score_\text{WeightedIsoformTmGCOligoScoring} + w_{\text{targeted_exons}}, & \text{if oligo is in a targeted exon} \\
            score_\text{WeightedIsoformTmGCOligoScoring}, & \text{otherwise}
        \end{cases}


    :param Tm_min: Minimum acceptable melting temperature.
    :type Tm_min: float
    :param Tm_opt: Optimal melting temperature.
    :type Tm_opt: float
    :param Tm_max: Maximum acceptable melting temperature.
    :type Tm_max: float
    :param GC_content_min: Minimum acceptable GC content.
    :type GC_content_min: float
    :param GC_content_opt: Optimal GC content.
    :type GC_content_opt: float
    :param GC_content_max: Maximum acceptable GC content.
    :type GC_content_max: float
    :param Tm_parameters: Parameters for calculating the melting temperature.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_salt_correction_parameters: Parameters for salt correction in Tm calculation (optional).
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
    :param targeted_exons: List of targeted exons.
    :type targeted_exons: list
    :param isoform_weight: Weight given to the isoform consensus in the scoring.
    :type isoform_weight: float
    :param Tm_weight: Weight given to the melting temperature in the scoring.
    :type Tm_weight: float
    :param GC_weight: Weight given to the GC content in the scoring.
    :type GC_weight: float
    :param targeted_exons_weight: Weight assigned to the targeted exons in the scoring (default: 2). Set this parameter to 4 to make the presence of the oligo in targeted exons critically important, 0 for not at all.
    :type targeted_exons_weight: float
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Parameters for chemical correction in Tm calculation (optional).
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
    :type Tm_chem_correction_parameters: dict, optional
    """

    def __init__(
        self,
        Tm_min: float,
        Tm_opt: float,
        Tm_max: float,
        GC_content_min: float,
        GC_content_opt: float,
        GC_content_max: float,
        Tm_parameters: dict,
        targeted_exons: list,
        isoform_weight: float,
        Tm_weight: float,
        GC_weight: float,
        targeted_exons_weight: float,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ):
        super().__init__(
            Tm_min,
            Tm_opt,
            Tm_max,
            GC_content_min,
            GC_content_opt,
            GC_content_max,
            Tm_parameters,
            isoform_weight,
            Tm_weight,
            GC_weight,
            Tm_salt_correction_parameters,
            Tm_chem_correction_parameters,
        )
        self.targeted_exons_weight = targeted_exons_weight
        self.targeted_exons = targeted_exons

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:
        exon_numbers = oligo_database.get_oligo_attribute_value(
            "exon_number", flatten=True, region_id=region_id, oligo_id=oligo_id
        )
        if exon_numbers is None:
            in_targeted_exons = False
        else:
            in_targeted_exons = any(item in self.targeted_exons for item in exon_numbers)
        return super().get_score(
            oligo_database, region_id, oligo_id, sequence_type
        ) + self.targeted_exons_weight * (not in_targeted_exons)
