############################################
# imports
############################################

from abc import ABC, abstractmethod

import pandas as pd

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoAttributes, OligoDatabase
from oligo_designer_toolsuite.utils import flatten_attribute_list

############################################
# Oligo Scoring Classes
############################################


class OligoScoringBase(ABC):
    """Abstract base class for scoring oligonucleotides based on specified criteria.
    It provides a framework for scoring oligos, with an abstract method `get_score` that must be implemented in subclasses.
    """

    def apply(self, oligo_database: OligoDatabase, region_id: str, sequence_type: _TYPES_SEQ):
        """Scores all oligonucleotides in the provided dictionary using a defined scoring function,
        updating the dictionary with scores and also returning scores as a pandas.Series for efficient data manipulation.

        :param oligo_database: OligoDatabase containing the oligonucleotides with their respective information (e.g. oligo sequence and oligo attributes).
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequences being scored, which must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :return: Tuple containing the updated dictionary of oligonucleotides and a pandas.Series with the computed scores.
        :rtype: (dict, pandas.Series)
        """
        """Scores all oligonucleotides in the provided database using a defined scoring function,
        updating the database with scores and also returning scores as a pandas.Series for efficient data manipulation.

        :param oligo_database: OligoDatabase containing the oligonucleotides with their respective information (e.g. oligo sequence and oligo attributes).
        :type oligo_database: OligoDatabase
        :param region_id: The identifier for the specific region to score oligonucleotides.
        :type region_id: str
        :param sequence_type: The type of sequences being scored, which must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :return: Tuple containing the updated dictionary of oligonucleotides and a pandas.Series with the computed scores.
        :rtype: (dict, pandas.Series)
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
    ):
        """Abstract method to compute the score of a given oligonucleotide sequence.
        This method needs to be implemented by all subclasses to define specific scoring logic.

        :param oligo_attributes: A dictionary containing attributes of the oligo.
        :type oligo_attributes: dict
        :param sequence_type: The type of sequence (e.g., 'oligo', 'target').
        :type sequence_type: _TYPES_SEQ
        :return: Computed score of the oligonucleotide.
        :rtype: float
        """


class GCOligoScoring(OligoScoringBase):
    """Scoring class that calculates oligo scores based on their GC content deviation from the optimal GC content.

    $score = |GC_{opt} - GC_{oligo}|$.

    :param GC_content_opt: Optimal GC content.
    :type GC_content_opt: float
    """

    def __init__(
        self,
        GC_content_opt: float,
    ):
        """Constructor for the GCOligoScoring class."""
        self.GC_content_opt = GC_content_opt

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ):
        """Calculates the GC content score for a given oligonucleotide sequence based
        on its deviation from the optimal GC content.
        Score: the lower the better.

        :param oligo_attributes: A dictionary containing attributes of the oligo.
        :type oligo_attributes: dict
        :param sequence_type: The type of sequence (e.g., 'oligo', 'target').
        :type sequence_type: _TYPES_SEQ
        :return: The calculated score based on Tm and GC content.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_attribute_value(
            attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        GC_content_oligo = OligoAttributes._calc_GC_content(sequence=sequence)
        return abs(self.GC_content_opt - GC_content_oligo)


class WeightedGCUtrScoring(OligoScoringBase):
    """Scoring class that evaluates oligos based on their GC content and whether they originate from untranslated regions (UTRs).

    This class assigns scores to oligos by considering their GC content and whether they are part of the 3' or 5' UTRs.
    The score is computed as a weighted sum of the difference from optimal GC content and the presence in UTR regions, with different weights assigned to each factor.

    $score = w_{GC} * |GC_{opt} - GC_{oligo}| + w_{UTR} * I_{UTR}$.

    :param GC_content_opt: Optimal percentage of guanine and cytosine.
    :type GC_content_opt: float
    :param GC_weight: Weight for the GC content component in the scoring, defaults to 1.
    :type GC_weight: float
    :param UTR_weight: Weight for the UTR component in the scoring, defaults to 10.
    :type UTR_weight: float
    """

    def __init__(
        self,
        GC_content_opt: float,
        GC_weight: float = 1,
        UTR_weight: float = 10,
    ):
        """Constructor for the WeightedGCUtrScoring class."""
        self.GC_content_opt = GC_content_opt
        self.GC_weight = GC_weight
        self.UTR_weight = UTR_weight

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ):
        """Calculates the GC content score for a given oligonucleotide sequence based
        on its deviation from the optimal GC content and whether it originates from the UTR.
        Score: the lower the better.

        :param oligo_attributes: Dictionary containing attributes of the oligo.
        :type oligo_attributes: dict
        :param sequence_type: The type of sequence being scored.
        :type sequence_type: _TYPES_SEQ
        :return: Calculated score for the oligo.
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


class WeightedTmGCOligoScoring(OligoScoringBase):
    """Scoring class that calculates oligo scores based on their melting temperature (Tm) and GC content (GC).

    $score =
    w_{Tm}[I_{Tm_{oligo} \ge Tm_{opt}}(\frac{|Tm_{oligo} - Tm_{opt}|}{Tm_{max} - Tm_{opt}}) + I_{Tm_{oligo} < Tm_{opt}}(\frac{|Tm_{oligo} - Tm_{opt}|}{Tm_{opt} - Tm_{min}})] +
    w_{GC}[I_{GC_{oligo} \ge GC_{opt}}\frac{|GC_{oligo} - GC_{opt}|}{GC_{max} - GC_{opt}} + I_{GC_{oligo} < GC_{opt}}\frac{|GC_{oligo} - GC_{opt}|}{GC_{opt} - GC_{min}}]$.

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
    :param Tm_parameters: Parameters for calculating melting temperature.
    :type Tm_parameters: dict
    :param Tm_salt_correction_parameters: Parameters for salt correction in Tm calculation, optional.
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Parameters for chemical correction in Tm calculation, optional.
    :type Tm_chem_correction_parameters: dict, optional
    :param Tm_weight: Weight factor for Tm deviations in score calculation.
    :type Tm_weight: float
    :param GC_weight: Weight factor for GC content deviations in score calculation.
    :type GC_weight: float
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
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
        Tm_weight: float = 1,
        GC_weight: float = 1,
    ):
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
    ):
        """Calculates the oligo score based on Tm and GC content deviations from their optimal values.
        An error term is added to the GC and Tm scores to normalize both properties and make them comparable for the scoring.
        Score: the lower the better.

        :param oligo_attributes: A dictionary containing attributes of the oligo.
        :type oligo_attributes: dict
        :param sequence_type: The type of sequence (e.g., 'oligo', 'target').
        :type sequence_type: _TYPES_SEQ
        :return: The calculated score based on Tm and GC content.
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

    def _generate_scoring_functions(self):
        """Sets up dynamic functions to calculate errors in Tm and GC content based on the deviation from optimal values.
        This method is intended for internal use to configure the scoring calculations upon initialization.
        """
        self.Tm_error = self._generate_error_function(self.Tm_min, self.Tm_opt, self.Tm_max)
        self.GC_error = self._generate_error_function(
            self.GC_content_min, self.GC_content_opt, self.GC_content_max
        )

    def _generate_error_function(self, min_val, opt_val, max_val):
        """Generates scoring functions for Tm and GC content based on provided parameters.

        :param min_val: The minimum value.
        :type min_val: float
        :param opt_val: The optimal value.
        :type opt_val: float
        :param max_val: The maximum value.
        :type max_val: float
        """
        dif_max = max_val - opt_val
        dif_min = opt_val - min_val
        if dif_max == dif_min:
            return lambda dif: abs(dif) / dif_max
        else:
            return lambda dif: abs(dif) / (dif_max if dif > 0 else dif_min)


class WeightedIsoformTmGCOligoScoring(WeightedTmGCOligoScoring):
    """Scoring class that calculates oligo scores based on their melting temperature (Tm), GC content (GC) and isoform consensus (IC).
    The isoform consensus indicated the percentage of the total number of isoforms of the genomic region that are covered by the oligo.

    $score =
    w_{Tm}[I_{Tm_{oligo} \ge Tm_{opt}}(\frac{|Tm_{oligo} - Tm_{opt}|}{Tm_{max} - Tm_{opt}}) + I_{Tm_{oligo} < Tm_{opt}}(\frac{|Tm_{oligo} - Tm_{opt}|}{Tm_{opt} - Tm_{min}})] +
    w_{GC}[I_{GC_{oligo} \ge GC_{opt}}\frac{|GC_{oligo} - GC_{opt}|}{GC_{max} - GC_{opt}} + I_{GC_{oligo} < GC_{opt}}\frac{|GC_{oligo} - GC_{opt}|}{GC_{opt} - GC_{min}}] +
    w_{IC}IC$.

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
    :param Tm_parameters: Parameters for calculating melting temperature.
    :type Tm_parameters: dict
    :param Tm_salt_correction_parameters: Parameters for salt correction in Tm calculation, optional.
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Parameters for chemical correction in Tm calculation, optional.
    :type Tm_chem_correction_parameters: dict, optional
    :param isoform_weight: Weight assigned to isoform consensus in scoring, defaults to 2.
    :type isoform_weight: float, optional
    :param Tm_weight: Weight factor for Tm deviations in score calculation.
    :type Tm_weight: float, optional
    :param GC_weight: Weight factor for GC content deviations in score calculation.
    :type GC_weight: float, optional
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
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
        isoform_weight: float = 2,
        Tm_weight: float = 1,
        GC_weight: float = 1,
    ):
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
    ):
        """Calculates the oligo score based on Tm and GC content deviations from their optimal values and the isoform consensus.
        An error term is added to the GC and Tm scores to normalize both properties and make them comparable for the scoring.
        Score: the lower the better.

        :param oligo_attributes: A dictionary containing attributes of the oligo.
        :type oligo_attributes: dict
        :param sequence_type: The type of sequence (e.g., 'oligo', 'target').
        :type sequence_type: _TYPES_SEQ
        :return: The calculated score based on Tm and GC content.
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
