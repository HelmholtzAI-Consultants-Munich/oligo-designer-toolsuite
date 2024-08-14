############################################
# imports
############################################

from abc import ABC, abstractmethod

import pandas as pd

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoAttributes, OligoDatabase

############################################
# Oligo Scoring Classes
############################################


class OligoScoringBase(ABC):

    def apply(
        self, oligo_database: OligoDatabase, region_id: str, sequence_type: _TYPES_SEQ
    ) -> tuple[OligoDatabase, pd.Series]:

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
        """ """


class GCOligoScoring(OligoScoringBase):

    def __init__(
        self,
        GC_content_opt: float,
    ) -> None:
        """Constructor for the GCOligoScoring class."""
        self.GC_content_opt = GC_content_opt

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:
        sequence = oligo_database.get_oligo_attribute_value(
            attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )
        GC_content_oligo = OligoAttributes._calc_GC_content(sequence=sequence)
        return abs(self.GC_content_opt - GC_content_oligo)


class WeightedGCUtrScoring(OligoScoringBase):

    def __init__(
        self,
        GC_content_opt: float,
        GC_weight: float = 1,
        UTR_weight: float = 10,
    ) -> None:
        """Constructor for the WeightedGCUtrScoring class."""
        self.GC_content_opt = GC_content_opt
        self.GC_weight = GC_weight
        self.UTR_weight = UTR_weight

    def get_score(
        self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: _TYPES_SEQ
    ) -> float:

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
        self.Tm_error = self._generate_error_function(self.Tm_min, self.Tm_opt, self.Tm_max)
        self.GC_error = self._generate_error_function(
            self.GC_content_min, self.GC_content_opt, self.GC_content_max
        )

    def _generate_error_function(self, min_val: float, opt_val: float, max_val: float):

        dif_max = max_val - opt_val
        dif_min = opt_val - min_val
        if dif_max == dif_min:
            return lambda dif: abs(dif) / dif_max
        else:
            return lambda dif: abs(dif) / (dif_max if dif > 0 else dif_min)


class WeightedIsoformTmGCOligoScoring(WeightedTmGCOligoScoring):

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
