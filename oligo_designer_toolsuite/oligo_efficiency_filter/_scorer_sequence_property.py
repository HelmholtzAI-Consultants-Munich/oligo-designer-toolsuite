############################################
# imports
############################################

from typing import Any

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import BaseScorer
from oligo_designer_toolsuite.oligo_property_calculator import calc_gc_content, calc_tm_nn

############################################
# Sequence Property Scorer Classes
############################################


class SequencePropertyScorer(BaseScorer):
    """
    Abstract scorer class for evaluating oligonucleotide properties based on deviations
    from optimal or desired sequence characteristics (e.g., GC content, melting temperature).
    """

    def __init__(self) -> None:
        """Constructor for the SequencePropertyScorer class."""

    def _normalize_deviation(self, val_dev: float, val_min: float, val_opt: float, val_max: float) -> float:
        """
        Normalize the deviation of a value from its optimum based on asymmetrical min/max ranges.

        :param val_dev: The raw deviation from the optimal value.
        :type val_dev: float
        :param val_min: The minimum acceptable value in the range.
        :type val_min: float
        :param val_opt: The optimal (target) value.
        :type val_opt: float
        :param val_max: The maximum acceptable value in the range.
        :type val_max: float
        :return: Normalized deviation in the range [0, 1].
        :rtype: float
        """
        dev_max = abs(val_max - val_opt)
        dev_min = abs(val_opt - val_min)

        if dev_max == dev_min:
            val_dev_norm = abs(val_dev) / dev_max

        else:
            val_dev_norm = abs(val_dev) / (dev_max if val_dev > 0 else dev_min)
        return val_dev_norm


class DeviationFromOptimalGCContentScorer(SequencePropertyScorer):
    """
    Scores oligonucleotides based on their absolute deviation from a defined optimal GC content.

    :param GC_content_opt: The optimal GC content to compare against.
    :type GC_content_opt: float
    :param score_weight: Weight applied to the GC content deviation score.
    :type score_weight: float
    """

    def __init__(self, GC_content_opt: float, score_weight: float) -> None:
        """Constructor for the DeviationFromOptimalGCContentScorer class."""

        self.GC_content_opt = GC_content_opt
        self.score_weight = score_weight

    def apply(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_id: str,
        sequence_type: str,
        **_: Any,
    ) -> float:
        """
        Calculate a score based on the absolute deviation of GC content from the optimal value.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the score is computed.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param kwargs: Additional keyword arguments.
        :type kwargs: dict[str, Any]
        :return: Weighted score based on GC content deviation.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_property_value(
            property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )

        if not isinstance(sequence, str):
            raise TypeError(f"Sequence must be a string, got {type(sequence)}")

        GC_content_oligo = calc_gc_content(sequence=sequence)
        GC_content_dev = abs(self.GC_content_opt - GC_content_oligo)

        score = self.score_weight * GC_content_dev

        return score


class DeviationFromOptimalTmScorer(SequencePropertyScorer):
    """
    Scores oligos based on the absolute deviation from an optimal melting temperature (Tm).

    :param Tm_opt: Optimal melting temperature for scoring.
    :type Tm_opt: float
    :param Tm_parameters: Parameters for Tm calculation using the nearest-neighbor model.
    :type Tm_parameters: dict
    :param score_weight: Weight applied to the Tm deviation score.
    :type score_weight: float
    :param Tm_salt_correction_parameters: Parameters for salt correction.
    :type Tm_salt_correction_parameters: dict
    :param Tm_chem_correction_parameters: Parameters for chemical correction.
    :type Tm_chem_correction_parameters: dict
    """

    def __init__(
        self,
        Tm_opt: float,
        Tm_parameters: dict,
        score_weight: float,
        Tm_salt_correction_parameters: dict | None = None,
        Tm_chem_correction_parameters: dict | None = None,
    ) -> None:
        """Constructor for the DeviationFromOptimalTmScorer class."""

        self.Tm_opt = Tm_opt
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters
        self.score_weight = score_weight

    def apply(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_id: str,
        sequence_type: str,
        **_: Any,
    ) -> float:
        """
        Calculate a score based on the absolute deviation of Tm from the optimal value.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the score is computed.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param kwargs: Additional keyword arguments.
        :type kwargs: dict[str, Any]
        :return: Weighted score based on Tm deviation.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_property_value(
            property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )

        if not isinstance(sequence, str):
            raise TypeError(f"Sequence must be a string, got {type(sequence)}")

        Tm_oligo = calc_tm_nn(
            sequence=sequence,
            Tm_parameters=self.Tm_parameters,
            Tm_salt_correction_parameters=self.Tm_salt_correction_parameters,
            Tm_chem_correction_parameters=self.Tm_chem_correction_parameters,
        )
        Tm_dev = abs(self.Tm_opt - Tm_oligo)

        score = self.score_weight * Tm_dev

        return score


class NormalizedDeviationFromOptimalGCContentScorer(SequencePropertyScorer):
    """
    Scores oligos based on the normalized deviation of their GC content from a specified optimal value,
    taking into account asymmetric acceptable GC content ranges.

    This scoring approach allows for non-linear penalization of deviations, depending on whether the GC content is above or
    below the optimum. The deviation is normalized using the distance between the optimal GC content and the provided minimum
    or maximum bounds, enabling flexible weighting of upstream or downstream deviations.

    :param GC_content_min: Minimum acceptable GC content.
    :type GC_content_min: float
    :param GC_content_opt: Optimal GC content.
    :type GC_content_opt: float
    :param GC_content_max: Maximum acceptable GC content.
    :type GC_content_max: float
    :param score_weight: Weight applied to the normalized deviation score.
    :type score_weight: float
    """

    def __init__(
        self,
        GC_content_min: float,
        GC_content_opt: float,
        GC_content_max: float,
        score_weight: float,
    ) -> None:
        """Constructor for the NormalizedDeviationFromOptimalGCContentScorer class."""

        self.GC_content_min = GC_content_min
        self.GC_content_opt = GC_content_opt
        self.GC_content_max = GC_content_max
        self.score_weight = score_weight

    def apply(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_id: str,
        sequence_type: str,
        **_: Any,
    ) -> float:
        """
        Calculate a score based on normalized deviation of GC content from optimal value.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the score is computed.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param kwargs: Additional keyword arguments.
        :type kwargs: dict[str, Any]
        :return: Weighted score based on normalized GC content deviation.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_property_value(
            property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )

        if not isinstance(sequence, str):
            raise TypeError(f"Sequence must be a string, got {type(sequence)}")

        GC_content_oligo = calc_gc_content(sequence=sequence)
        GC_content_dev = GC_content_oligo - self.GC_content_opt
        GC_content_dev_norm = self._normalize_deviation(
            GC_content_dev, self.GC_content_min, self.GC_content_opt, self.GC_content_max
        )

        score = self.score_weight * GC_content_dev_norm

        return score


class NormalizedDeviationFromOptimalTmScorer(SequencePropertyScorer):
    """
    Scores oligos based on the normalized deviation of their melting temperature (Tm) from a specified optimal value,
    taking into account asymmetric acceptable Tm ranges.

    This scoring approach allows for non-linear penalization of deviations, depending on whether the Tm is above or
    below the optimum. The deviation is normalized using the distance between the optimal Tm and the provided minimum
    or maximum bounds, enabling flexible weighting of upstream or downstream deviations.

    :param Tm_min: Minimum acceptable Tm.
    :type Tm_min: float
    :param Tm_opt: Optimal Tm.
    :type Tm_opt: float
    :param Tm_max: Maximum acceptable Tm.
    :type Tm_max: float
    :param Tm_parameters: Parameters for Tm calculation using the nearest-neighbor model.
    :type Tm_parameters: dict
    :param score_weight: Weight applied to the normalized deviation score.
    :type score_weight: float
    :param Tm_salt_correction_parameters: Salt correction parameters.
    :type Tm_salt_correction_parameters: dict
    :param Tm_chem_correction_parameters: Chemical correction parameters.
    :type Tm_chem_correction_parameters: dict

    """

    def __init__(
        self,
        Tm_min: float,
        Tm_opt: float,
        Tm_max: float,
        Tm_parameters: dict,
        score_weight: float,
        Tm_salt_correction_parameters: dict | None = None,
        Tm_chem_correction_parameters: dict | None = None,
    ) -> None:
        """Constructor for the NormalizedDeviationFromOptimalTmScorer class."""

        self.Tm_min = Tm_min
        self.Tm_opt = Tm_opt
        self.Tm_max = Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters
        self.score_weight = score_weight

    def apply(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_id: str,
        sequence_type: str,
        **_: Any,
    ) -> float:
        """
        Calculate a score based on normalized deviation of Tm from the optimal value.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the score is computed.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param kwargs: Additional keyword arguments.
        :type kwargs: dict[str, Any]
        :return: Weighted score based on normalized Tm deviation.
        :rtype: float
        """
        sequence = oligo_database.get_oligo_property_value(
            property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
        )

        if not isinstance(sequence, str):
            raise TypeError(f"Sequence must be a string, got {type(sequence)}")

        Tm_oligo = calc_tm_nn(
            sequence=sequence,
            Tm_parameters=self.Tm_parameters,
            Tm_salt_correction_parameters=self.Tm_salt_correction_parameters,
            Tm_chem_correction_parameters=self.Tm_chem_correction_parameters,
        )
        Tm_dev = Tm_oligo - self.Tm_opt
        Tm_dev_norm = self._normalize_deviation(Tm_dev, self.Tm_min, self.Tm_opt, self.Tm_max)

        score = self.score_weight * Tm_dev_norm

        return score
