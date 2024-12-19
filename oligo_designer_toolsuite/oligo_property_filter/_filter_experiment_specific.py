############################################
# imports
############################################

from Bio.SeqUtils import Seq

from oligo_designer_toolsuite.database import OligoAttributes
from oligo_designer_toolsuite.oligo_property_filter import PropertyFilterBase

############################################
# Padlock Filter Classes
############################################


class PadlockArmsFilter(PropertyFilterBase):

    def __init__(
        self,
        arm_length_min: int,
        arm_Tm_dif_max: float,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ) -> None:
        """Constructor for the PadlockArmsFilter class."""
        super().__init__()
        if arm_Tm_max <= arm_Tm_min:
            raise ValueError("Tm_max is lower that Tm_min!")
        self.arm_length_min = arm_length_min
        self.arm_Tm_dif_max = arm_Tm_dif_max
        self.arm_Tm_min = arm_Tm_min
        self.arm_Tm_max = arm_Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters

    def apply(self, sequence: Seq) -> bool:
        _, _, ligation_site = OligoAttributes()._calc_padlock_arms(
            sequence,
            self.arm_length_min,
            self.arm_Tm_dif_max,
            self.arm_Tm_min,
            self.arm_Tm_max,
            self.Tm_parameters,
            self.Tm_salt_correction_parameters,
            self.Tm_chem_correction_parameters,
        )

        if ligation_site:
            return True
        else:
            return False


class DetectionOligoFilter(PropertyFilterBase):

    def __init__(
        self,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        min_thymines: int,
        arm_length_min: int,
        arm_Tm_dif_max: float,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ) -> None:
        """Constructor for the DetectionOligoFilter class."""
        super().__init__()
        if arm_Tm_max <= arm_Tm_min:
            raise ValueError("Tm_max is lower that Tm_min!")
        self.detect_oligo_length_min = detect_oligo_length_min
        self.detect_oligo_length_max = detect_oligo_length_max
        self.min_thymines = min_thymines
        self.arm_length_min = arm_length_min
        self.arm_Tm_dif_max = arm_Tm_dif_max
        self.arm_Tm_min = arm_Tm_min
        self.arm_Tm_max = arm_Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters

    def apply(self, sequence: Seq) -> bool:
        _, _, ligation_site = OligoAttributes()._calc_padlock_arms(
            sequence=sequence,
            arm_length_min=self.arm_length_min,
            arm_Tm_dif_max=self.arm_Tm_dif_max,
            arm_Tm_min=self.arm_Tm_min,
            arm_Tm_max=self.arm_Tm_max,
            Tm_parameters=self.Tm_parameters,
            Tm_salt_correction_parameters=self.Tm_salt_correction_parameters,
            Tm_chem_correction_parameters=self.Tm_chem_correction_parameters,
        )
        if not ligation_site:
            return False

        detect_oligo_even, _, _ = OligoAttributes()._calc_detect_oligo(
            sequence=sequence,
            ligation_site=ligation_site,
            detect_oligo_length_min=self.detect_oligo_length_min,
            detect_oligo_length_max=self.detect_oligo_length_max,
            min_thymines=self.min_thymines,
        )
        if detect_oligo_even:
            return True
        else:
            return False
