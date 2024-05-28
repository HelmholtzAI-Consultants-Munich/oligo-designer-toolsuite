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
    """A filter designed to evaluate padlock probe arms for oligonucleotide sequences. It ensures that the arms
    meet specific criteria including minimum arm length, maximum temperature difference between arms, and individual
    arm melting temperatures within specified limits. Additionally, it can adjust for salt and chemical conditions
    affecting melting temperatures.

    :param arm_length_min: Minimum length for each arm of the padlock probe.
    :type arm_length_min: int
    :param arm_Tm_dif_max: Maximum allowed difference in melting temperature (Tm) between the two arms.
    :type arm_Tm_dif_max: float
    :param arm_Tm_min: Minimum acceptable melting temperature for each arm.
    :type arm_Tm_min: float
    :param arm_Tm_max: Maximum acceptable melting temperature for each arm.
    :type arm_Tm_max: float
    :param Tm_parameters: Parameters for the nearest-neighbor thermodynamic model to calculate Tm.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_salt_correction_parameters: Optional parameters for salt correction of Tm calculations.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Optional parameters for chemical correction of Tm calculations.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
    :type Tm_chem_correction_parameters: dict, optional
    """

    def __init__(
        self,
        arm_length_min: int,
        arm_Tm_dif_max: float,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ):
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

    def apply(self, sequence: Seq):
        """Applies the padlock arms filter to a DNA sequence.
        Applies the filter to evaluate if a given sequence is suitable for padlock probes based on
        arm length and melting temperature criteria.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
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
    """A filter designed to evaluate if detection oligos can be designed for oligonucleotides with ligation sites.
    It ensures that the detection oligo meets specific criteria including minimum and maximum length, and minimum
    number of Thymins in the sequence.

    This is dependent on the padlock arm criteria, i.e. it first evaluates if padlock probe arms can be designed
    and then uses the calculated ligation site as input for the detection oligo design.

    :param detect_oligo_length_min: The minimum length of the detection oligo.
    :type detect_oligo_length_min: int
    :param detect_oligo_length_max: The maximum length of the detection oligo.
    :type detect_oligo_length_max: int
    :param min_thymines: The minimum number of thymines in the detection oligo.
    :type min_thymines: int
    :param arm_length_min: Minimum length for each arm of the padlock probe.
    :type arm_length_min: int
    :param arm_Tm_dif_max: Maximum allowed difference in melting temperature (Tm) between the two arms.
    :type arm_Tm_dif_max: float
    :param arm_Tm_min: Minimum acceptable melting temperature for each arm.
    :type arm_Tm_min: float
    :param arm_Tm_max: Maximum acceptable melting temperature for each arm.
    :type arm_Tm_max: float
    :param Tm_parameters: Parameters for the nearest-neighbor thermodynamic model to calculate Tm.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_salt_correction_parameters: Optional parameters for salt correction of Tm calculations.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
    :type Tm_salt_correction_parameters: dict, optional
    :param Tm_chem_correction_parameters: Optional parameters for chemical correction of Tm calculations.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
    :type Tm_chem_correction_parameters: dict, optional
    """

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
    ):
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

    def apply(self, sequence: Seq):
        """Applies the detection oligo filter to a DNA sequence.
        Applies the filter to evaluate if a given sequence is suitable for designing detection oligos based on
        detection oligo length and number of Thymines criteria.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: True if the sequence meets the criteria, False otherwise.
        :rtype: bool
        """
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
