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
    """
    A filter class for evaluating the suitability of padlock probe arms based on specific thermodynamic criteria.

    The `PadlockArmsFilter` class checks whether a DNA sequence meets the criteria for forming stable padlock probe arms.
    This includes evaluating the minimum arm length and the melting temperature (Tm) difference between the arms.

    :param arm_length_min: The minimum required length for each padlock arm.
    :type arm_length_min: int
    :param arm_Tm_dif_max: The maximum allowed difference in melting temperature between the two arms.
    :type arm_Tm_dif_max: float
    :param arm_Tm_min: The minimum allowed melting temperature for the padlock arms.
    :type arm_Tm_min: float
    :param arm_Tm_max: The maximum allowed melting temperature for the padlock arms.
    :type arm_Tm_max: float
    :param Tm_parameters: Parameters for calculating the melting temperature.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
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
        """
        Calculate the melting temperatures and the ligation site for the padlock probe arms,
        determining if the sequence can form a valid padlock probe based on the specified parameters.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence meets the criteria for padlock probe arms, `False` otherwise.
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
    """
    A filter class for evaluating the suitability of detection oligonucleotides based on specific length, thymine content, and thermodynamic criteria.

    The `DetectionOligoFilter` class checks whether a DNA sequence meets the criteria for forming valid detection oligonucleotides.
    It evaluates the minimum and maximum lengths, thymine content, and melting temperature (Tm) differences between the padlock probe arms.

    :param detect_oligo_length_min: The minimum required length for the detection oligo.
    :type detect_oligo_length_min: int
    :param detect_oligo_length_max: The maximum allowed length for the detection oligo.
    :type detect_oligo_length_max: int
    :param min_thymines: The minimum required number of thymines in the detection oligo.
    :type min_thymines: int
    :param arm_length_min: The minimum required length for each padlock arm.
    :type arm_length_min: int
    :param arm_Tm_dif_max: The maximum allowed difference in melting temperature between the two arms.
    :type arm_Tm_dif_max: float
    :param arm_Tm_min: The minimum allowed melting temperature for the padlock arms.
    :type arm_Tm_min: float
    :param arm_Tm_max: The maximum allowed melting temperature for the padlock arms.
    :type arm_Tm_max: float
    :param Tm_parameters: Parameters for calculating the melting temperature.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
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
        """
        Evaluate if the sequence can form stable padlock arms and
        siutable detection oligos based on length and thymine content.

        :param sequence: The nucleotide sequence.
        :type sequence: Seq
        :return: `True` if the sequence meets the criteria for detection oligonucleotides, `False` otherwise.
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
