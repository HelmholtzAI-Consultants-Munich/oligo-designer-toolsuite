############################################
# imports
############################################

from Bio.SeqUtils import Seq
from Bio.SeqUtils import MeltingTemp as mt

from . import PropertyFilterBase

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

    ###TODO: move this function to utils as it is also used in the Tm filter once database refactor is merged
    def _get_Tm(self, sequence: Seq):
        """Internal method to calculate the melting temperature of a sequence.

        :param sequence: The DNA sequence for which Tm is calculated.
        :type sequence: Seq
        :return: The calculated melting temperature.
        :rtype: float
        """
        Tm = mt.Tm_NN(sequence, **self.Tm_parameters)
        if self.Tm_salt_correction_parameters is not None:
            Tm += mt.salt_correction(**self.Tm_salt_correction_parameters, seq=sequence)
        if self.Tm_chem_correction_parameters is not None:
            Tm = mt.chem_correction(Tm, **self.Tm_chem_correction_parameters)
        return round(Tm, 4)

    def _find_arms(self, sequence: Seq):
        """Internal method to identify the optimal ligation site in a DNA sequence for padlock probes by ensuring the arms formed are
        of sufficient length and their melting temperatures (Tm) are within specified constraints.

        The method iteratively adjusts the ligation site to find arm lengths and Tm values that meet the criteria
        of minimum arm length and Tm differences within the specified maximum. The process stops once suitable
        arms are found or if no configuration meets the criteria.

        :param sequence: The DNA sequence to analyze for padlock probe arm suitability.
        :type sequence: Seq
        :return: A tuple containing a boolean indicating if suitable arms were found, the Tm of the first and second
                 arm, the difference in Tm between arms, and the final ligation site position.
        :rtype: (bool, float, float, float, int)
        """

        len_sequence = len(sequence)
        ligation_site = len_sequence // 2

        arms_long_enough = (ligation_site >= self.arm_length_min) and (
            (len_sequence - ligation_site) >= self.arm_length_min
        )
        Tm_found = False
        sign_factor = 1  # switch between positive and negative shift
        shift = 1  # distance of ligation site shift

        while arms_long_enough and not Tm_found:
            Tm_arm1 = self._get_Tm(sequence[:ligation_site])
            Tm_arm2 = self._get_Tm(sequence[ligation_site:])
            Tm_dif = round(abs(Tm_arm2 - Tm_arm1), 2)
            Tm_found = (
                (Tm_dif <= self.arm_Tm_dif_max)
                and (self.arm_Tm_min <= Tm_arm1 <= self.arm_Tm_max)
                and (self.arm_Tm_min <= Tm_arm2 <= self.arm_Tm_max)
            )
            if not Tm_found:
                ligation_site += sign_factor * shift
                sign_factor *= -1
                shift += 1
                arms_long_enough = (ligation_site >= self.arm_length_min) and (
                    (len_sequence - ligation_site) >= self.arm_length_min
                )

        return Tm_found, Tm_arm1, Tm_arm2, Tm_dif, ligation_site

    def apply(self, sequence: Seq):
        """
        Applies the padlock arms filter to a DNA sequence.
        Applies the filter to evaluate if a given sequence is suitable for padlock probes based on
        arm length and melting temperature criteria.

        :param sequence: The DNA sequence to be checked.
        :type sequence: Seq
        :return: A tuple indicating if the sequence passes the filter and a dictionary containing
                 detailed results (arm Tm values, Tm difference, and ligation site) if the condition is met.
        :rtype: (bool, dict)
        """
        Tm_found, Tm_arm1, Tm_arm2, Tm_dif, ligation_site = self._find_arms(sequence)

        if Tm_found:
            return True, {
                "arm1_Tm": Tm_arm1,
                "arm2_Tm": Tm_arm2,
                "arms_Tm_dif": Tm_dif,
                "ligation_site": ligation_site,
            }
        else:
            return False, {}
