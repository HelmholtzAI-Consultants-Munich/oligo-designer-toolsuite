from Bio.SeqUtils import MeltingTemp as mt

from . import PropertyFilterBase


class PadlockArms(PropertyFilterBase):
    """Filters the sequences by arms melting temperature.

    :param min_arm_length: minimum arm length
    :type min_arm_length: int
    :param max_arm_Tm_dif: maximum difference between meltin temperature of arms
    :type max_arm_Tm_dif: float
    :param arm_Tm_min: minimum melting temperature
    :type arm_Tm_min: float
    :param arm_Tm_max: maximum melting temperature
    :type arm_Tm_max: float
    :param Tm_parameters: parameters to compute the melting temperature, for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_correction_parameters: parameters to correct the melting temperature,for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_correction_parameters: dict
    """

    def __init__(
        self,
        min_arm_length,
        max_arm_Tm_dif,
        arm_Tm_min,
        arm_Tm_max,
        Tm_parameters,
        Tm_correction_parameters,
    ) -> None:
        """Initialize the class"""

        super().__init__()
        self.min_arm_length = min_arm_length
        self.max_arm_Tm_dif = max_arm_Tm_dif
        self.arm_Tm_min = arm_Tm_min
        self.arm_Tm_max = arm_Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_correction_parameters = Tm_correction_parameters

    def __get_Tm(self, sequence):
        """Computes the melting temperature of the sequence.

        :param sequence: sequence for which the melting temperature is computed
        :type sequence: str
        :return: melting temperature
        :rtype: float
        """

        Tm = mt.Tm_NN(sequence, **self.Tm_parameters)
        Tm_corrected = round(mt.chem_correction(Tm, **self.Tm_correction_parameters), 2)
        return Tm_corrected

    def __find_arms(self, sequence):
        """Find the ligation site for the two padlock probe arms according melting temperature constraints
            The search starts in the center of the probe and shifts alternating 1 more nucleotide to the right and to
            the left.

        :param sequence: sequence for which the arms are computed
        :type sequence: str
        """

        probe_length = len(sequence)
        ligation_site = probe_length // 2

        arms_long_enough = (ligation_site >= self.min_arm_length) and (
            (probe_length - ligation_site) >= self.min_arm_length
        )
        Tm_found = False
        sign_factor = 1  # switch between positive and negative shift
        shift = 1  # distance of ligation site shift

        while arms_long_enough and not Tm_found:
            Tm_arm1 = self.__get_Tm(sequence[:ligation_site])
            Tm_arm2 = self.__get_Tm(sequence[ligation_site:])
            Tm_dif = round(abs(Tm_arm2 - Tm_arm1), 2)
            Tm_found = (
                (Tm_dif <= self.max_arm_Tm_dif)
                and (self.arm_Tm_min <= Tm_arm1 <= self.arm_Tm_max)
                and (self.arm_Tm_min <= Tm_arm2 <= self.arm_Tm_max)
            )
            if not Tm_found:
                ligation_site += sign_factor * shift
                sign_factor *= -1
                shift += 1
                arms_long_enough = (ligation_site >= self.min_arm_length) and (
                    (probe_length - ligation_site) >= self.min_arm_length
                )

        if Tm_found:
            arm_features = {
                "melt_temp_arm1": Tm_arm1,
                "melt_temp_arm2": Tm_arm2,
                "dif_melt_temp_arms": Tm_dif,
                "ligation_site": ligation_site,
            }
        else:
            arm_features = {}
        return Tm_found, arm_features

    def apply(self, sequence):
        """Applies the filter to the sequence.

        :param sequence: sequence for which the filter is applied
        :type sequence: str
        :return: Tm_found, arm_features
        :rtype: bool, dict
        """

        Tm_found, arm_features = self.__find_arms(sequence)
        return Tm_found, arm_features