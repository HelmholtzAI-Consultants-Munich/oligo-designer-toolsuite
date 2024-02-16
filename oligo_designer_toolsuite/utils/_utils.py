############################################
# imports
############################################

import csv

from Bio.SeqUtils import Seq
from Bio.SeqUtils import MeltingTemp as mt

############################################
# Collection of utility functions
############################################


def check_if_dna_sequence(seq: str, valid_characters={"A", "C", "T", "G"}) -> bool:
    return all(char.upper() in valid_characters for char in seq)


def check_if_list(obj):
    """
    Check if the input is a list, otherwise ensure that the input is converted to a list.

    :param obj: The input object.
    :type obj: Any
    :return: A list containing the input object or the input object itself if it's already a list.
    :rtype: list
    """
    if obj:
        obj = [obj] if not isinstance(obj, list) else obj
    return obj


def check_tsv_format(file):
    """Check if a given TSV file has any content.

    :param file: The path to the TSV file.
    :type file: str
    :return: True if the TSV file has content, False otherwise.
    :rtype: bool
    """
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)


def get_TmNN(
    sequence: Seq,
    Tm_parameters: dict,
    Tm_salt_correction_parameters: dict = None,
    Tm_chem_correction_parameters: dict = None,
):
    """Internal method to calculate the melting temperature of a sequence.

    :param sequence: The DNA sequence for which Tm is calculated.
    :type sequence: Seq
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
    :return: The calculated melting temperature.
    :rtype: float
    """
    Tm = mt.Tm_NN(sequence, **Tm_parameters)
    if Tm_salt_correction_parameters is not None:
        Tm += mt.salt_correction(**Tm_salt_correction_parameters, seq=sequence)
    if Tm_chem_correction_parameters is not None:
        Tm = mt.chem_correction(Tm, **Tm_chem_correction_parameters)
    return round(Tm, 4)
