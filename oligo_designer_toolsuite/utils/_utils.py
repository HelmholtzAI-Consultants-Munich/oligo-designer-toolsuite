############################################
# imports
############################################

import csv

from Bio.SeqUtils import Seq
from Bio.SeqUtils import MeltingTemp as mt

############################################
# Collection of utility functions
############################################


def check_if_dna_sequence(seq: str, valid_characters=["A", "C", "T", "G"]):
    """Verifies if a given sequence consists solely of valid DNA nucleotide characters (A, C, T, G).
    This function is case-insensitive and validates the entire sequence against a set of allowed characters.

    :param seq: The DNA sequence to validate.
    :type seq: str
    :param valid_characters: A list of characters representing valid nucleotides. Defaults to {"A", "C", "T", "G"}.
    :type valid_characters: list
    :return: True if the sequence is valid DNA, False otherwise.
    :rtype: bool
    """
    return all(char.upper() in valid_characters for char in seq)


def check_if_key_exists(nested_dict: dict, key: str):
    """Recursively check if a given key exists in a nested dictionary. The function searches through all levels of the nested dictionary to find the key.

    :param nested_dict: The nested dictionary to search for the key.
    :type nested_dict: dict
    :param key: The key to search for within the nested dictionary.
    :type key: str
    :return: True if the key is found anywhere within the nested dictionary, False otherwise.
    :rtype: bool
    """
    if isinstance(nested_dict, dict):
        if key in nested_dict:
            return True
        else:
            for value in nested_dict.values():
                if check_if_key_exists(value, key):
                    return True
    return False


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


def check_tsv_format(file: str):
    """Check if a given TSV file has any content.

    :param file: The path to the TSV file.
    :type file: str
    :return: True if the TSV file has content, False otherwise.
    :rtype: bool
    """
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)


# TODO: move to oligo database attributes class
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
