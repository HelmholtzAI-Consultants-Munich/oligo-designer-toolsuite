from itertools import combinations

import numpy as np
from scipy.spatial.distance import hamming

def make_barcode(raw_barcode: list, n_bits: int):
    """Creates the actual barcode form a list containig hte inidices of the "one" bits.

    :param raw_barcode: list of the "one" bits in the barcode.
    :type raw_barcode: list
    :param n_bits: Number of bits contained in each barcode.
    :type n_bits: int, optional
    """

    barcode = np.zeros(n_bits, dtype=np.int8)
    for i in raw_barcode:
        barcode[i] = 1
    return barcode

def generate_codebook(n_regions: int, n_bits: int = 16, min_hamming_dist: int = 4, hamming_weight: int = 4):
    """This funtion generates the codebook containig a collection of valid barcodes that fulfill the following constraints:
    - The barcodes have all the same `hamming_weight`, i.e. the number of bits containing 1.
    - The Hamming distance between every pair of valid barcodes is higher than the given threshold `min_hamming_dist`.
    The defaults parameters follow the MHD4 standard.
    
    :param n_regions: Number of regions that need to be assigned to a barcode.
    :type n_regions: int
    :param n_bits: Number of bits contained in each barcode, defaults to 16
    :type n_bits: int, optional
    :param min_hamming_dist: Minimum distanced between two valid barcodes, defaults to 4
    :type min_hamming_dist: int, optional
    :param hamming_weight: Number of bits containing one in each barcode, defaults to 4
    :type hamming_weight: int, optional
    :return: Collection of all the valid barcodes ganarated
    :rtype: list
    """

    codebook = []
    for raw_barcode in combinations(iterable=range(n_bits), r=hamming_weight):
        new_barcode = make_barcode(raw_barcode=raw_barcode)
        # check if the barcode passes the requerements
        add_new_barcode = True
        for barcode in codebook:
            hamming_dist = hamming(new_barcode, barcode) * n_bits
            if hamming_dist < min_hamming_dist:
                add_new_barcode = False
                break
        if add_new_barcode:
            codebook.append(new_barcode)
    if len(codebook) < n_regions:
        raise ValueError(f"The number of valid barcodes ({len(codebook)}) is lower than the number of regions({n_regions}). Consider increasing the number of bits.")
    return codebook