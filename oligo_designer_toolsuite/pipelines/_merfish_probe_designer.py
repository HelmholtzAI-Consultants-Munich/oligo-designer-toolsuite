from itertools import combinations

import numpy as np
import pandas as pd
from scipy.spatial.distance import hamming

from oligo_designer_toolsuite.database import OligoDatabase

def make_barcode(raw_barcode: list, n_bits: int):
    """Creates the actual barcode from a list containing the indices of the "one" bits.

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
    """This function generates the codebook containing a collection of valid barcodes that fulfill the following constraints:
    - The barcodes have all the same `hamming_weight`, i.e. the number of bits containing 1.
    - The Hamming distance between every pair of valid barcodes is higher than the given threshold `min_hamming_dist`.
    The default parameters follow the MHD4 standard.
    
    :param n_regions: Number of regions that need to be assigned to a barcode.
    :type n_regions: int
    :param n_bits: Number of bits contained in each barcode, defaults to 16
    :type n_bits: int, optional
    :param min_hamming_dist: Minimum distance between two valid barcodes, defaults to 4
    :type min_hamming_dist: int, optional
    :param hamming_weight: Number of bits containing one in each barcode, defaults to 4
    :type hamming_weight: int, optional
    :return: Collection of all the valid barcodes generated
    :rtype: list
    """

    codebook = []
    for raw_barcode in combinations(iterable=range(n_bits), r=hamming_weight):
        new_barcode = make_barcode(raw_barcode=raw_barcode)
        # check if the barcode passes the requirements
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


def create_readout_probes_table(readout_probes: OligoDatabase, channels_ids: list, n_bits: int):
    assert len(readout_probes.database) >= n_bits, f"There are less readout probes ({len(readout_probes.database)}) than bits ({n_bits})."
    table = pd.DataFrame(columns=["bit", "channel", "readout_probe_id", "readout_probe_sequence"], index=list(range(n_bits)))
    n_channels = len(channels_ids)
    channel = 0
    for i, (readout_probe_id, readout_probe_features) in enumerate(readout_probes.database.items()):
        table.iloc[i] = [i, channels_ids[channel], readout_probe_id, readout_probe_features["oligo"]]
        channel = (channel + 1) % n_channels
        if i >= n_bits -1:
            break

    return table