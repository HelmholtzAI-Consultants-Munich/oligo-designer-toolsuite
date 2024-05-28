from itertools import product

import numpy as np
import pandas as pd

from oligo_designer_toolsuite.database import OligoDatabase

class Barcode():
    """Class to generate channel-wise barcodes, i.e. barcodes where each pseudocolor belongs to the same channel.
    To allow for error corrections an additional barcode round is added and the barcodes are generated to be robust against a sigle deletion.
    For example for a barcode with (i,j,k) pseudocolors, to the additional barcode round we will assig the `(i+j+k) mod n_pseudocolors` pseoudocolor.

    :param pseudocolors: List of the pseoudocolors contained in the barcode.
    :type pseudocolors: list
    :param channel: channel to wich the barcode belongs to.
    :type channel: int
    :param n_pseudocolors: Total number of pseudocolors.
    :type n_pseudocolors: int
    :param n_channels: Total number fo channels
    :type n_channels: int
    """
    def __init__(self, pseudocolors: list, channel: int, n_pseudocolors: int, n_channels: int):
        """init method
        """
        self.pseudocolors = pseudocolors + [sum(pseudocolors)%n_pseudocolors]
        print("next")
        print(self.pseudocolors)
        self.n_pseudocolors = n_pseudocolors
        self.n_channels = n_channels
        self.channel = channel

    def to_binary(self):
        """The method tranforms the barcode in a binary barcode.

        :return: Binary barcode.
        :rtype: np.Array
        """
        assert self.n_pseudocolors > max(self.pseudocolors), f"The number of pseudocolor is {self.n_pseudocolors}, while the barcode contains {max(self.pseudocolors)} pseudocolors."
        assert self.n_channels > self.channel, f"The number of channles is {self.n_channels}, while the barcode contains {self.channel} channels."
        n_barcode_rounds = len(self.pseudocolors)
        barcode = np.zeros(self.n_channels*self.n_pseudocolors*n_barcode_rounds, dtype=np.int8)
        for i, pseudocolor in enumerate(self.pseudocolors):
            barcode[i*self.n_pseudocolors*self.n_channels + self.n_channels*pseudocolor + self.channel] = 1
        return barcode

def generate_codebook(n_regions: int, n_barcode_rounds: int = 4, n_pseudocolors: int = 20, n_channels: int = 3):
    """This function generates the codebook containig a collection of barcodes originated form the given barcode rounds, pseudocolors and channels.
    Specifically, each barcode is generated such that all the readout probes encoding for one region is associated to the same channel.
    To allow for error corrections the barcodes are generated to be robust against a sigle deletion. To do so, the pseudocolor assignoed to the last barcode round
    is deterministically derived from the other pseoudocolors.
    
    :param n_regions: Number of regions that need to be assigned to a barcode.
    :type n_regions: int
    :param n_barcode_rounds: Number of barcode runs, defaults to 4.
    :type n_barcode_rounds: int, optional
    :param n_pseudocolors: Number of pseudocolors, defaults to 20.
    :type n_pseudocolors: int, optional
    :param n_channels: Number of channels, defaults to 3.
    :type n_channels: int, optional
    :return: Collection of all the valid barcodes generated
    :rtype: list
    """
    codebook = []
    codebook_size = n_channels * (n_pseudocolors**(n_barcode_rounds - 1))
    if codebook_size < n_regions:
        raise ValueError(f"The number of valid barcodes ({codebook_size}) is lower than the number of regions ({n_regions}). Consider increasing the number of psudocolors or barcoding rounds.")
    for pseudocolors in product(range(n_pseudocolors), repeat=n_barcode_rounds - 1):
        pseudocolors = list(pseudocolors)
        for channel in range(n_channels):
            barcode = Barcode(pseudocolors=pseudocolors, channel=channel, n_pseudocolors=n_pseudocolors, n_channels=n_channels)
            codebook.append(barcode.to_binary())
    return codebook    


def create_readout_probes_table(readout_probes : OligoDatabase, channels_ids: list, n_barcode_rounds: int = 4, n_pseudocolors: int = 20):
    """This function associates the readout probes with the bit of the barcodes and assigns them to a specific barcoding round, pseudocolor and channel.

    :param readout_probes: Oligo databse containing all the readout probes.
    :type readout_probes: OligoDatabase
    :param channels_ids: Names of the available channles.
    :type channels_ids: list
    :param n_barcode_rounds: Number of barcode runs, defaults to 4
    :type n_barcode_rounds: int, optional
    :param n_pseudocolors: Number of pseudocolors, defaults to 20
    :type n_pseudocolors: int, optional
    :return: table containig the information for each readout probe.
    :rtype: pd.DataFrame
    """
    n_channels = len(channels_ids)
    n_bits = n_barcode_rounds*n_pseudocolors*n_channels
    assert len(readout_probes.database) >= n_bits, f"There are less readout probes ({len(readout_probes.database)}) than bits ({n_bits})."
    table = pd.DataFrame(columns=["bit", "barcode_round", "pseudocolor", "channel", "readout_probe_id", "readout_probe_sequence"], index=list(range(n_bits)))
    barcode_round = 0
    pseudocolor = 0
    channel = 0
    for i, (oligo_id, oligo_features) in enumerate(readout_probes.database.items()):
        table.iloc[i] = [i, barcode_round, pseudocolor, channels_ids[channel], oligo_id, oligo_features["oligo"]]
        channel = (channel + 1) % n_channels
        if channel == 0:
            pseudocolor = (pseudocolor + 1) % n_pseudocolors
            if pseudocolor == 0:
                barcode_round = (barcode_round + 1) % n_barcode_rounds
        if i >= n_bits-1:
            break
    return table
