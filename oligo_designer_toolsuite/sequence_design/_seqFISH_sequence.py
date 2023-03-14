import numpy as np 
from Bio.Seq import Seq

class SeqfishProbesCreator:
    """
    This class is used to assemble probes using primary probes, readout probes and barcodes, that were designed for each gene.
    """
    def __init__(
        self,
    ):
        pass
        
    def create_probes(self, oligos_DB, readout_probes, barcodes):
        """"
        Function for assembling probes.
        :param oligos_DB: database of oligos; dictionary, were each key corresponds to one gene
        :type oligos_DB: OligoDatabase
        :param readout_probes: list of readout probes
        :type readout_probes: list of Seq
        :param barcodes: barcodes for each gene from the database (each barcode consists of 4 digits)
        :type barcodes: dict (str : list of 4 int)
        :return: oligo_database, where "sequence" are assembled sequences
        :rtype: OligoDatabase
        """
        for i in oligos_DB.keys():
            barcode = barcodes[i]
            s1 = readout_probes[barcode[0]]
            s2 = readout_probes[barcode[1]]
            s3 = readout_probes[barcode[2]]
            s4 = readout_probes[barcode[3]]
            left = str(s1) + str(s2)
            right = str(s3) + str(s4)
            for j in oligos_DB[i].keys():
                seq = str(oligos_DB[i][j]['sequence'])
                seq = left + seq
                seq = seq + right
                oligos_DB[i][j]['sequence'] = Seq(seq)
        return oligos_DB
