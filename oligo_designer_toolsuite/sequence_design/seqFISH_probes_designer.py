import numpy as np 
from Bio.Seq import Seq

class SeqfishProbesCreator:

     def __init__(
        self,

    ):
        self.n = 1
        


     def create_probes(self, oligos_DB, readout_probes, barcodes):
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
