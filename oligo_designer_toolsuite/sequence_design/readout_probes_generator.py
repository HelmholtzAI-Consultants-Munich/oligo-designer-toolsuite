import numpy as np 
from oligo_designer_toolsuite.oligo_property_filter import GCContent, ProhibitedSequences
from Bio.Seq import Seq
import os

class ReadoutProbes:
    
    def __init__(
        self,
        length: int,
        number_probes:int, 
        GC_min: int, 
        GC_max: int,
        number_consecutive: int,
        random_seed:int,
        blast_filter ,
        reference_DB
    ):
    
        self.length = length
        self.num_probes = number_probes
        self.GC_content_filter = GCContent(GC_content_min = GC_min, GC_content_max = GC_max)
        self.proh_seq_filter = ProhibitedSequences(num_consecutive=number_consecutive)
        self.lib = ['A', 'C', 'G', 'T']
        self.blast_filter = blast_filter
        #self.ref = os.path.basename(reference_DB.file_reference_DB)
        self.seed = random_seed
    
    
    def create_probes(self):
        probes = list()
        while len(probes) < self.num_probes:
            sub_dict = dict() # dict "probe" : "seq"
            ind = 0
            while ind < self.num_probes:
                seq_letters = np.random.choice(self.lib, self.length, replace = True)
                seq = ''.join(seq_letters)
                res_proh, _ = self.proh_seq_filter.apply(seq)
                res_GC, _ = self.GC_content_filter.apply(seq)
                if res_proh and res_GC:
                    sub_dict ['prob' + str(ind) ]= {"sequence" : Seq(seq)}
                    ind += 1
            #res = self.blast_filter._run_blast(sub_dict, 'gene', self.ref)
            for i in res.keys():
                probes.append(res[i]['sequence'])

        return probes[:self.num_probes]

    