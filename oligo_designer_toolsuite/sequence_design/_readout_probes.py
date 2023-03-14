from abc import ABC, abstractmethod
import numpy as np 
from oligo_designer_toolsuite.oligo_property_filter import GCContent, ConsecutiveRepeats
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

class ReadoutProbesBase(ABC):
    """"
    Abstract class to create readout probes 
    
    """

    @abstractmethod
    def create_readout_probes(self):
        """
        Function, that creates readout probes
        :return: list of readout sequences, that fulfil experiment contraints
        :rtype: list of Seq() 
        """

    def write_readout_probes(self, readout_probes, filename="readout_probes.fna"):
        output_file = os.path.join(self.dir_output, filename)
        output_into_file = list()
        for i in range(0, len(readout_probes)):
            record = SeqRecord(readout_probes[i], id="readout_"+str(i+1), name="",description="")
            output_into_file.append(record)
        with open(output_file, "w") as f:
            SeqIO.write(output_into_file, f, "fasta")

# TODO: Use database isntead of dict, documentation
class SeqFishReadoutProbes(ReadoutProbesBase):
    """"
    Class to create readout probes in SeqFISH+ experiment
    :param length: length of readout probe
    :type length: int
    :param number_probes: number of probes, that should be created (correspond to the num of pseudocolors in the experiment)
    :type number_probes: int
    :param GC_min: minimum GC-content of the probe
    :type GC_min: int
    :param GC_max: maximum GC-content of the probe
    :type GC_max: int
    :param number_consecutive: min num of consecutive nucleotides, that are not allowed in the probe
    :type number_consecutive: int
    :param blast_filter: Blast filter (typically this filter was already created during specificity filtering)
    :type blast_filter: BlastnFilter
    :param reference_DB: reference database, that was created for Blast specificity filtering
    :type reference_DB: ReferenceDatabase
    """
    def __init__(
        self,
        length: int,
        number_probes:int, 
        GC_min: int, 
        GC_max: int,
        number_consecutive: int,
        random_seed:int,
        dir_output: str,
        blast_filter ,
        reference_DB
    ):
    
        self.length = length
        self.num_probes = number_probes
        self.GC_content_filter = GCContent(GC_content_min = GC_min, GC_content_max = GC_max)
        self.proh_seq_filter = ConsecutiveRepeats(num_consecutive=number_consecutive)
        self.lib = ['A', 'C', 'G', 'T']
        self.blast_filter = blast_filter
        self.ref = os.path.basename(reference_DB.file_fasta) # good idea?
        self.seed = random_seed
        self.dir_output = dir_output
    
    
    def create_readout_probes(self):
        """
        Function, that creates readout probes
        :return: list of readout sequences, that fulfil experiment contraints
        :rtype: list of Seq() 
        """
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
            res = self.blast_filter._run_blast(sub_dict, 'gene', self.ref)
            for i in res.keys():
                probes.append(res[i]['sequence'])

        return probes[:self.num_probes]

    