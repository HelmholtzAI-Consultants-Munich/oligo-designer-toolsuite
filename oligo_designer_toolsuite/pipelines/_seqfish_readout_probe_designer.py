import os
import sys
sys.path.append(os.path.dirname(os.getcwd()))
from oligo_designer_toolsuite.sequence_design._readout_probes_generator import ReadoutProbes
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class SeqFISHReadoutProbeDesigner:
    """
    This class is designed to generate readout probes for SeqFISH+ experiment
    :param config: config file (dictionary) with all configurations for the experiment
    :type config: dict
    :param blastn: Blast filter, that was used during specificity filtering
    :type blastn: Blastn (SpecificityFilterBase)
    :param reference: reference DB, that was used during specificity filtering
    :type reference: ReferenceDatabase
    :param dir_output: output directory (used to save readout prbes as a separate file)
    :type dir_output: str
    """
    def __init__(self, config, blast, ref, dir_output) -> None:
        self.config = config
        self.blastn = blast
        self.reference = ref
        self.dir_output = dir_output

    def create_readout_probes(self):
        readout_generator = ReadoutProbes(length=self.config["length_readout"],  number_probes = self.config["num_pseudocolors"], GC_min = self.config["GC_min_readout"], GC_max= self.config["GC_max_readout"] ,number_consecutive = self.config["number_consecutive_readout"], random_seed = 0, blast_filter = self.blastn, reference_DB = self.reference)
        readout_probes = readout_generator.create_probes()
        output_file = os.path.join(self.dir_output, "readout_probes.fna")
        output_into_file = list()
        for i in range(0, len(readout_probes)):
            record = SeqRecord(readout_probes[i], id="readout_"+str(i+1), name="",description="")
            output_into_file.append(record)
        with open(output_file, "w") as f:
            SeqIO.write(output_into_file, f, "fasta")
        return readout_probes
        

        
               
