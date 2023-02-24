import random

import numpy as np
import yaml
from Bio import SeqIO

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_filter import GCContent, ConsecutiveRepeats

from Bio.Seq import Seq
import os


class ReadoutProbes:

    def __init__(
            self,
            length: int,
            number_probes: int,
            GC_min: int,
            GC_max: int,
            blast_filter,
            reference_DB,
            config_path,
            use_default_readouts=False
    ):

        with open(config_path, 'r') as yaml_file:
            self.config = yaml.safe_load(yaml_file)

        self.readout_oligo_config = self.config["readout_oligo"]
        self.length = length
        self.num_probes = number_probes
        self.GC_content_filter = GCContent(GC_content_min=GC_min, GC_content_max=GC_max)
        self.repeat_num_max = ConsecutiveRepeats(3)
        self.lib = ['A', 'C', 'G', 'T']
        self.blast_filter = blast_filter
        self.ref = os.path.basename(reference_DB.file_reference_DB)
        self.use_default_readouts = use_default_readouts
        self.oligo_25nt_path = os.path.join("data", "bc25mer.240k.fasta")
        self.oligo_25nt_dict = SeqIO.to_dict(SeqIO.parse(self.oligo_25nt_path, "fasta"))
        self.oligo_30nt_dict = {}
        self.oligo_30nt_DB = None

        oligo_30nt_dict = {}
        for key, value in self.oligo_25nt_dict.items():
            sequence = "".join(random.choices("ATCG", k=length - len(value)))
            if random.choice([True, False]):
                # Add the sequence to the end of the value
                self.oligo_30nt_dict[key] = value + sequence
            else:
                # Add the sequence to the start of the value
                self.oligo_30nt_dict[key] = sequence + value
        self.oligo_30nt_path = os.path.join(self.readout_oligo_config["oligo_output"], "bc30mer_oligo.fasta")

        with open(self.oligo_30nt_path, "w") as handle:
            for name, seq in self.oligo_30nt_dict.items():
                handle.write(">" + name + "\n")
                handle.write(str(seq) + "\n")
            print("write bc30mer_oligo done!")


        self.default_readouts = [
            "CGCAACGCTTGGGACGGTTCCAATCGGATC",
            "CGAATGCTCTGGCCTCGAACGAACGATAGC",
            "ACAAATCCGACCAGATCGGACGATCATGGG",
            "CAAGTATGCAGCGCGATTGACCGTCTCGTT",
            "GCGGGAAGCACGTGGATTAGGGCATCGACC",
            "AAGTCGTACGCCGATGCGCAGCAATTCACT",
            "CGAAACATCGGCCACGGTCCCGTTGAACTT",
            "ACGAATCCACCGTCCAGCGCGTCAAACAGA",
            "CGCGAAATCCCCGTAACGAGCGTCCCTTGC",
            "GCATGAGTTGCCTGGCGTTGCGACGACTAA",
            "CCGTCGTCTCCGGTCCACCGTTGCGCTTAC",
            "GGCCAATGGCCCAGGTCCGTCACGCAATTT",
            "TTGATCGAATCGGAGCGTAGCGGAATCTGC",
            "CGCGCGGATCCGCTTGTCGGGAACGGATAC",
            "GCCTCGATTACGACGGATGTAATTCGGCCG",
            "GCCCGTATTCCCGCTTGCGAGTAGGGCAAT",
        ]
        self.oligo_readout_DB = ReadoutProbes.create_readout_DB()
    def create_readout_DB(self):
        oligo_database = OligoDatabase(
            file_fasta=self.oligo_30nt_path,
            oligo_length_min=self.readout_oligo_config["oligo_length_min"],
            oligo_length_max=self.readout_oligo_config["oligo_length_max"],
            n_jobs=1,
            dir_output=self.readout_oligo_config["oligo_DB_output"],
        )
        return oligo_database.create_database()

    def create_readouts(self):

        # step 1
        # Generate 5 random sequences of length seq_length and assign them randomly to the selected key-value pairs
        # sequences in self.oligo_readout_DB

        # step2
        # remove probs with significant homology to members of the transcriptome

        # step3
        # blast each potential readout probe against the previous build readout probs library
        # remove probs contains s contiguous stretch of homology > 14nt

        # step4
        # select subset of possible readout probs and built Blast library
        # remove probs contain a region of homology to another longer than 10 nt


        # Return the modified dictionary
        pass

    def get_default_readouts(self):
        '''
        This function returns 16 validated readout probs from the merfish paper
        '''

        return self.default_readouts
