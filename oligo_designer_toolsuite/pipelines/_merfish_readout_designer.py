import random

import numpy as np
from Bio import SeqIO

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
            use_default_readouts=False
    ):

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
        self.oligo_20nt_dict = {k: v[:-5] for k, v in self.oligo_25nt_dict.items()}  # step 1

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

    def create_readouts(self, num_seq, seq_length=5):
        # Randomly select num_seq key-value pairs from the dictionary
        selected_pairs = dict(random.sample(self.oligo_25nt_dict.items(), num_seq))

        # step 1
        # Generate 5 random sequences of length seq_length and assign them randomly to the selected key-value pairs
        for key, value in selected_pairs.items():
            sequence = "".join(random.choices("ATCG", k=seq_length))
            if random.choice([True, False]):
                # Add the sequence to the end of the value
                selected_pairs[key] = value + sequence
            else:
                # Add the sequence to the start of the value
                selected_pairs[key] = sequence + value

        # step2
        # remove probs with significant homology to members of the transcriptome

        # step3
        # blast each potential readout probe against the previous build readout probs library
        # remove probs contains s contiguous stretch of homology > 14nt

        # step4
        # select subset of possible readout probs and built Blast library
        # remove probs contain a region of homology to another longer than 10 nt


        # Return the modified dictionary
        return selected_pairs

    def get_default_readouts(self):
        '''
        This function returns 16 validated readout probs from the merfish paper
        '''

        return self.default_readouts
