import random

import numpy as np
import yaml
from Bio import SeqIO

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_filter import GCContent, ConsecutiveRepeats

from Bio.Seq import Seq
import os
from pathlib import Path


class ReadoutProbes:

    def __init__(
            self,
            config,
            dir_output,
            file_transcriptome,
            region_generator

    ):
        self.config =config
        self.dir_output = os.path.join(dir_output, "readout_probes")
        self.file_transcriptome=file_transcriptome
        self.region_generator=region_generator

        length =self.config["readout_oligo"]["oligo_length_max"]
        self.readout_oligo_config = self.config["readout_oligo"]
        self.oligo_25nt_path = self.config["readout_oligo"]["file_bc25mer"]#os.path.join("data", "bc25mer.240k.fasta")
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

        oligo_database = OligoDatabase(
            file_fasta=self.oligo_30nt_path,
            oligo_length_min=self.readout_oligo_config["oligo_length_min"],
            oligo_length_max=self.readout_oligo_config["oligo_length_max"],
            n_jobs=1,
            dir_output=self.readout_oligo_config["oligo_DB_output"],
        )
        self.oligo_readout_DB =  oligo_database.create_database()

    def create_readouts(self):
        self.oligo_readout_DB.write_database(filename="merfish_readout_probes_init.txt")
        #select a subset of the database

    




 
        # blast each potential readout probe against the previous build readout probs library
        dir_specificity1 = os.path.join(self.config["dir_output"], "specificity_temporary1")
        blast_filter1 = Blastn(
            dir_specificity=dir_specificity1,
            word_size=self.config["readout_blast_setup"]['blast1_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
         # create reference DB with fasta file
        reference_database1 = ReferenceDatabase(
            file_fasta=self.oligo_30nt_path)
        reference_database1.load_fasta_into_database()
        specificity_filter1 = SpecificityFilter(filters=[blast_filter1],
                                                write_regions_with_insufficient_oligos=self.config[
                                                    "write_removed_genes"])
        oligo_database = specificity_filter1.apply(oligo_database=oligo_database,
                                                   reference_database=reference_database1, n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="readout_database_specificity_filter1.txt")

        # remove probs with significant homology to members of the transcriptome
        dir_specificity2 = os.path.join(self.dir_output, "specificity_temporary2") # folder where the temporary files will be written

        reference_database2 = ReferenceDatabase(
            file_fasta = self.file_transcriptome,
            files_source = self.region_generator.files_source,
            species = self.region_generator.species,
            annotation_release = self.region_generator.annotation_release,
            genome_assembly = self.region_generator.genome_assembly,
            dir_output=self.dir_output
        )
        reference_database2.load_fasta_into_database()
        blast_filter2 = Blastn(
            dir_specificity=dir_specificity, 
            word_size=self.config["readout_blast_setup"]["blast2_word_size"],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        # initialize the specificity filter class
        specificity_filter2 = SpecificityFilter(filters=[blast_filter2], write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        oligo_database = specificity_filter2.apply(oligo_database=oligo_database, reference_database=reference_database2, n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="readout_database_specificity_filter2.txt")




        # Return the readout probes
        pass

    def get_default_readouts(self):
        '''
        This function returns 16 validated readout probs from the merfish paper
        '''

        return self.default_readouts
