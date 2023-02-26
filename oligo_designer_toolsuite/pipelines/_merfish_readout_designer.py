import random

import numpy as np
import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_filter import GCContent, ConsecutiveRepeats
from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilter,
    Blastn,
)
from oligo_designer_toolsuite.database import ReferenceDatabase

from Bio.Seq import Seq
import os
from pathlib import Path



class ReadoutProbes:

    def __init__(
            self,
            config,
            dir_output,
            file_transcriptome,
            region_generator,
            primer_fasta_file=None

    ):
        self.config =config
        self.dir_output = os.path.join(dir_output, "readout_probes")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)
        self.file_transcriptome=file_transcriptome
        self.region_generator=region_generator
        self.primer_fasta_file=primer_fasta_file

        length =self.config["readout_oligo"]["oligo_length_max"]
        self.readout_oligo_config = self.config["readout_oligo"]

       
        self.oligo_25nt_path = self.config["readout_oligo"]["file_bc25mer"]
        self.oligo_25nt_dict = {rec.id: rec.seq for rec in SeqIO.parse(self.oligo_25nt_path, "fasta")}
        self.oligo_30nt_dict = {}

        #select num_seq random sequences 
        num_seq=1000
        sequences = list(self.oligo_25nt_dict.keys())
        selected_sequences =random.sample(sequences, num_seq)
        for key in selected_sequences:
            value=self.oligo_25nt_dict[key]
            sequence = "".join(random.choices("ATCG", k=length - len(value)))
            if random.choice([True, False]):
                # Add the sequence to the end of the value
                self.oligo_30nt_dict[key] = value + sequence
            else:
                # Add the sequence to the start of the value
                self.oligo_30nt_dict[key] = sequence + value
        self.oligo_30nt_path = os.path.join(self.readout_oligo_config["oligo_output"], "bc30mer_oligo.fna")
        output = []
        for name, seq in self.oligo_30nt_dict.items():
            output.append(SeqRecord(seq, name, "", ""))
        with open(self.oligo_30nt_path, "w") as handle:
            SeqIO.write(output, handle, "fasta")


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

        self.oligo_database = OligoDatabase(
            file_fasta=self.oligo_30nt_path,
            oligo_length_min=self.readout_oligo_config["oligo_length_min"],
            oligo_length_max=self.readout_oligo_config["oligo_length_max"],
            n_jobs=1,
            dir_output=self.dir_output,
        )
        self.oligo_database.create_database()

    def create_readouts(self, num_readouts):

         # blast each potential readout probe against the previous build primer probs library
        dir_specificity = os.path.join(self.dir_output, "specificity_temporary1")
        blast_filter = Blastn(
            dir_specificity=dir_specificity,
            word_size=self.config["readout_blast_setup"]['blast1_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        reference_database = ReferenceDatabase(
            file_fasta=self.primer_fasta_file)
        reference_database.load_fasta_into_database()
        specificity_filter = SpecificityFilter(filters=[blast_filter],
                                                write_regions_with_insufficient_oligos=self.config[
                                                    "write_removed_genes"])
        self.oligo_database = specificity_filter.apply(oligo_database=self.oligo_database,
                                                   reference_database=reference_database, n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database = self.oligo_database.write_database(filename="readout_database_specificity_filter1.txt")
        
        # blast each potential readout probe against the previous built readout probes library
        dir_specificity1 = os.path.join(self.dir_output, "specificity_temporary1")
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
        self.oligo_database = specificity_filter1.apply(oligo_database=self.oligo_database,
                                                   reference_database=reference_database1, n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database =self.oligo_database.write_database(filename="readout_database_specificity_filter1.txt")

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
            dir_specificity=dir_specificity2, 
            word_size=self.config["readout_blast_setup"]["blast2_word_size"],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        # initialize the specificity filter class
        specificity_filter2 = SpecificityFilter(filters=[blast_filter2], write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        self.oligo_database = specificity_filter2.apply(oligo_database=self.oligo_database, reference_database=reference_database2, n_jobs=self.config["n_jobs"])

        file_database = self.oligo_database.write_database(filename="readout_database_full.txt")

        
        #return a specified number of readout probes (num_readouts)
        readout_oligos_dict = {}
        readout_genes = list(self.oligo_database.database.keys())[0:num_readouts]
        readout_oligo_ids = [list(self.oligo_database.database[gene].keys())[0] for gene in readout_genes]
        for gene, oligo_id in zip(readout_genes, readout_oligo_ids):
            readout_oligos_dict[gene] = str(self.oligo_database.database[gene][oligo_id]["sequence"])
        return list(readout_oligos_dict.values())

    def get_default_readouts(self):
        '''
        This function returns 16 validated readout probs from the merfish paper
        '''

        return self.default_readouts
