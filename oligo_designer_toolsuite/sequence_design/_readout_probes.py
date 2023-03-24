import os
import random
from abc import ABC, abstractmethod

import numpy as np
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_property_filter import (ConsecutiveRepeats,
                                                            GCContent,
                                                            PropertyFilter)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    Blastn, ExactMatches, SpecificityFilter)

from oligo_designer_toolsuite.utils._sequence_design import generate_random_sequence


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

    def _get_readout_probes_list(self, n_readout_probes):
        """Rerurns readout probes list or db?"""
        #return a specified number of readout probes (num_readouts)
        readout_oligos_dict = {}
        readout_genes = list(self.oligo_database.database.keys())[0:n_readout_probes]
        readout_oligo_ids = [list(self.oligo_database.database[gene].keys())[0] for gene in readout_genes]
        
        for gene, oligo_id in zip(readout_genes, readout_oligo_ids):
            readout_oligos_dict[gene] = str(self.oligo_database.database[gene][oligo_id]["sequence"])
        return list(readout_oligos_dict.values())

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
        reference_DB, 
        sequence_alphabet: list[str] = ['A', 'C', 'G', 'T']
    ):
    
        self.length = length
        self.num_probes = number_probes
        self.property_filters = [
            GCContent(GC_content_min = GC_min, GC_content_max = GC_max),
            ConsecutiveRepeats(num_consecutive=number_consecutive),
        ]
        self.specificity_filters = [blast_filter]

        self.sequence_alphabet = sequence_alphabet
        self.ref = os.path.basename(reference_DB.file_fasta)
        self.seed = random_seed
        self.dir_output = dir_output
        self.readout_database = None
    
    
    def create_readout_probes(self, property_filter, specificity_filter):
        """
        Function, that creates readout probes
        :return: list of readout sequences, that fulfil experiment contraints
        :rtype: list of Seq()  
        """

        self.readout_database = OligoDatabase(file_fasta=None, dir_output=self.dir_output)

        self.readout_database.create_random_database(
            self.length * 100,
            self.num_probes,
            sequence_alphabet = self.sequence_alphabet)
        
        
        # property_filter = PropertyFilter(self.property_filters)
        self.readout_database = property_filter.apply(self.readout_database)

        # specificity_filter = SpecificityFilter(self.specificity_filters)
        self.readout_database = specificity_filter.apply(self.readout_database, self.ref)


class MerfishReadoutProbes(ReadoutProbesBase):

    def __init__(
            self,
            readout_oligo_length,
            oligo_25nt_file,
            dir_output,
            file_transcriptome,
            region_generator,
            num_sequences = 1000,
            primer_fasta_file=None

    ):
        """
        This class is used to design the readout probes.

        :param dir_output: output directory
        :type dir_output: str
        :param file_transcriptome: directory of the fasta file for the transcriptome
        :type file_transcriptome: str
        :param region_generator: region generator used to create the file_transcriptome
        :type region_generator: CustomGenomicRegionGenerator
        :param primer_fasta_file: directory of the fasta file where previously generated primers are stored
        :type primer_fasta_file: str
        """
        self.dir_output = dir_output
        self.file_transcriptome=file_transcriptome
        self.region_generator=region_generator
        self.primer_fasta_file=primer_fasta_file

        oligo_25nt_sequences = [rec.seq for rec in SeqIO.parse(oligo_25nt_file, "fasta")]
        oligo_25nt_sequences = random.sample(oligo_25nt_sequences, num_sequences)

        augmented_sequences = []
        for sequence in oligo_25nt_sequences:
            random_sequence = generate_random_sequence(readout_oligo_length - len(sequence))
            append_beginning = random.choice([True, False])
            if append_beginning:
                augmented_sequences.append(random_sequence + sequence)
            else:
                augmented_sequences.append(sequence + random_sequence)


        oligo_30_nt_sequences = {'bc30mer_' + str(i) : sequence for i, sequence in enumerate(augmented_sequences)}
        self.readout_database = OligoDatabase(file_fasta=None, dir_output=self.dir_output)

        self.readout_database.create_database_from_sequences(oligo_30_nt_sequences)

    def create_readout_probes(self, specificity_filter, specificity_filter_1, specificity_filter_2):
        '''
        Function to create the readout probes
        param num_readouts: number of readout probes that should be created
        type num_readouts: int
        '''
        if (self.primer_fasta_file is not None):
            # blast each potential readout probe against the previous build primer probs library
            # dir_specificity = os.path.join(self.dir_output, "specificity_temporary1")
            # exact_matches = ExactMatches(dir_specificity=dir_specificity)
            # blast_filter = Blastn(
            #     dir_specificity=dir_specificity,
            #     word_size=self.config["readout_blast_setup"]['blast1_word_size'],
            #     percent_identity=self.config["percent_identity"],
            #     coverage=self.config["coverage"],
            #     strand=self.config["strand"],
            # )
            # filters = [exact_matches, blast_filter]
            reference_database = ReferenceDatabase(
                file_fasta=self.primer_fasta_file)
            # reference_database.load_fasta_into_database()
            # specificity_filter = SpecificityFilter(filters=filters,
            #                                         write_regions_with_insufficient_oligos=self.config[
            #                                             "write_removed_genes"])
            self.oligo_database = specificity_filter.apply(oligo_database=self.oligo_database,
                                                       reference_database=reference_database, n_jobs=self.config["n_jobs"])
            # if self.config["write_intermediate_steps"]:
            #     file_database = self.oligo_database.write_database(filename="readout_database_specificity_filter1.txt")
        
        # blast each potential readout probe against the previous built readout probes library
        # dir_specificity1 = os.path.join(self.dir_output, "specificity_temporary1")
        # blast_filter1 = Blastn(
        #     dir_specificity=dir_specificity1,
        #     word_size=self.config["readout_blast_setup"]['blast1_word_size'],
        #     percent_identity=self.config["percent_identity"],
        #     coverage=self.config["coverage"],
        #     strand=self.config["strand"],
        # )
        #  # create reference DB with fasta file
        reference_database1 = ReferenceDatabase(
            file_fasta=self.oligo_30nt_path)
        # reference_database1.load_fasta_into_database()
        # specificity_filter1 = SpecificityFilter(filters=[blast_filter1],
        #                                         write_regions_with_insufficient_oligos=self.config[
        #                                             "write_removed_genes"])
        self.oligo_database = specificity_filter_1.apply(oligo_database=self.oligo_database,
                                                   reference_database=reference_database1, n_jobs=self.config["n_jobs"])
        # if self.config["write_intermediate_steps"]:
        #     file_database =self.oligo_database.write_database(filename="readout_database_specificity_filter1.txt")

        # remove probs with significant homology to members of the transcriptome
        # dir_specificity2 = os.path.join(self.dir_output, "specificity_temporary2") # folder where the temporary files will be written

        reference_database2 = ReferenceDatabase(
            file_fasta = self.file_transcriptome,
            files_source = self.region_generator.files_source,
            species = self.region_generator.species,
            annotation_release = self.region_generator.annotation_release,
            genome_assembly = self.region_generator.genome_assembly,
            dir_output=self.dir_output
        )
        reference_database2.load_fasta_into_database()
        # blast_filter2 = Blastn(
        #     dir_specificity=dir_specificity2, 
        #     word_size=self.config["readout_blast_setup"]["blast2_word_size"],
        #     percent_identity=self.config["percent_identity"],
        #     coverage=self.config["coverage"],
        #     strand=self.config["strand"],
        # )
        # # initialize the specificity filter class
        # specificity_filter2 = SpecificityFilter(filters=[blast_filter2], write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        self.oligo_database = specificity_filter_2.apply(oligo_database=self.oligo_database, reference_database=reference_database2, n_jobs=self.config["n_jobs"])

        file_database = self.oligo_database.write_database(filename="readout_database_full.txt")

        
        

    def create_default_readout_probes(self, dir_output):
        '''
        This function returns 16 validated readout probs from the merfish paper
        '''
        self.readout_database = OligoDatabase(file_fasta=None, dir_output=self.dir_output)
        
        default_readout_probes = {
            'merfish_readout_probes' : [
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
        ]}

        self.readout_database.create_database_from_sequences(default_readout_probes)
