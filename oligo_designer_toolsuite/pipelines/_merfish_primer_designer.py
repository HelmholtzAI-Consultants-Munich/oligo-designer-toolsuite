import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import yaml
from oligo_designer_toolsuite.oligo_property_filter import (
    PropertyFilter,
    MaskedSequences,
    GCContent,
    MeltingTemperatureNN,
    PadlockArms,
    ConsecutiveRepeats,
    GCClamp
)
from oligo_designer_toolsuite.database import ReferenceDatabase, OligoDatabase

from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilter,
    Blastn,
)
import random
from pathlib import Path


class PrimerProbes:

    def __init__(
            self,
            num_seq,
            config,
            dir_output,
            file_transcriptome,
            region_generator

    ):
        """
        This class is used to design the primer probes.

        :param num_seq: number of primer probes which should be created
        :type num_seq: int
        :param config: config file
        :type config: file pointer
        :param dir_output: output directory
        :type dir_output: str
        :param file_transcriptome: directory of the fasta file for the transcriptome
        :type file_transcriptome: str
        :param region_generator: region generator used to create the file_transcriptome
        :type region_generator: CustomGenomicRegionGenerator
        """
        self.num_seq = num_seq
        self.length = 25  # need to trim the sequece
       
        self.config = config
        self.primer_oligo_config = self.config["primer_oligo"]

        self.dir_output = os.path.join(dir_output, "primer_probes")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.file_transcriptome = file_transcriptome
        self.region_generator = region_generator
        self.config_param = self.config["primers_setup"].copy()
        self.GC_min = self.config_param["GC_content_min"]
        self.GC_max = self.config_param["GC_content_max"]
        self.GC_content_filter = GCContent(GC_content_min=self.GC_min, GC_content_max=self.GC_max)
        self.repeat_num_max = ConsecutiveRepeats(3)
        self.lib = ['A', 'C', 'G', 'T']
        self.T7promoter = "TAATACGACTCACTATAG"
    
        self.oligo_25nt_path = self.primer_oligo_config["file_bc25mer"]

        self.oligo_25nt_dict = {rec.id: rec.seq for rec in SeqIO.parse(self.oligo_25nt_path, "fasta")}
        #self.oligo_20nt_dict = {k: v[:-5] for k, v in self.oligo_25nt_dict.items()}
        #select n random sequences 
        n=num_seq*50 #start with 50 times the number of required primers
        if (n>240000):
            n=240000
        sequences = list(self.oligo_25nt_dict.keys())
        selected_sequences =random.sample(sequences, n)
        self.oligo_20nt_dict={}
        for key in selected_sequences:
            value=self.oligo_25nt_dict[key]
            self.oligo_20nt_dict [key]= value[:-5] 

        self.oligo_20nt_path = os.path.join(self.primer_oligo_config["oligo_output"], "primer_oligo.fasta")
        with open(self.oligo_20nt_path, "w") as handle:
            for name, seq in self.oligo_20nt_dict.items():
                handle.write(">" + name + "\n")
                handle.write(str(seq) + "\n")
            print("write bc20mer_oligo done!")

        self.consecutive_repeats = ConsecutiveRepeats(self.config_param["max_repeats_nt"])
        self.GC_clamp = GCClamp(self.config_param["GC_clamp_n"])

        # create the list of filters
        # self.melting_temperature

        self.filters = [self.GC_content_filter, self.consecutive_repeats, self.GC_clamp]

        # specificty filters

        # first blast against human transcriptome
        self.dir_specificity1 = os.path.join(self.dir_output, "specificity_temporary1")
        self.blast_filter1 = Blastn(
            dir_specificity=self.dir_specificity1,
            word_size=self.config["primers_blast_setup"]['blast1_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        self.reference_database1 = ReferenceDatabase(
            file_fasta=self.file_transcriptome,
            files_source=self.region_generator.files_source,
            species=self.region_generator.species,
            annotation_release=self.region_generator.annotation_release,
            genome_assembly=self.region_generator.genome_assembly,
            dir_output=self.dir_output
        )
        # self.reference_database1.load_fasta_into_database()

        # second blast against 3' end of other primers
        self.dir_specificity2 = os.path.join(self.dir_output, "specificity_temporary2")
        self.blast_filter2 = Blastn(
            dir_specificity=self.dir_specificity2,
            word_size=self.config["primers_blast_setup"]['blast2_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )

        # third blast against 3' end of T7 promoter - Trim T7 to blast word size
        self.dir_specificity3 = os.path.join(self.dir_output, "specificity_temporary3")
        self.blast_filter3 = Blastn(
            dir_specificity=self.dir_specificity3,
            word_size=self.config["primers_blast_setup"]['blast3_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )

        T7_dict = dict()
        T7_dict['T7promoter'] = {
            "sequence": Seq(self.T7promoter[-self.config["primers_blast_setup"]['blast3_word_size']:])}
        # create reference DB with fasta file
        fasta_reference_database3 = self.blast_filter3._create_fasta_file(T7_dict, self.dir_specificity3, 'T7')
        self.reference_database3 = ReferenceDatabase(
            file_fasta=fasta_reference_database3)
        self.reference_database3.load_fasta_into_database()

        self.oligo_database = OligoDatabase(
            file_fasta=self.oligo_20nt_path,
            oligo_length_min=self.primer_oligo_config["oligo_length_min"],
            oligo_length_max=self.primer_oligo_config["oligo_length_max"],
            n_jobs=1,
            dir_output=self.dir_output,
        )
        self.oligo_database.create_database()


    def create_primer(self):
        property_filter = PropertyFilter(filters=self.filters,
                                         write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # property filter
        print("Property filter...Start")
        oligo_database = property_filter.apply(oligo_database=self.oligo_database, n_jobs=self.config["n_jobs"])
        print("Property filter...Done")

        print("Specifity filter 1...Start")
        # specifity filter 1
        specificity_filter1 = SpecificityFilter(filters=[self.blast_filter1],
                                                write_regions_with_insufficient_oligos=self.config[
                                                    "write_removed_genes"])
        print("Specifity filter 1...Apply")
        oligo_database = specificity_filter1.apply(oligo_database=oligo_database,
                                                   reference_database=self.reference_database1,
                                                   n_jobs=self.config["n_jobs"])
        print("Specifity filter 1...Done")
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="primer_database_specificity_filter_1.txt")

        # specificity filter 2
        print("Specifity filter 2...Start")
        # create reference db with trimmed primers

        trimmed_primers = {}

        primer_genes = list(oligo_database.database.keys())
        primer_oligo_ids = [list(oligo_database.database[gene].keys())[0] for gene in primer_genes]
        for gene, oligo_id in zip(primer_genes, primer_oligo_ids):
            # Get 3' end sequences
            trimmed_primers[gene] = str(oligo_database.database[gene][oligo_id]["sequence"][-self.config["primers_blast_setup"]['blast2_word_size']:])

        # create reference DB with fasta file

        fasta_reference_database2 = os.path.join(self.dir_specificity2,"oligos_primers.fna")

        with open(fasta_reference_database2, "w") as handle:
            for name, seq in trimmed_primers.items():
                handle.write(">" + name + "\n")
                handle.write(str(seq) + "\n")
            print("write trimmed primer done!")

        reference_database2 = ReferenceDatabase(
            file_fasta=fasta_reference_database2)
        reference_database2.load_fasta_into_database()
        specificity_filter2 = SpecificityFilter(filters=[self.blast_filter2],
                                                write_regions_with_insufficient_oligos=self.config[
                                                    "write_removed_genes"])
        oligo_database = specificity_filter2.apply(oligo_database=oligo_database,
                                                   reference_database=reference_database2, n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="primer_database_specificity_filter2.txt")
        print("Specifity filter 2...Done")
        # specificity filter 3
        print("Specifity filter 3...Start")
        specificity_filter3 = SpecificityFilter(filters=[self.blast_filter3],
                                                write_regions_with_insufficient_oligos=self.config[
                                                    "write_removed_genes"])
        oligo_database = specificity_filter3.apply(oligo_database=oligo_database,
                                                   reference_database=self.reference_database3,
                                                   n_jobs=self.config["n_jobs"])
        print("Specifity filter 3...Done")


        #save fasta file of primers
        primer_file_database = os.path.join(self.dir_output, "primers.fna")
        output = []
        
        keys = list(oligo_database.database.keys())[0:2*self.num_seq]
        for key in keys:
            oligo=list(oligo_database.database[key])[0]
            seq=oligo_database.database[key][oligo]["sequence"]
            output.append(SeqRecord(seq, key, "", ""))
         
        with open(primer_file_database, "w") as handle:
            SeqIO.write(output, handle, "fasta")

        print("Writing Primer1...Start")
        primer1_oligos_dict = {}
        primer1_genes = list(oligo_database.database.keys())[0:self.num_seq]
        primer1_oligo_ids = [list(oligo_database.database[gene].keys())[0] for gene in primer1_genes]
        for gene, oligo_id in zip(primer1_genes, primer1_oligo_ids):
            primer1_oligos_dict[gene] = str(oligo_database.database[gene][oligo_id]["sequence"])
        print("Writing Primer1...Done")

        print("Writing Primer2...Start")
        primer2_oligos_dict = {}
        primer2_genes = list(oligo_database.database.keys())[self.num_seq + 1: (self.num_seq * 2)+1]
        primer2_oligo_ids = [list(oligo_database.database[gene].keys())[0] for gene in primer2_genes]
        for gene, oligo_id in zip(primer2_genes, primer2_oligo_ids):
            primer2_seq = str(oligo_database.database[gene][oligo_id]["sequence"].reverse_complement())
            primer2_seq = primer2_seq + self.T7promoter
            primer2_oligos_dict[gene] = primer2_seq
        print("Writing Primer2...Done")
        return list(primer1_oligos_dict.values()), list(
            primer2_oligos_dict.values()), primer_file_database  # maybe take half of them for primer1 half for primer2?
