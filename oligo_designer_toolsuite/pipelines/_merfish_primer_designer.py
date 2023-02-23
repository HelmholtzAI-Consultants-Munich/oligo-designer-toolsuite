import numpy as np
from Bio import SeqIO
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
from oligo_designer_toolsuite.database import ReferenceDatabase

from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilter,
    Blastn,
)


class PrimerProbes:

    def __init__(
            self,
            config_path,
            file_transcriptome,
            region_generator
          
    ):
        self.length = 25  # need to trim the sequece
        #self.config_file = os.path.join(os.getcwd(), "tutorials", "configs", "merfish_probe_designer_test.yaml")
        with open(config_path, 'r') as yaml_file:
            self.config = yaml.safe_load(yaml_file)
        self.file_transcriptome = file_transcriptome
        self.region_generator= region_generator
        self.config_param = self.config["primers_setup"].copy()
        self.GC_min = self.config_param["GC_content_min"]
        self.GC_max = self.config_param["GC_content_max"]
        self.GC_content_filter = GCContent(GC_content_min=self.GC_min, GC_content_max=self.GC_max)
        self.repeat_num_max = ConsecutiveRepeats(3)
        self.lib = ['A', 'C', 'G', 'T']
        self.T7promoter = "TAATACGACTCACTATAG"
        #self.blast_filter = blast_filter
        #self.ref = os.path.basename(reference_DB.file_reference_DB)
        self.oligo_25nt_path = os.path.join(os.getcwd(),"oligo_designer_toolsuite","experiment_specific","data", "bc25mer.240k.fasta")
        self.oligo_25nt_dict = SeqIO.to_dict(SeqIO.parse(self.oligo_25nt_path, "fasta"))
        self.oligo_20nt_dict = {k: v[:-5] for k, v in self.oligo_25nt_dict.items()}  # step 1


        self.oligo_20nt_DB = None




        self.Tm_params = self.config["Tm_parameters"]["shared"].copy()
        self.Tm_correction_param = self.config["Tm_correction_parameters"]["shared"].copy()
        self.melting_temperature = MeltingTemperatureNN(
            #min_arm_length=self.config["primers_setup"]["min_arm_length"],
            #max_arm_Tm_dif=self.config["primers_setup"]["max_arm_Tm_dif"],
            Tm_max=self.config_param["Tm_max"],
            Tm_min=self.config_param["Tm_min"],
            Tm_parameters=self.Tm_params,
            Tm_correction_parameters=self.Tm_correction_param,
        )

        self.consecutive_repeats = ConsecutiveRepeats(self.config_param["Repeat_AA_max"])
        self.GC_clamp = GCClamp(self.config_param["GC_clamp_n"])

        # create the list of filters
        self.filters = [self.GC_content_filter, self.melting_temperature, self.consecutive_repeats, self.GC_clamp]


        self.oligo_primer_DB = None
        # specificty filters

        

        #first blast against human transcriptome
        self.dir_specificity1 = os.path.join(self.config["dir_output"], "specificity_temporary1")
        self.blast_filter1 = Blastn(
            dir_specificity= self.dir_specificity1,
            word_size=self.config["primers_blast_setup"]['blast1_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
        )
        self.reference_database1= ReferenceDatabase(
            file_fasta = self.file_transcriptome,
            files_source = self.region_generator.files_source,
            species = self.region_generator.species,
            annotation_release = self.region_generator.annotation_release,
            genome_assembly = self.region_generator.genome_assembly,
            dir_output=self.dir_output
        )
        self.reference_database1.load_fasta_into_database()


        #second blast against 3' end of other primers
        self.dir_specificity2 = os.path.join(self.config["dir_output"], "specificity_temporary2")
        self.blast_filter2 = Blastn(
            dir_specificity= self.dir_specificity2,
            word_size=self.config["primers_blast_setup"]['blast2_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
        )
        
        

        #third blast against 3' end of T7 promoter - Trim T7 to blast word size
        self.dir_specificity3 = os.path.join(self.config["dir_output"], "specificity_temporary3")
        self.blast_filter3 = Blastn(
            dir_specificity= self.dir_specificity3,
            word_size=self.config["primers_blast_setup"]['blast3_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
        )
        
        T7_dict=dict()
        T7_dict['T7promoter'] = {"sequence": Seq(self.T7promoter[-self.config["primers_blast_setup"]['blast3_word_size']:])}
        #create reference DB with fasta file
        fasta_reference_database3 = self.blast_filter2._create_fasta_file(T7_dict, self.dir_specificity3, 'T7')
        self.reference_database3= ReferenceDatabase(
            file_fasta = fasta_reference_database3)
        self.reference_database3.load_fasta_into_database()
        

    def create_primer1(self):
        primer1 = list()

        property_filter = PropertyFilter(filters=self.filters,
                                         write_genes_with_insufficient_oligos=self.config["write_removed_genes"])
        # property filter
        oligo_database = property_filter.apply(oligo_database=self.oligo_20nt_dict, n_jobs=self.config["n_jobs"])

        
        #specifity filter 1
        specificity_filter1 = SpecificityFilter(filters=[self.blast_filter1], write_genes_with_insufficient_oligos=self.config["write_removed_genes"])
        oligo_database = specificity_filter1.apply(oligo_database=oligo_database, reference_database=self.reference_database1, n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="primer_database_specificity_filter_1.txt")

        # specificity filter 2
            #create reference db with trimmed primers
        trimmed_primers=  {k: v[-self.config["primers_blast_setup"]['blast2_word_size']:] for k, v in oligo_database.items()}
            # create reference DB with fasta file
        fasta_reference_database2 = self.blast_filter2._create_fasta_file(trimmed_primers, self.dir_specificity2, 'primers')
        reference_database2= ReferenceDatabase(
            file_fasta = fasta_reference_database2)
        reference_database2.load_fasta_into_database()
        specificity_filter2 = SpecificityFilter(filters=[self.blast_filter2], write_genes_with_insufficient_oligos=self.config["write_removed_genes"])
        oligo_database = specificity_filter2.apply(oligo_database=oligo_database, reference_database=reference_database2, n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="primer_database_specificity_filter2.txt")

        #specificity filter 3
        specificity_filter3 = SpecificityFilter(filters=[self.blast_filter3], write_genes_with_insufficient_oligos=self.config["write_removed_genes"])
        oligo_database = specificity_filter3.apply(oligo_database=oligo_database, reference_database=self.reference_database3, n_jobs=self.config["n_jobs"])
        
        
        primer1_file_database = oligo_database.write_database(filename="primer1_database.txt")
        
        

        

        return oligo_database[1:self.num_probes], primer1_file_database

    def create_primer2(self):
        pass

    def creat_primer_probs(self):
        pass
