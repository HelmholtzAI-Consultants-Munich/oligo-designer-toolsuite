import os
import random
import warnings
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import (CustomGenomicRegionGenerator,
                                               EnsemblGenomicRegionGenerator,
                                               NcbiGenomicRegionGenerator,
                                               OligoDatabase,
                                               ReferenceDatabase)
from oligo_designer_toolsuite.oligo_property_filter import (
    ConsecutiveRepeats, GCClamp, GCContent, HairpinSecondaryStructure,
    MaskedSequences, MeltingTemperatureNN, PadlockArms, PropertyFilter)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGenerator, padlock_heuristic_selection)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    Blastn, ExactMatches, SpecificityFilter)


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
            primer2_seq = self.T7promoter + primer2_seq
            primer2_oligos_dict[gene] = primer2_seq
        print("Writing Primer2...Done")
        return list(primer1_oligos_dict.values()), list(
            primer2_oligos_dict.values()), primer_file_database  # maybe take half of them for primer1 half for primer2?


class TargetProbes:
    """
    This class is used to design the target probes.

    :param config: config file
    :type config: file pointer
    :param dir_output: output directory
    :type dir_output: str
    :param file_transcriptome: directory of the fasta file for the transcriptome
    :type file_transcriptome: str
    :param region_generator: region generator used to create the file_transcriptome
    :type region_generator: CustomGenomicRegionGenerator
    """

    def __init__(
            self,
            config,
            dir_output,
            file_transcriptome,
            region_generator
    ):
        self.config= config
        self.dir_output = os.path.join(dir_output, "target_probes")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)
        self.file_transcriptome = file_transcriptome
        self.region_generator = region_generator


    def create_target(self):

        # creating the database
        # define the database class
        oligo_database = OligoDatabase(
            file_fasta = self.file_transcriptome,
            oligo_length_min = self.config["target_oligo"]["oligo_length_min"],
            oligo_length_max = self.config["target_oligo"]["oligo_length_max"],
            min_oligos_per_region = self.config["min_oligos_per_gene"],
            files_source = self.region_generator.files_source,
            species = self.region_generator.species,
            annotation_release = self.region_generator.annotation_release,
            genome_assembly = self.region_generator.genome_assembly,
            n_jobs = self.config["n_jobs"],
            dir_output=self.dir_output
        )

        # read the genes file
        if self.config["file_genes"] is None:
            warnings.warn(
                "No file containing the genes was provided, all the genes are ussed to generate the probes. This chioce can use a lot of resources."
            )
            genes = None
        else:
            with open(self.config["file_genes"]) as handle:
                lines = handle.readlines()
                genes = [line.rstrip() for line in lines]
                
        # generate the oligo sequences from gene transcripts
        oligo_database.create_database(region_ids=genes) 




        #Property filters
        
        # the melting temperature params need to be preprocessed
        Tm_params = self.config["Tm_parameters"]["shared"].copy()
        Tm_params.update(self.config["Tm_parameters"]["property_filter"])
        Tm_params["nn_table"] = getattr(mt, Tm_params["nn_table"])
        Tm_params["tmm_table"] = getattr(mt, Tm_params["tmm_table"])
        Tm_params["imm_table"] = getattr(mt, Tm_params["imm_table"])
        Tm_params["de_table"] = getattr(mt, Tm_params["de_table"])

        Tm_correction_param = self.config["Tm_correction_parameters"]["shared"].copy()
        Tm_correction_param.update(self.config["Tm_correction_parameters"]["property_filter"])
        melting_temperature = MeltingTemperatureNN(
            Tm_min=self.config["targets_setup"]["Tm_min"], 
            Tm_max=self.config["targets_setup"]["Tm_max"], 
            Tm_parameters=Tm_params, 
            Tm_chem_correction_parameters=Tm_correction_param
        )
        consecutive_repeats = ConsecutiveRepeats(self.config["targets_setup"]["max_repeats_nt"])
        gc_content = GCContent(
            GC_content_min=self.config["targets_setup"]["GC_content_min"], 
            GC_content_max=self.config["targets_setup"]["GC_content_max"]    
        )
        secondary_structure = HairpinSecondaryStructure(
            T=self.config["targets_setup"]["internal_secondary_structures_T"],
             DG=self.config["targets_setup"]["internal_secondary_structures_threshold_deltaG"]
        )
        # create the list of filters
        filters = [gc_content, melting_temperature, consecutive_repeats,secondary_structure]

        # initialize the property filter class
        property_filter = PropertyFilter(filters=filters, write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        oligo_database = property_filter.apply(oligo_database=oligo_database, n_jobs=self.config["n_jobs"])
        # write the intermediate result in a file
        # write the intermediate result in a file
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="oligo_database_property_filter.txt")
 
        
        #Specificity filters to remove probes with more than 1 RNA species target 
        dir_specificity = os.path.join(self.dir_output, "specificity_temporary") # folder where the temporary files will be written

        reference_database = ReferenceDatabase(
            file_fasta = self.file_transcriptome,
            files_source = self.region_generator.files_source,
            species = self.region_generator.species,
            annotation_release = self.region_generator.annotation_release,
            genome_assembly = self.region_generator.genome_assembly,
            dir_output=self.dir_output
        )
        reference_database.load_fasta_into_database()

        # intialize the filter classes
        exact_matches = ExactMatches(dir_specificity=dir_specificity)
        blastn = Blastn(
            dir_specificity=dir_specificity, 
            word_size=self.config["targeting_sequences_setup"]["word_size"],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        filters = [exact_matches, blastn]
        # initialize the specificity filter class
        specificity_filter = SpecificityFilter(filters=filters, write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        oligo_database = specificity_filter.apply(oligo_database=oligo_database, reference_database=reference_database, n_jobs=self.config["n_jobs"])




        #Specificity filter to remove cross hybridization targets
        targets_fasta= oligo_database.write_fasta_from_database(filename = 'target_probes_init')
        reference_database2 = ReferenceDatabase(file_fasta =targets_fasta)
        reference_database2.load_fasta_into_database()
        dir_specificity2 = os.path.join(self.dir_output, "specificity_temporary2") # folder where the temporary files will be written
        # intialize the filter classes
        blastn2 = Blastn(
            dir_specificity=dir_specificity2, 
            word_size=self.config["targeting_sequences_setup"]["word_size"],
            percent_identity=self.config["targeting_sequences_setup"]["percent_identity_ch"],
            coverage=self.config["coverage"],
            strand='minus', 
        )
        filters2 = [blastn2]
        # initialize the specificity filter class
        specificity_filter = SpecificityFilter(filters=filters2, write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        oligo_database = specificity_filter.apply(oligo_database=oligo_database, reference_database=reference_database2, n_jobs=self.config["n_jobs"])



        # write the result

        file_database = oligo_database.write_database(filename="merfish_target_probes.txt")
    

        # return  target probes 
        return oligo_database, file_database




