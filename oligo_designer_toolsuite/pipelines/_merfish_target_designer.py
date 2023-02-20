import numpy as np
from oligo_designer_toolsuite.oligo_property_filter import (
    PropertyFilter,
    MaskedSequences,
    GCContent,
    MeltingTemperature,
    ConsecutiveRepeats,
    Secondary_struct
)
from Bio.SeqUtils import MeltingTemp as mt
from oligo_designer_toolsuite.database import CustomOligoDB, NcbiOligoDB, EnsemblOligoDB
from Bio.Seq import Seq
import os
import yaml

from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilter,
    ExactMatches,
    Blastn,
)

class TargetProbes:

    def __init__(
            self,
            config_file
    ):
        self.config_file= config_file
        with open(self.config_file, 'r') as yaml_file:
            self.config = yaml.safe_load(yaml_file)
        self.dir_output = os.path.join(os.path.dirname(os.getcwd()), self.config["dir_output"]) 

    def create_target(self):

        # define the database class
        if self.config["source"] == "ncbi":
            # dowload the fasta files formthe NCBI server
            oligo_database = NcbiOligoDB(
                oligo_length_min=self.config["oligo_length_min"],
                oligo_length_max=self.config["oligo_length_max"],
                species=self.config["species"],
                annotation_release=self.config["annotation_release"],
                n_jobs=self.config["n_jobs"],
                dir_output=self.dir_output,
                min_oligos_per_gene=self.config["min_oligos_per_gene"],
                )
        elif self.config["source"] == "ensembl":
            # dowload the fasta files formthe NCBI server
            oligo_database = EnsemblOligoDB(
                oligo_length_min=self.config["oligo_length_min"],
                oligo_length_max=self.config["oligo_length_max"],
                species=self.config["species"],
                annotation_release=self.config["annotation_release"],
                n_jobs=self.config["n_jobs"],
                dir_output=self.dir_output,
                min_oligos_per_gene=self.config["min_oligos_per_gene"],
                )
        elif self.config["source"] == "custom":
            # use already dowloaded files
            oligo_database = CustomOligoDB(
                oligo_length_min=self.config["oligo_length_min"],
                oligo_length_max=self.config["oligo_length_max"],
                species=self.config["species"],
                genome_assembly=self.config["genome_assembly"],
                annotation_release=self.config["annotation_release"],
                files_source=self.config["files_source"],
                annotation_file=self.config["annotation_file"],
                sequence_file=self.config["sequence_file"],
                n_jobs=self.config["n_jobs"],
                dir_output=self.dir_output,
                min_oligos_per_gene=self.config["min_oligos_per_gene"],
                )
        else:
            raise ValueError("Annotation source not supported!") 

        # read the genes file
        with open(self.config["file_genes"]) as handle:
            lines = handle.readlines()
            genes = [line.rstrip() for line in lines]
            
        #generate the oligo sequences from gene transcripts
        oligo_database.create_oligos_DB(genes=genes, region='transcripts')

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
        melting_temperature = MeltingTemperature(
            Tm_min=self.config["targets_setup"]["Tm_min"], 
            Tm_max=self.config["targets_setup"]["Tm_max"], 
            Tm_parameters=Tm_params, 
            Tm_correction_parameters=Tm_correction_param
        )
        consecutive_repeats = ConsecutiveRepeats(self.config["targets_setup"]["max_repeats_AA"])
        gc_content = GCContent(
            GC_content_min=self.config["targets_setup"]["GC_content_min"], 
            GC_content_max=self.config["targets_setup"]["GC_content_max"]    
        )
        secondary_structure =Secondary_struct(
            T=self.config["targets_setup"]["internal_secondary_structures_T"],
             DG=self.config["targets_setup"]["internal_secondary_structures_threshold_deltaG"]
        )
        # create the list of filters
        filters = [gc_content, melting_temperature, consecutive_repeats,secondary_structure]

        # initialize the property filter class
        property_filter = PropertyFilter(filters=filters, write_genes_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        oligo_database = property_filter.apply(oligo_database=oligo_database, n_jobs=self.config["n_jobs"])
        # write the intermediate result in a file
        if self.config["write_intermediate_steps"]:
            oligo_database.write_oligos_DB(format=self.config["file_format"], dir_oligos_DB="property_filter")
        
        #Specificity filters to remove probes with more than 1 RNA species target (no cross-hybridization targets with Tm>72)
        dir_specificity = os.path.join(self.dir_output, "specificity_temporary") # folder where the temporary files will be written

        # generate the reference
        reference_database = CustomReferenceDB(
            species=oligo_database.species,
            genome_assembly=oligo_database.genome_assembly,
            annotation_release=oligo_database.annotation_release,
            files_source=oligo_database.files_source,
            annotation_file=oligo_database.annotation_file,
            sequence_file=oligo_database.sequence_file,
            dir_output=self.dir_output
        )
        reference_database.create_reference_DB(block_size=self.config["block_size"]) # leave the standard parameters
        # filter reference database by melting temperature
        melting_temperature2 = MeltingTemperature(
            Tm_min=self.config["targets_setup"]["cross_hybridization_targets_Tm_min"], 
            Tm_max=self.config["targets_setup"]["cross_hybridization_targets_Tm_max"], 
            Tm_parameters=Tm_params, 
            Tm_correction_parameters=Tm_correction_param
        )
        property_filter2 = PropertyFilter(filters=[melting_temperature2], write_genes_with_insufficient_oligos=self.config["write_removed_genes"])
        reference_database = property_filter2.apply(oligo_database=reference_database, n_jobs=self.config["n_jobs"])

        # intialize the filter classes
        exact_matches = ExactMatches(dir_specificity=dir_specificity)
        blastn = Blastn(
            dir_specificity=dir_specificity, 
            word_size=self.config["targeting_sequences_setup"]["word_size"],
            percent_identity=self.config["targeting_sequences_setup"]["percent_identity"],
            coverage=self.config["coverage"]
        )
        filters = [exact_mathces, blastn]
        # initialize the specificity filter class
        specificity_filter = SpecificityFilter(filters=filters, write_genes_with_insufficient_oligos=self.config["write_removed_genes"])
        # filte r the database
        oligo_database = specificity_filter.apply(oligo_database=oligo_database, reference_database=reference_database, n_jobs=self.config["n_jobs"])
        # write the intermediate result
        if self.config["write_intermediate_steps"]:
            oligo_database.write_oligos_DB(format=self.config["file_format"], dir_oligos_DB="specificity_filter")

        

        
        

        #write the final results
        oligo_database.write_oligos_DB(format=self.config["file_format"], dir_oligos_DB="target_sequences")

        # return  target probes 
        return oligo_database




