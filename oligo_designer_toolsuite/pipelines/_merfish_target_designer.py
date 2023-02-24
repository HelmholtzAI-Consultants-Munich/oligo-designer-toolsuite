import warnings

import numpy as np
from oligo_designer_toolsuite.oligo_property_filter import (
    PropertyFilter,
    MaskedSequences,
    GCContent,
    MeltingTemperatureNN,
    ConsecutiveRepeats,
    Secondary_struct
)
from Bio.SeqUtils import MeltingTemp as mt
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.database import ReferenceDatabase
from oligo_designer_toolsuite.database import CustomGenomicRegionGenerator, NcbiGenomicRegionGenerator, EnsemblGenomicRegionGenerator
from Bio.Seq import Seq
import os
import yaml
from pathlib import Path

from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilter,
    ExactMatches,
    Blastn,
)
from oligo_designer_toolsuite.oligo_efficiency import(
    PadlockOligoScoring,
    PadlockSetScoring,
)
from oligo_designer_toolsuite.oligo_selection import OligosetGenerator, padlock_heuristic_selection

class TargetProbes:

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
        property_filter = PropertyFilter(filters=filters, write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        oligo_database = property_filter.apply(oligo_database=oligo_database, n_jobs=self.config["n_jobs"])
        # write the intermediate result in a file
        # write the intermediate result in a file
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="oligo_database_property_filter.txt")
 
        
        #Specificity filters to remove probes with more than 1 RNA species target (no cross-hybridization targets with Tm>72)
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

        # filter reference database by melting temperature
        melting_temperature2 = MeltingTemperatureNN(
            Tm_min=self.config["targets_setup"]["cross_hybridization_targets_Tm_min"], 
            Tm_max=self.config["targets_setup"]["cross_hybridization_targets_Tm_max"], 
            Tm_parameters=Tm_params, 
            Tm_chem_correction_parameters=Tm_correction_param
        )
        property_filter2 = PropertyFilter(filters=[melting_temperature2], write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        #reference_database = property_filter2.apply(oligo_database=reference_database, n_jobs=self.config["n_jobs"])
        #NOTE: property filters don' work for reference database

        # intialize the filter classes
        exact_matches = ExactMatches(dir_specificity=dir_specificity)
        blastn = Blastn(
            dir_specificity=dir_specificity, 
            word_size=self.config["targeting_sequences_setup"]["word_size"],
            percent_identity=self.config["targeting_sequences_setup"]["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        filters = [exact_matches, blastn]
        # initialize the specificity filter class
        specificity_filter = SpecificityFilter(filters=filters, write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        # filter the database
        oligo_database = specificity_filter.apply(oligo_database=oligo_database, reference_database=reference_database, n_jobs=self.config["n_jobs"])
        

        # # write the intermediate result
        # if self.config["write_intermediate_steps"]:
        #     oligo_database.write_database(filename="merfish_target_probes.txt")

        # # initialize the scoring classes
        # oligos_scoring = PadlockOligoScoring(
        #     Tm_min=self.config["targets_setup"]["Tm_min"],
        #     Tm_opt=self.config["Tm_opt"],
        #     Tm_max=self.config["targets_setup"]["Tm_max"],
        #     GC_content_min=self.config["targets_setup"]["GC_content_min"],
        #     GC_content_opt=self.config["GC_content_opt"],
        #     GC_content_max=self.config["targets_setup"]["GC_content_max"],
        #     Tm_weight=self.config["Tm_weight"],
        #     GC_weight=self.config["GC_weight"],
        # )
        # set_scoring = PadlockSetScoring()

        # # initialize the oligoset generator class
        # oligoset_generator = OligosetGenerator(
        #     oligoset_size=self.config["oligoset_size"], 
        #     min_oligoset_size=self.config["min_oligoset_size"],
        #     oligos_scoring=oligos_scoring,
        #     set_scoring=set_scoring,
        #     heurustic_selection=padlock_heuristic_selection,
        #     write_regions_with_insufficient_oligos=self.config["write_removed_genes"]
        # )

        # # generate the oligoset
        # oligo_database = oligoset_generator.apply(oligo_database=oligo_database, n_sets=self.config["n_sets"], n_jobs=self.config["n_jobs"])

        # write the result

        file_database = oligo_database.write_database(filename="merfish_target_probes.txt")
    

        # return  target probes 
        return oligo_database, file_database




