import os
import sys
from copy import deepcopy
import datetime
import logging
import time

from oligo_designer_toolsuite.database import *
from oligo_designer_toolsuite.oligo_property_filter import *
from oligo_designer_toolsuite.oligo_specificity_filter import *
import yaml
from Bio.SeqUtils import MeltingTemp as mt



# use both blast and bowtie + AI filter (with parmeters as similar as possible for comparison) and check which one filters less and how many oligos does the AI filter save.


### 0. define the logging file

timestamp = datetime.now()
file_logger = f"log_test_AI_filters_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt"
logging.getLogger("test_AI_filters")
logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
    handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
)


### 1. generate the oligos

logging.info("Start of the pipeline")
dir_package = os.path.dirname(os.getcwd())
sys.path.append(dir_package)

config_file = os.path.join(dir_package, "data/configs/test_AI_filters.yaml") # custom config


with open(config_file, 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)

dir_output = os.path.join(dir_package, config["dir_output"]) # create the complete path for the output directory
config["file_genes"] = os.path.join(dir_package, config["file_genes"])  # create the complete path for genes file

# only needed when custom config used
config["source_params"]["file_annotation"] = os.path.join(dir_package, config["source_params"]["file_annotation"]) 
config["source_params"]["file_sequence"] = os.path.join(dir_package, config["source_params"]["file_sequence"]) 

metadata = {"files_source" : config["source"], "source_params" : config["source_params"]}

# If the Ncbi config file is selected
region_generator = NcbiGenomicRegionGenerator(
    taxon=config["source_params"]["taxon"],
    species=config["source_params"]["species"], 
    annotation_release=config["source_params"]["annotation_release"], 
    dir_output=dir_output
)
file_transcriptome = region_generator.generate_transcript_reduced_representation(include_exon_junctions=True, exon_junction_size=2*config["probe_length_max"])

# define the database class
oligo_database = OligoDatabase(
    min_oligos_per_region = config["min_probes_per_gene"],
    metadata = metadata,
    n_jobs = 2,
    dir_output=dir_output
)

# read the genes file
with open(config["file_genes"]) as handle:
    lines = handle.readlines()
    genes = [line.rstrip() for line in lines]
    
#generate the oligo sequences from gene transcripts
oligo_database.create_database(file_fasta = file_transcriptome, oligo_length_min = config["probe_length_min"], oligo_length_max = config["probe_length_max"],region_ids = genes)
logging.info("Oligos generated")

### 2. property filters

# the melting temperature params need to be preprocessed
Tm_params = config["Tm_parameters_probe"]
Tm_params["nn_table"] = getattr(mt, Tm_params["nn_table"])
Tm_params["tmm_table"] = getattr(mt, Tm_params["tmm_table"])
Tm_params["imm_table"] = getattr(mt, Tm_params["imm_table"])
Tm_params["de_table"] = getattr(mt, Tm_params["de_table"])

Tm_chem_correction_param = config["Tm_chem_correction_param_probe"]

# initialize the filters clasees
masked_sequences = MaskedSequences()
gc_content = GCContent(GC_content_min=config["GC_content_min"], GC_content_max=config["GC_content_max"])
melting_temperature = MeltingTemperatureNN(
    Tm_min=config["Tm_min"], 
    Tm_max=config["Tm_max"], 
    Tm_parameters=Tm_params, 
    Tm_chem_correction_parameters=Tm_chem_correction_param
)
padlock_arms = PadlockArms(
    min_arm_length=config["min_arm_length"],
    max_arm_Tm_dif=config["max_arm_Tm_dif"],
    arm_Tm_min=config["arm_Tm_min"],
    arm_Tm_max=config["arm_Tm_max"],
    Tm_parameters=Tm_params,
    Tm_chem_correction_parameters=Tm_chem_correction_param,
)
# create the list of filters
filters = [masked_sequences, gc_content, melting_temperature, padlock_arms]

# initialize the property filter class
property_filter = PropertyFilter(filters=filters)
# filter the database
oligo_database = property_filter.apply(oligo_database=oligo_database, n_jobs=config["n_jobs"])
# write the intermediate result in a file
if config["write_intermediate_steps"]:
    file_database = oligo_database.write_database(filename="probe_database_property_filter.txt")
logging.info("Oligos filtered")


### 3. specificity filters (blat or bowtie + AI filters)
# TODO: what about a comparision without any ai filters (to understand the time overhead)
for ai_filter_thershold in config["ai_filter_thershold"]:
    start = time.time()
    oligo_database_2 = deepcopy(oligo_database) # apply the 2 specificity filters to 2 separate 
    dir_specificity = os.path.join(dir_output, "specificity_temporary") # folder where the temporary files will be written

    reference = ReferenceDatabase(
        file_fasta = file_transcriptome,
        metadata = metadata,
        dir_output=dir_output
        )

    # Blast filter
    exact_matches = ExactMatches(dir_specificity=dir_specificity)
    blastn = Blastn(
        dir_specificity=dir_specificity, 
        word_size=config["blast_word_size"],
        percent_identity=config["blast_percent_identity"],
        coverage=config["blast_coverage"],
        strand="plus",
        ai_filter="hybridization_probability",
        ai_filter_threshold=ai_filter_thershold
    )
    filters = [exact_matches, blastn]

    # initialize the specificity filter class
    specificity_filter = SpecificityFilter(filters=filters)
    # filter the database
    oligo_database = specificity_filter.apply(oligo_database=oligo_database, reference_database=reference, n_jobs=config["n_jobs"])

    logging.info(f"Blast Ai filtering done in {time.time() - start} seconds.")
    logging.info(f"The average saved oligos with threshold {ai_filter_thershold} is {blastn.filtered['Difference'].mean()}")
    #write hte dataframe 
    blastn.filtered.to_csv(os.path.join(config["output"], f"balstn_{ai_filter_thershold}"))


    #Bowtie filter
    start = time.time()
    # same as before, but with a the bowtie filter this time
    bowtie = Bowtie(dir_specificity=dir_specificity,
                    num_mismatches=config['num_mismatches'],
                    ai_filter="hybridization_probability",
                    ai_filter_threshold=ai_filter_thershold)

    filters = [exact_matches, bowtie]

    # initialize the specificity filter class
    specificity_filter = SpecificityFilter(filters=filters)
    # filter the database
    oligo_database2 = specificity_filter.apply(oligo_database=oligo_database2, reference_database=reference, n_jobs=config["n_jobs"])

    logging.info(f"Blast Ai filtering done in {time.time() - start} seconds.")
    logging.info(f"The average saved oligos with threshold {ai_filter_thershold} is {bowtie.filtered['Difference'].mean()}")
    #write hte dataframe 
    bowtie.filtered.to_csv(os.path.join(config["output"], f"bowtie_{ai_filter_thershold}"))
