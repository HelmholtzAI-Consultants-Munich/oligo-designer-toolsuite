############################################
# imports
############################################

import os
import time
from pathlib import Path

from oligo_designer_toolsuite.utils import get_config, print_config, rm_intermediate_files
from oligo_designer_toolsuite.datamodule import DataModule
from oligo_designer_toolsuite.probefilter import ProbeFilter
from oligo_designer_toolsuite.get_overlap_matrix import get_overlap_matrix
from oligo_designer_toolsuite.nonoverlapping_sets import get_nonoverlapping_sets
from oligo_designer_toolsuite.padlock_design import design_padlocks



############################################
# functions
############################################

def download_and_filter(config, dir_output, logging, download_only=False):
    dir_output = os.path.join(dir_output, "")
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    logging.info('Results will be saved to: {}'.format(dir_output))

    datamodule = DataModule(config, logging, dir_output)
    
    t = time.time()
    datamodule.load_annotations()
    if download_only:
        logging.info("Download finished. Pipeline interrupted. Set download_only to False to run the full pipeline.")
        return
    datamodule.load_genes()
    datamodule.load_transcriptome()
    datamodule.load_probes()
    t = (time.time() - t)/60

    probefilter = ProbeFilter(config, logging, dir_output, datamodule.file_transcriptome_fasta, datamodule.genes)
    del datamodule # free memory

    logging.info('Time to load annotations, genes, transcriptome and probes: {} min'.format(t))
    print('Time to load annotations, genes, transcriptome and probes: {} min \n'.format(t))

    t = time.time()
    probefilter.filter_probes_by_exactmatch()
    t = (time.time() - t)/60

    logging.info('Time to filter with extact matches: {} min'.format(t))
    print('Time to filter with extact matches: {} min \n'.format(t))

    t = time.time()
    probefilter.run_blast_search()
    t = (time.time() - t)/60

    logging.info('Time to run Blast search: {} min'.format(t))
    print('Time to run Blast search: {} min \n'.format(t))

    t = time.time()
    probefilter.filter_probes_by_blast_results()
    t = (time.time() - t)/60

    logging.info('Time to filter with Blast results: {} min'.format(t))
    print('Time to filter with Blast results: {} min \n'.format(t))

    probefilter.log_statistics()

def find_non_overlapping_sets(config, dir_output, logging):
    t = time.time()

    get_overlap_matrix(os.path.join(dir_output,"probes"),os.path.join(dir_output,"overlap"))
    t = (time.time() - t)/60    
    
    logging.info('Time to compute overlap matrices: {} min'.format(t))
    print('Time to compute overlap matrices: {} min \n'.format(t))
    
    t = time.time()

    get_nonoverlapping_sets(config,os.path.join(dir_output,"probes"),os.path.join(dir_output,"overlap"),
                            os.path.join(dir_output,"probesets"),n_sets=100)
    t = (time.time() - t)/60

    logging.info('Time to find nonoverlapping probe sets: {} min'.format(t))
    print('Time to find nonoverlapping probe sets: {} min \n'.format(t))

def design_padlock_probes(config, dir_output, logging):    
    t = time.time()    
    design_padlocks(config,os.path.join(dir_output,"probes"),os.path.join(dir_output,"probesets"),
                    os.path.join(dir_output,"padlock_probes"))
    t = (time.time() - t)/60

    logging.info('Time to design padlock probes: {} min'.format(t))
    print('Time to design padlock probes: {} min \n'.format(t))    


