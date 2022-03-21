############################################
# imports
############################################

import os
import time
import argparse
import logging
from datetime import datetime
from pathlib import Path

timestamp = datetime.now()
logging.basicConfig(filename='log_probe_design_{}-{}-{}-{}-{}.txt'.format(timestamp.year, timestamp.month, timestamp.day, timestamp.hour, timestamp.minute), level=logging.NOTSET)

from src.utils import get_config, print_config, rm_intermediate_files
from src.datamodule import DataModule
from src.probefilter import ProbeFilter
from src.get_overlap_matrix import get_overlap_matrix
from src.nonoverlapping_sets import get_nonoverlapping_sets
from src.padlock_design import design_padlocks


############################################
# functions
############################################

def args():
    '''Argument parser for command line arguments.

    :return: Command line arguments with their values.
    :rtype: dict
    '''
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('-c','--config',help='path to config yaml file',type=str,required=True)
    args_parser.add_argument('-o','--output',help='path of output folder',type=str,required=True)
    args_parser.add_argument('-d','--download_only',help='dry run that only downloads the necessary files',action='store_true',required=False)
    return args_parser.parse_args()

def probe_pipeline(config, dir_output, download_only=False):
    '''Pipeline of probe designer. Sets up all required directories; loads annotations, genes and probes; 
    and filters probes based on sequence properties and blast aligment search results.

    :param config: User-defined parameters, where keys are the parameter names and values are the paremeter values.
    :type config: dict
    :param dir_output: User-defined output directory.
    :type dir_output: string
    :param download_only: Whether to interupt the pipeline after download of the necessary files.
    :type download_only: bool
    '''
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
    
    t = time.time()
    get_overlap_matrix(os.path.join(dir_output,"probes"),os.path.join(dir_output,"overlap"))
    t = (time.time() - t)/60    
    
    logging.info('Time to compute overlap matrices: {} min'.format(t))
    print('Time to compute overlap matrices: {} min \n'.format(t))
    
    t = time.time()
    get_nonoverlapping_sets(config,os.path.join(dir_output,"probes"),os.path.join(dir_output,"overlap"),
                            os.path.join(dir_output,"probesets"),n_sets=100
                           )
    t = (time.time() - t)/60

    logging.info('Time to find nonoverlapping probe sets: {} min'.format(t))
    print('Time to find nonoverlapping probe sets: {} min \n'.format(t))
    
    t = time.time()    
    design_padlocks(config,os.path.join(dir_output,"probes"),os.path.join(dir_output,"probesets"),
                    os.path.join(dir_output,"padlock_probes"))
    t = (time.time() - t)/60

    logging.info('Time to design padlock probes: {} min'.format(t))
    print('Time to design padlock probes: {} min \n'.format(t))    


############################################
# main
############################################
    
if __name__ == '__main__':
    '''Main function of probe designer.
    '''
    # get comman line arguments
    parameters = args()

    dir_output = parameters.output
    config = get_config(parameters.config)
    download_only = parameters.download_only
    print_config(config, logging)

    logging.info('#########Start Pipeline#########')
    print('#########Start Pipeline#########')

    t_pipeline = time.time()
    probe_pipeline(config, dir_output, download_only=download_only)
    t_pipeline = (time.time() - t_pipeline)/60

    logging.info('Time Pipeline: {} min'.format(t_pipeline))
    logging.info('#########End Pipeline#########')

    print('Time Pipeline: {} min \n'.format(t_pipeline))
    print('#########End Pipeline#########')

    #rm_intermediate_files(dir_output)
