############################################
# imports
############################################

import os
import time
import argparse
import logging
from pathlib import Path
from datetime import datetime

timestamp = datetime.now()
file_logger = f'log_padlock_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt'
logging.getLogger('padlock_probe_designer')
logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s", level=logging.NOTSET, handlers=[logging.FileHandler(file_logger), logging.StreamHandler()])

from oligo_designer_toolsuite.utils import get_config, print_config
from oligo_designer_toolsuite.annotation_loader import AnnotationLoader
from oligo_designer_toolsuite.probe_filter import ProbeFilter
from oligo_designer_toolsuite.probesets_generator import ProbesetsGenerator
from oligo_designer_toolsuite.probe_sequence_designer import ProbeSequenceDesigner


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


def download_annotations(config, dir_output, download_only=False):

    t = time.time()

    annotations = AnnotationLoader(config, dir_output)
    annotations.load_annotations()

    if not download_only:
        annotations.load_genes()
        annotations.load_transcriptome()
        annotations.load_probes()
    
    t = (time.time() - t)/60
    logging.info('Time to load annotations, genes, transcriptome and probes: {} min'.format(t))
    
    return annotations

    
def filter_probes(config, datamodule, dir_output, dir_annotations = None):

    probefilter = ProbeFilter(config, dir_output, datamodule.file_transcriptome_fasta, datamodule.genes, dir_annotations)
    
    t = time.time()
    probefilter.filter_probes_by_exactmatch()
    t = (time.time() - t)/60
    logging.info('Time to filter with extact matches: {} min'.format(t))

    t = time.time()
    probefilter.run_blast_search()
    t = (time.time() - t)/60
    logging.info('Time to run Blast search: {} min'.format(t))

    t = time.time()
    probefilter.filter_probes_by_blast_results()
    t = (time.time() - t)/60
    logging.info('Time to filter with Blast results: {} min'.format(t))


def generate_probe_sets(config, dir_output, dir_probes = None):

    probesets = ProbesetsGenerator(config, dir_output, dir_probes)

    t = time.time()
    probesets.get_overlap_matrix()
    t = (time.time() - t)/60    
    logging.info('Time to compute overlap matrices: {} min'.format(t))
    
    t = time.time()
    probesets.get_probe_sets(n_sets=100)
    t = (time.time() - t)/60
    logging.info('Time to find nonoverlapping probe sets: {} min'.format(t))


def design_padlock_probes(config, dir_output, dir_probes = None, dir_probesets = None):  

    sequencedesigner = ProbeSequenceDesigner(config, dir_output, dir_probes, dir_probesets)
    
    t = time.time()   
    sequencedesigner.design_padlocks()
    t = (time.time() - t)/60
    logging.info('Time to design padlock probes: {} min'.format(t))


def main():
    '''Main function of probe designer.
    '''
    # get comman line arguments
    parameters = args()

    dir_output = os.path.join(parameters.output)
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    logging.info('Results will be saved to: {}'.format(dir_output))

    config = get_config(parameters.config)
    print_config(config, logging)

    download_only = parameters.download_only

    logging.info('#########Start Pipeline#########')
    t_pipeline = time.time()

    annotations = download_annotations(config, dir_output, download_only)
    if download_only:
        logging.info("Download finished. Pipeline interrupted. Set download_only to False to run the full pipeline.")
        return

    filter_probes(config, annotations, dir_output)
    del annotations # free memory

    generate_probe_sets(config, dir_output)
    design_padlock_probes(config, dir_output)
    
    t_pipeline = (time.time() - t_pipeline)/60
    logging.info('Time Pipeline: {} min'.format(t_pipeline))
    logging.info('#########End Pipeline#########')
