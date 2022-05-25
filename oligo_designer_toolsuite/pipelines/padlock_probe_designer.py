############################################
# imports
############################################

import os
import time
import argparse
from pathlib import Path

import logging
logger = logging.getLogger()
logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler()])

from datetime import datetime
timestamp = datetime.now()

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

    annotations = AnnotationLoader(config, dir_output)
    annotations.load_annotations()

    if not download_only:
        annotations.load_genes()
        annotations.load_transcriptome()
        annotations.load_probes()
    
    return annotations

    
def filter_probes(config, datamodule, dir_output, dir_annotations = None):

    probefilter = ProbeFilter(config, dir_output, datamodule.file_transcriptome_fasta, datamodule.genes, dir_annotations)
    probefilter.filter_probes_by_exactmatch()
    probefilter.run_blast_search()
    probefilter.filter_probes_by_blast_results()


def generate_probe_sets(config, dir_output, dir_probes = None):

    probesets = ProbesetsGenerator(config, dir_output, dir_probes)
    probesets.get_overlap_matrix()
    probesets.get_probe_sets(n_sets=100)


def design_padlock_probes(config, dir_output, dir_probes = None, dir_probesets = None):

    sequencedesigner = ProbeSequenceDesigner(config, dir_output, dir_probes, dir_probesets)
    sequencedesigner.design_padlocks()


def main():
    '''Main function of probe designer.
    '''
    # setup logging for standalone pipeline
    fh = logging.FileHandler(f'log_padlock_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt')
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)

    # get comman line arguments
    parameters = args()
    config = get_config(parameters.config)
    print_config(config)
    download_only = parameters.download_only


    logging.info('#########Start Pipeline#########')

    ### Create output directories
    dir_output = os.path.join(parameters.output)
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    logging.info('Results will be saved to: {}'.format(dir_output))

    ### Load annotations and probe lists
    t = time.time()

    annotations = download_annotations(config, dir_output, download_only)
    logging.info('Download gene annotation from {} and save as gene gtf: {}'.format(config['annotation_source'], annotations.file_gene_gtf))
    logging.info('Downloaded genome annotation from {} and save as genome fasta: {}'.format(config['annotation_source'], annotations.file_genome_fasta))
    
    if download_only:
        logging.info('Download finished. Pipeline interrupted. Set download_only to False to run the full pipeline.')
        return

    if config['file_genes'] is None:
        logging.info('Loaded gene list from {} annotation.'.format(config['annotation_source']))
    else:
        logging.info('Loaded gene list from {}.'.format(config['file_genes']))

    logging.info('Design probes for {} genes.'.format(len(annotations.genes)))

    t = (time.time() - t)/60
    logging.info('Time to load annotations, genes, transcriptome and probes: {} min'.format(t))


    ### Specificity Filter
    t = time.time()

    filter_probes(config, annotations, dir_output)
    del annotations # free memory

    t = (time.time() - t)/60
    logging.info('Time to filter with extact matches and blastn results: {} min'.format(t))


    ### Get ranked independent probe sets
    t = time.time()

    generate_probe_sets(config, dir_output)

    t = (time.time() - t)/60
    logging.info('Time to find nonoverlapping probe sets: {} min'.format(t))


    ### Design padlock probes
    t = time.time()

    design_padlock_probes(config, dir_output)

    t = (time.time() - t)/60
    logging.info('Time to design padlock probes: {} min'.format(t))
    

    logging.info('#########End Pipeline#########')


if __name__ == "__main__":
    main() 
