############################################
# imports
############################################

from datetime import datetime
timestamp = datetime.now()

import logging
logging.basicConfig(filename='log_probe_design_{}-{}-{}-{}-{}.txt'.format(timestamp.year, timestamp.month, timestamp.day, timestamp.hour, timestamp.minute), level=logging.NOTSET)

import os
import time
import argparse

import src.utils as utils
import src.load_annotations as la
import src.get_probes as gp
import src.filter_probes as fp

############################################
# pipeline
############################################

def args():
    
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('-c','--config',help='path to config yaml file',type=str,required=True)
    return args_parser.parse_args()


############################################

def probe_pipeline(config):
    
    dir_output = config['dir_output']
    logging.info('Results will be saved to: {}'.format(dir_output))
    get_annotation = config['get_annotation']

    #############
    ### download annotations
    #############
    if get_annotation:
        source = config['source']
        logging.info('Download annotation from: {}'.format(source))
        
        if source == 'ncbi':
            # if source ncbi we need a chromosome mapping
            mapping = la.download_chr_mapping(dir_output)
        else:
            mapping = None
        
        if source == 'ensemble':
            release_ensemble = config['release_ensemble']
            logging.info('Using Ensemble release {}'.format(release_ensemble))
        else:
            release_ensemble = None

        print('Download gene annotation.')
        file_gene_gtf = la.download_gene_gtf(source, release_ensemble, mapping, dir_output)

        print('Download genome annotation.')
        file_genome_fasta = la.download_genome_fasta(source, release_ensemble, mapping, dir_output)
        
        logging.info('File gene gtf: {}'.format(file_gene_gtf))
        logging.info('File genome fasta: {}'.format(file_genome_fasta))
 
    else:
        file_gene_gtf = config['file_gene_gtf']
        file_genome_fasta = config['file_genome_fasta']

    #############
    ### get gene list
    #############
    print('Get gene list.')    
    exon_annotation = utils.load_exon_annotation(file_gene_gtf)

    if config['gene_list']:
        genes = config['gene_list']
    else:
        genes = utils.get_gene_list(exon_annotation)
    logging.info('Probes for {} genes will be designed.'.format(len(genes)))

    batch_size = config['batch_size']
    number_batchs = int(len(genes) / batch_size) + (len(genes) % batch_size > 0)
    logging.info('{} processes will run in parallel with {} genes in one batch'.format(number_batchs, batch_size))

    print('Get mapping.')
    t = time.time()
    utils.get_gene_transcrip_exon_mapping(exon_annotation, dir_output)
    logging.info('Time to get mapping: {} min'.format((time.time() - t)/60))
    

    #############
    ### get probes
    #############
    probe_length = config['probe_length']
    GC_content_min = config['GC_content_min']
    GC_content_max = config['GC_content_max']
    Tm_parameters = config['Tm_parameters']
    Tm_min = config['Tm_min']
    Tm_max = config['Tm_max']

    dir_output_probes = os.path.join(dir_output, 'probes/')
    if os.path.isdir(dir_output_probes) == False:
        os.mkdir(dir_output_probes)
    
    print('Get probes.')
    t = time.time()
    gp.get_probes(number_batchs, batch_size, exon_annotation, genes, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_probes)
    logging.info('Time to get probes: {} min'.format((time.time() - t)/60))


    #############
    ### filter probes
    #############
    num_threads_blast = 4
    word_size = config['word_size'] 
    coverage = config['coverage']
    percent_identity = config['percent_identity']
    min_probes_per_gene = config['min_probes_per_gene']


    print('Run Blast search.')
    t = time.time()
    fp.run_blast_search(number_batchs, word_size, num_threads_blast, file_gene_gtf, file_genome_fasta, dir_output_probes)
    logging.info('Time to run Blast search: {} min'.format((time.time() - t)/60))


    dir_output_results = os.path.join(dir_output, 'results/')
    if os.path.isdir(dir_output_results) == False:
        os.mkdir(dir_output_results)

    print('Filter probes with Blast results')
    t = time.time()
    fp.filter_probes_by_blast_results(number_batchs, probe_length, percent_identity, coverage, min_probes_per_gene, dir_output_probes, dir_output_results)
    logging.info('Time to filter with Blast results: {} min'.format((time.time() - t)/60))
    
    
############################################
    
if __name__ == '__main__':
    
    # get comman line arguments
    parameters = args()
    config = utils.get_config(parameters.config) 
    utils.print_config(config, logging)


    logging.info('#########Start Pipeline#########')

    t = time.time()
    probe_pipeline(config)
    logging.info('Time for whole pipeline: {} min'.format((time.time() - t)/60))

    logging.info('#########End Pipeline#########')
