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

def args():
    
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('-c','--config',help='path to config yaml file',type=str,required=True)
    args_parser.add_argument('-o','--output',help='path of output folder',type=str,required=True)
    return args_parser.parse_args()


############################################
# pipeline
############################################

def probe_pipeline(config, dir_output):
    
    logging.info('Results will be saved to: {}'.format(dir_output))
    get_annotation = config['get_annotation']


    #############
    ### download annotations
    #############

    dir_output_annotations = os.path.join(dir_output, 'annotations/')
    if os.path.isdir(dir_output_annotations) == False:
        os.mkdir(dir_output_annotations)

    if get_annotation:
        source = config['source']
        logging.info('Download annotation from: {}'.format(source))
        
        if source == 'ncbi':
            ftp_gene = config['ftp_gene']
            ftp_genome = config['ftp_genome']
            ftp_chr_mapping = config['ftp_chr_mapping']

            # if source ncbi we need a chromosome mapping
            mapping = la.download_chr_mapping(ftp_chr_mapping, dir_output_annotations)

        elif source == 'ensemble':
            ftp_gene = config['ftp_gene']
            ftp_genome = config['ftp_genome']

            # if source ensemble we don't need a chromosome mapping
            mapping = None

        else:
            raise ValueError('Error: unknown source "{}"'.format(source))

        file_gene_gtf = la.download_gene_gtf(source, mapping, ftp_gene, dir_output_annotations)
        file_genome_fasta = la.download_genome_fasta(source, mapping, ftp_genome, dir_output_annotations)
        
        logging.info('File gene gtf: {}'.format(file_gene_gtf))
        logging.info('File genome fasta: {}'.format(file_genome_fasta))
 
    else:
        file_gene_gtf = config['file_gene_gtf']
        file_genome_fasta = config['file_genome_fasta']


    #############
    ### get gene list
    #############
 
    exon_annotation = utils.load_exon_annotation(file_gene_gtf)

    if config['file_genes']:
        file_genes = config['file_genes']
        genes = utils.read_gene_list(file_genes)
        print(genes)
    else:
        genes = utils.get_gene_list(exon_annotation)
    logging.info('Probes for {} genes will be designed.'.format(len(genes)))

    number_batchs = config['number_batchs']
    batch_size = int(len(genes) / number_batchs) + (len(genes) % number_batchs > 0)

    logging.info('{} processes will run in parallel with {} genes in one batch'.format(number_batchs, batch_size))


    #############
    ### get probes
    #############

    probe_length = config['probe_length']
    GC_content_min = config['GC_content_min']
    GC_content_max = config['GC_content_max']
    Tm_parameters = utils.get_Tm_parameters(config['Tm_parameters'])
    Tm_min = config['Tm_min']
    Tm_max = config['Tm_max']
    
    t = time.time()
    gp.get_probes(number_batchs, batch_size, exon_annotation, genes, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_annotations)
    logging.info('Time to get probes: {} min'.format((time.time() - t)/60))


    #############
    ### filter probes
    #############

    dir_output_blast = os.path.join(dir_output, 'blast/')
    if os.path.isdir(dir_output_blast) == False:
        os.mkdir(dir_output_blast)

    dir_output_probes = os.path.join(dir_output, 'probes/')
    if os.path.isdir(dir_output_probes) == False:
        os.mkdir(dir_output_probes)

    num_threads_blast = config['num_threads']
    word_size = config['word_size'] 
    percent_identity = config['percent_identity']
    coverage = config['coverage'] / 100
    ligation_site = config['ligation_site']
    min_probes_per_gene = config['min_probes_per_gene']

    t = time.time()
    fp.filter_probes_by_GC_Tm_exactmatch(number_batchs, GC_content_min, GC_content_max, Tm_min, Tm_max, dir_output_annotations)
    logging.info('Time to filter with CG vontent, melting temperature and extact matches: {} min'.format((time.time() - t)/60))

    t = time.time()
    fp.run_blast_search(number_batchs, word_size, percent_identity, num_threads_blast, file_gene_gtf, file_genome_fasta, dir_output_annotations, dir_output_blast)
    logging.info('Time to run Blast search: {} min'.format((time.time() - t)/60))

    
    t = time.time()
    fp.filter_probes_by_blast_results(number_batchs, probe_length, coverage, ligation_site, min_probes_per_gene, dir_output_annotations, dir_output_blast, dir_output_probes)
    logging.info('Time to filter with Blast results: {} min'.format((time.time() - t)/60))
    
    
############################################
    
if __name__ == '__main__':
    
    # get comman line arguments
    parameters = args()

    config = utils.get_config(parameters.config)
    utils.print_config(config, logging)

    dir_output = parameters.output


    logging.info('#########Start Pipeline#########')
    probe_pipeline(config, dir_output)
    logging.info('#########End Pipeline#########')
