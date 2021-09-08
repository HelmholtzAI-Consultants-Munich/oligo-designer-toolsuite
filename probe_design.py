############################################
# imports
############################################

import time
import argparse
import logging
import pickle
import copy

import src.utils as utils
import src.load_exome as le
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
    get_annotation = config['get_annotation']

    #############
    ### get exon sequences
    #############
    if get_annotation:
        source = config['source']
        release_ensemble = config['release_ensemble']
        
        if source == 'ncbi':
            # if source ncbi we need a chromosome mapping
            mapping = le.download_chr_mapping(dir_output)
        else:
            mapping = None
        
        if source == 'ensemble':
            logging.info('Using Ensemble release {}'.format(release_ensemble))

        logging.info('Download gene annotation from: {}'.format(source))
        file_gene_annotation = le.download_gene_gtf(dir_output, source, release_ensemble)
        logging.info('File gene annotation: {}'.format(file_gene_annotation))

        logging.info('Get exome annotation')
        file_exon_annotation = le.get_exons_annotation(file_gene_annotation, source, mapping)
        logging.info('File exon annotation: {}'.format(file_exon_annotation))

        logging.info('Download genome annotation from: {}'.format(source))
        file_genome_sequence = le.download_genome_fasta(dir_output, source, release_ensemble, mapping)
        logging.info('File genome sequence: {}'.format(file_genome_sequence))

        logging.info('Get exome fasta')
        file_exon_sequence = le.get_exons_fasta(file_exon_annotation, file_genome_sequence)
        logging.info('File exon sequence: {}'.format(file_exon_sequence))

    else:
        file_gene_annotation = config['file_gene_annotation']
        file_exon_annotation = config['file_exon_annotation']
        file_genome_sequence = config['file_genome_sequence']
        file_exon_sequence = config['file_exon_sequence']


    #############
    ### get probes
    #############
    probe_length = config['probe_length']
    min_probes_per_gene = config['min_probes_per_gene']
    GC_content_min = config['GC_content_min']
    GC_content_max = config['GC_content_max']
    Tm_min = config['Tm_min']
    Tm_max = config['Tm_max']

    t = time.time()
    mapping_exon_to_gene = gp.get_mapping_gene_exon(dir_output, file_gene_annotation, get_annotation)
    logging.info('Time to process gene-exon mapping dicts: {} s'.format(time.time() - t))


    t = time.time()
    mapping_gene_to_probes = gp.get_mapping_gene_probes(dir_output, file_exon_sequence, mapping_exon_to_gene, probe_length, GC_content_min, GC_content_max, Tm_min, Tm_max, get_annotation)
    logging.info('Time to process gene-probe mapping dict: {} s'.format(time.time() - t))


    t = time.time()
    mapping_gene_to_probes_unique = gp.get_unique_probes(dir_output, mapping_gene_to_probes, copy.deepcopy(mapping_gene_to_probes), get_annotation)
    logging.info('Time to process unique probes: {} s'.format(time.time() - t))

    
    removed_genes = []
    for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_unique.items():
        if len(mapping_probe_to_exon) < min_probes_per_gene:
            removed_genes.append(gene_id)
    logging.info('{} genes were filtered out because they had < {} probes per gene.'.format(len(removed_genes), min_probes_per_gene))


    #############
    ### filter probes
    #############
    word_size = config['word_size']
    coverage = config['coverage']
    percent_identity = config['percent_identity']
    num_threads = int(config['num_threads'])

    t = time.time()
    mapping_gene_to_probes_blastn = fp.filter_probes_with_blastn(dir_output, file_exon_sequence, mapping_gene_to_probes_unique, mapping_exon_to_gene, min_probes_per_gene, probe_length, word_size, coverage, percent_identity, num_threads)
    logging.info('Time to process blastn filter: {} s'.format(time.time() - t))

    #############
    ### save results
    #############
    utils.write_output(dir_output, mapping_gene_to_probes_blastn)
    

    
############################################
    
if __name__ == '__main__':
    
    # get comman line arguments
    parameters = args()
    config = utils.get_config(parameters.config) 

    logging.basicConfig(filename='{}logging.txt'.format(config['dir_output']), level=logging.NOTSET)
    utils.print_config(config, logging)


    logging.info('#########Start Pipeline#########')

    t = time.time()
    probe_pipeline(config)
    logging.info('Time for whole pipeline: {} s'.format(time.time() - t))

    logging.info('#########End Pipeline#########')
