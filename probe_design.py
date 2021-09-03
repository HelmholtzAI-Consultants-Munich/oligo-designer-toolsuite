############################################
# imports
############################################

import time
import argparse
import logging
import pickle

import src.utils as utils
import src.load_exome as lp_exome
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

def main():
    
    
    # get comman line arguments
    parameters = args()
    config = utils.get_config(parameters.config) 
    dir_output = config['dir_output']
    
    logging.basicConfig(filename='{}logging.log'.format(dir_output), level=logging.INFO)
    utils.print_config(config, logging)
    logging.info('#########Start Pipeline#########')


    ### get exon sequences
    if config['get_gene_and_exon_annotation_file']:
        source = config['source']
        release_ensemble = config['release_ensemble']
        
        if source == 'ncbi':
            # if source ncbi we need a chromosome mapping
            mapping = lp_exome.download_chr_mapping(dir_output)
        else:
            mapping = None
        if source == 'ensemble':
            logging.info('Using Ensemble release {}'.format(release_ensemble))

        logging.info('\nDownload gene annotation from: {}'.format(source))
        file_gene_annotation = lp_exome.download_gene_gtf(dir_output, source, release_ensemble)
        logging.info('\nFile gene annotation: {}'.format(file_gene_annotation))

        logging.info('\nGet exome annotation')
        file_exon_annotation = lp_exome.get_exons_annotation(file_gene_annotation, source, mapping)
        logging.info('File exon annotation: {}'.format(file_exon_annotation))
    
    else:
        file_gene_annotation = config['file_gene_annotation']
        file_exon_annotation = config['file_exon_annotation']


    if config['get_genome_annotation_file']:
        source = config['source']
        release_ensemble = config['release_ensemble']

        logging.info('\nDownload genome annotation from: {}'.format(source))
        file_genome_sequence = lp_exome.download_genome_fasta(dir_output, source, release_ensemble, mapping)
        logging.info('File genome sequence: {}'.format(file_genome_sequence))

    else:
        file_genome_sequence = config['file_genome_sequence']

    logging.info('\nGet exome fasta')
    file_exon_sequence = lp_exome.get_exons_fasta(file_exon_annotation, file_genome_sequence)
    logging.info('File exon sequence: {}'.format(file_exon_sequence))



    ### probe design
    probe_length = config['probe_length']
    min_probes_per_gene = config['min_probes_per_gene']
    GC_content_min = config['GC_content_min']
    GC_content_max = config['GC_content_max']
    Tm_min = config['Tm_min']
    Tm_max = config['Tm_max']

    t = time.time()
    file_mapping_exon_to_gene = gp.get_mapping_gene_exon(dir_output, file_gene_annotation)
    logging.info('Time to process gene-exon mapping dicts: {} s'.format(time.time() - t))

    mapping_exon_to_gene = pickle.load(open(file_mapping_exon_to_gene,'rb'))

    t = time.time()
    file_mapping_gene_to_probes = gp.get_mapping_gene_probes(dir_output, file_exon_sequence, mapping_exon_to_gene, probe_length, GC_content_min, GC_content_max, Tm_min, Tm_max)
    logging.info('Time to process gene-probe mapping dict: {} s'.format(time.time() - t))

    mapping_gene_to_probes = pickle.load(open(file_mapping_gene_to_probes,'rb'))
    mapping_gene_to_probes_unique = pickle.load(open(file_mapping_gene_to_probes,'rb'))

    t = time.time()
    file_mapping_gene_to_probes_unique = gp.get_unique_probes(dir_output, mapping_gene_to_probes, mapping_gene_to_probes_unique)
    logging.info('Time to process unique probes: {} s'.format(time.time() - t))

    mapping_gene_to_probes_unique = pickle.load(open(file_mapping_gene_to_probes_unique,'rb'))
    removed_genes = []
    for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_unique.items():
        if len(mapping_probe_to_exon) < min_probes_per_gene:
            removed_genes.append(gene_id)
    logging.info('{} genes were filtered out because they had < {} probes per gene.'.format(len(removed_genes), min_probes_per_gene))



    ### probe filtering
    coverage = config['coverage']
    percent_identity = config['percent_identity']
    num_threads = int(config['num_threads'])

    t = time.time()
    file_mapping_gene_to_probes_blastn = fp.filter_probes_with_blastn_2(dir_output, file_exon_sequence, mapping_gene_to_probes_unique, mapping_exon_to_gene, min_probes_per_gene, probe_length, coverage, percent_identity, num_threads)
    logging.info('Time to process blastn filter: {} s'.format(time.time() - t))
    
    #t = time.time()
    #file_mapping_gene_to_probes_hamming = fp.filter_probes_with_hamming_distance(dir_output, mapping_gene_to_probes_unique, min_probes_per_gene, probe_length, percent_identity)
    #logging.info('Time to process hamming distance filter: {} s'.format(time.time() - t))


    logging.info('#########End Pipeline#########')

    
############################################
    
if __name__ == '__main__':
    
    main()