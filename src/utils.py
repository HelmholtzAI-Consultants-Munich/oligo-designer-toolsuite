############################################
# imports
############################################

import os
import yaml
import gzip
import shutil
import time

import pybedtools
import gtfparse


############################################
# helper functions
############################################

def get_config(config):
    """
    Loads config file in yaml format. 
    Parameters
    ----------
        config: string
            File path to yaml config file.
    Returns
    -------
        yaml: dict
            User-defined parameters, where keys are the parameter names and values are the paremeter values.
    """
    with open(config, 'r') as ymlfile:
        return yaml.safe_load(ymlfile)


############################################
    
def print_config(config, logger):
    """
    Logs formatted config parameters as <parameter_name>: <parameter_value>. 
    Parameters
    ----------
        config: string
            File path to yaml config file.
        logging: logger
            Logger object to store important information.
    Returns
    -------
        --- none ---
    """
    logger.info('#########Parameter settings#########')
    for item, value in config.items(): 
        logger.info("{}: {}".format(item, value))

    
############################################
    
def decompress_gzip(file_gzip):
    """
    ... 
    Parameters
    ----------
        file_gzip: string
            
    Returns
    -------
        file_output: string
            
    """

    file_output = file_gzip.split('.gz')[0]
    with gzip.open(file_gzip, 'rb') as f_in:
        with open(file_output, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return file_output


############################################

def processes_wait(jobs):

    while True:
        all_finished = True
        
        for i, job in enumerate(jobs):
            if job.ready():
                print("Process %s has finished" % i)
            else:
                all_finished = False

        if all_finished:
            break

        time.sleep(1)


############################################

def load_exon_annotation(file_gene_gtf):
    gene_annotation = gtfparse.read_gtf(file_gene_gtf)
    exon_annotation = gene_annotation.loc[gene_annotation['feature'] == 'exon']
    exon_annotation = exon_annotation.assign(source='unknown')
    if not 'exon_id' in exon_annotation.columns:
        exon_annotation['exon_id'] = exon_annotation['transcript_id'] + '_exon' + exon_annotation['exon_number']

    return exon_annotation


############################################

def load_transcriptome_annotation(file_gene_gtf):
    gene_annotation = gtfparse.read_gtf(file_gene_gtf)
    transcriptome_annotation = gene_annotation.loc[gene_annotation['feature'] == 'transcript']
    transcriptome_annotation = transcriptome_annotation.assign(source='unknown')

    return transcriptome_annotation


############################################

def get_fasta(file_gtf, file_genome_fasta, file_fasta):
 
    annotation = pybedtools.BedTool(file_gtf)
    genome_sequence = pybedtools.BedTool(file_genome_fasta)

    annotation = annotation.sequence(fi=genome_sequence, s=True, name=True)
    annotation.save_seqs(file_fasta)


############################################

def get_gene_list(annotation):
    genes = sorted(list(annotation['gene_id'].unique()))
    return genes


############################################

def get_transcript_list(annotation):
    transcripts = sorted(list(annotation['transcript_id'].unique()))
    return transcripts


############################################

def get_exon_list(annotation):
    exons = sorted(list(annotation['exon_id'].unique()))
    return exons


############################################

def get_gene_transcrip_exon_mapping(exon_annotation, dir_output):
    file_mapping = os.path.join(dir_output, 'mapping_gene_transcript_exon.txt')

    with open(file_mapping, 'w') as handle:
        handle.write('gene_id\ttranscript_id\texon_id\n')

        for gene in get_gene_list(exon_annotation):
            exon_annotation_gene = exon_annotation[exon_annotation['gene_id'] == gene]
            for transcript in get_transcript_list(exon_annotation_gene):
                exon_annotation_transcript = exon_annotation_gene[exon_annotation_gene['transcript_id'] == transcript]
                for exon in get_exon_list(exon_annotation_transcript):
                    handle.write('{}\t{}\t{}\n'.format(gene, transcript, exon))


############################################