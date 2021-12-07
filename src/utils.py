############################################
# imports
############################################

import os
import yaml
import gzip
import shutil

import pybedtools
import gtfparse

from Bio.SeqUtils import MeltingTemp as mt

############################################
# helper functions
############################################

def get_config(config):
    """
    Loads config file in yaml format. 
    Parameters
    ----------
        config: string
            Path to yaml config file.
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
    Prints formatted config parameters as <parameter_name>: <parameter_value> to logging file. 
    Parameters
    ----------
        config: string
            Path to yaml config file.
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
    Decompress zip files.
    Parameters
    ----------
        file_gzip: string
            Path to zipped file.
            
    Returns
    -------
        file_output: string
            Path to unzipped file.
            
    """
    file_output = file_gzip.split('.gz')[0]
    with gzip.open(file_gzip, 'rb') as f_in:
        with open(file_output, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return file_output


############################################

def load_exon_annotation(file_gene_gtf):
    """
    Load exon annotation from gtf file. 
    Parameters
    ----------
        file_gene_gtf: string
            Path to gtf file with gene annotation.
    Returns
    -------
        exon_annotation: pandas.DataFrame
            Dataframe with exon annotation.
            
    """
    gene_annotation = gtfparse.read_gtf(file_gene_gtf)
    exon_annotation = gene_annotation.loc[gene_annotation['feature'] == 'exon']
    exon_annotation = exon_annotation.assign(source='unknown')
    if not 'exon_id' in exon_annotation.columns:
        exon_annotation['exon_id'] = exon_annotation['transcript_id'] + '_exon' + exon_annotation['exon_number']

    return exon_annotation


############################################

def load_transcriptome_annotation(file_gene_gtf):
    """
    Load transcriptome annotation from gtf file. 
    Parameters
    ----------
        file_gene_gtf: string
            Path to gtf file with gene annotation.
    Returns
    -------
        transcriptome_annotation: pandas.DataFrame
            Dataframe with transcript annotation.
    """
    gene_annotation = gtfparse.read_gtf(file_gene_gtf)
    transcriptome_annotation = gene_annotation.loc[gene_annotation['feature'] == 'transcript']
    transcriptome_annotation = transcriptome_annotation.assign(source='unknown')

    return transcriptome_annotation


############################################

def get_fasta(file_gtf, file_genome_fasta, file_fasta):
    """
    Get sequence for regions annotated in input gft file using genome fasta file. 
    Parameters
    ----------
        file_gtf: string
            Path to gtf file with annotated genomic regions.
        file_genome_fasta: string
            Path to fasta file with genome sequence.
        file_fasta: string
            Path to fasta file where retrieved sequences are written to.
    Returns
    -------
        --- none ---
    """
    annotation = pybedtools.BedTool(file_gtf)
    genome_sequence = pybedtools.BedTool(file_genome_fasta)

    annotation = annotation.sequence(fi=genome_sequence, s=True, name=True)
    annotation.save_seqs(file_fasta)


############################################

def read_gene_list(file_genes):
    """
    Read list of genes from text file. Gene names must be provided in ensemble or NCBI annotation format. 
    Parameters
    ----------
        file_genes: string
            Path to text file with gene names.
    Returns
    -------
        genes: list
            List of gene names.
    """
    with open(file_genes) as file:
        lines = file.readlines()
        genes = [line.rstrip() for line in lines]
    return genes


############################################

def get_Tm_parameters(Tm_parameters):
    """
    Convert config parameters for 'table' attributes of MeltingTemp function into MeltingTemp attributes. 
    Parameters
    ----------
        Tm_parameters: dict
            Dictionary with parameters for MeltingTemp function.
    Returns
    -------
        Tm_parameters: dict
            Dictionary with parameters for MeltingTemp function.
    """
    # get the attributed for the parameter
    Tm_parameters['nn_table'] = getattr(mt, Tm_parameters['nn_table'])
    Tm_parameters['tmm_table'] = getattr(mt, Tm_parameters['tmm_table'])
    Tm_parameters['imm_table'] = getattr(mt, Tm_parameters['imm_table'])
    Tm_parameters['de_table'] = getattr(mt, Tm_parameters['de_table'])

    return Tm_parameters