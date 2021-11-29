############################################
# imports
############################################

import os
import yaml
import gzip
import shutil

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

def get_gene_transcrip_exon_mapping(exon_annotation, dir_output):
    """
    Get a maaping for exons: to which transcript does the exon belong to, to which gene does the transcript belong to. 
    Parameters
    ----------
        exon_annotation: pandas.DataFrame
            Dataframe with exon annotation.
        dir_output: string
            Path to output directory for mapping csv file.
    Returns
    -------
        --- none ---
    """
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

def get_gene_list(annotation):
    """
    Retrive list of unique gene identifiers from annotation table.
    Parameters
    ----------
        annotation: pandas.DataFrame
            Dataframe with genomic annotation.
    Returns
    -------
        genes: list
            List of unique genes.
    """
    genes = sorted(list(annotation['gene_id'].unique()))
    return genes


############################################

def get_transcript_list(annotation):
    """
    Retrive list of unique transcript identifiers from annotation table.
    Parameters
    ----------
        annotation: pandas.DataFrame
            Dataframe with genomic annotation.
    Returns
    -------
        transcripts: list
            List of unique transcripts.
    """
    transcripts = sorted(list(annotation['transcript_id'].unique()))
    return transcripts


############################################

def get_exon_list(annotation):
    """
    Retrive list of unique exon identifiers from annotation table. 
    Parameters
    ----------
        annotation: pandas.DataFrame
            Dataframe with genomic annotation.
    Returns
    -------
        exons: list
            List of unique exons.
    """
    exons = sorted(list(annotation['exon_id'].unique()))
    return exons


############################################