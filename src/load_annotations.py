
############################################
# imports
############################################

import os
import itertools

import pandas as pd

from ftplib import FTP
from Bio import SeqIO

import src.utils as utils


############################################
# functions
############################################


def ftp_download(ftp_link, directory, file_pattern, dir_output):
    """
    ...
    Parameters
    ----------
        ftp_link: 
            
        directory: 
            
        file_pattern: 
            
        dir_output:
    Returns
    -------
        file_output: 
            
    """

    ftp = FTP(ftp_link)
    ftp.login() #login to ftp server
    ftp.cwd(directory) #move to directory

    allfiles = [] 
    ftp.retrlines('LIST ', allfiles.append) 

    for file in allfiles:
        if file_pattern in file: 
            if " -> " in file: # if file is a symbolic link
                file_download = file.split(" -> ")[1]
                file_output = file_download.split('/')[-1]
                ftp.retrbinary('RETR ' + file_download, open(os.path.join(dir_output, file_output), 'wb').write)
            else:
                file_download = file.split(" ")[-1]
                file_output = file_download
                ftp.retrbinary('RETR ' + file_download, open(os.path.join(dir_output, file_output), 'wb').write)

    ftp.quit()

    return file_output


############################################

def download_chr_mapping(dir_output):
    """
    ...
    Parameters
    ----------
        dir_output: 
            
    Returns
    -------
        mapping: 
            
    """

    ftp_link = 'ftp.ncbi.nlm.nih.gov'
    directory = 'genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/'
    file_pattern = 'assembly_report.txt'

    file_output = ftp_download(ftp_link, directory, file_pattern, dir_output)
    file_mapping = os.path.join(dir_output, file_output)

    # skip comment lines but keep last comment line for header
    with open(file_mapping) as input:
        *_comments, names = itertools.takewhile(lambda line: line.startswith('#'), input)
        names = names[1:].split()

    assembly_report = pd.read_table(file_mapping, names=names, sep="\t", comment='#')

    mapping_chromosome = assembly_report[assembly_report['Sequence-Role'] == 'assembled-molecule']
    mapping_chromosome = pd.Series(mapping_chromosome['Sequence-Name'].values, index=mapping_chromosome['RefSeq-Accn']).to_dict()

    mapping_scaffolds = assembly_report[assembly_report['Sequence-Role'] != 'assembled-molecule']
    mapping_scaffolds = pd.Series(mapping_scaffolds['GenBank-Accn'].values, index=mapping_scaffolds['RefSeq-Accn']).to_dict()

    mapping = dict(mapping_chromosome)
    mapping.update(mapping_scaffolds)

    return mapping


############################################

def download_gene_gtf(source, release_ensemble, mapping, dir_output):
    """
    ...
    Parameters
    ----------
        source: 
            
        release_ensemble: 

        mapping:
            
        dir_output: 
            
    Returns
    -------
        file_gene_gtf: 
            
    """

    if source == 'ensemble':
        ftp_link = 'ftp.ensembl.org'
        directory = 'pub/release-{}/gtf/homo_sapiens/'.format(release_ensemble)
        file_pattern = '{}.gtf'.format(release_ensemble)
    elif source == 'ncbi':
        ftp_link = 'ftp.ncbi.nlm.nih.gov'
        directory = '/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/'
        file_pattern = 'genomic.gtf'
    else:
        raise ValueError('Error: unknown source "{}"'.format(source))
    
    file_output = ftp_download(ftp_link, directory, file_pattern, dir_output)
    file_gene_gtf = utils.decompress_gzip(os.path.join(dir_output, file_output))

    if source == 'ncbi':
        process_ncbi_gene_gtf(dir_output, file_gene_gtf, mapping)
        
    return file_gene_gtf


############################################

def process_ncbi_gene_gtf(dir_output, file_gene_gtf, mapping):
    """
    ...
    Parameters
    ----------
        dir_output: 
            
        file_gene_gtf: 

        mapping:
            
    Returns
    -------
        --- None --- 
            
    """
    file_tmp = os.path.join(dir_output, 'temp.gtf')
    output = open(file_tmp, 'w')

    # write comment lines to new file
    with open(file_gene_gtf) as input:
        *_comments, names = itertools.takewhile(lambda line: line.startswith('#'), input)
        output.write(names)

    # read gtf file without comment lines
    gene_annotation = pd.read_table(file_gene_gtf, names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'], sep='\t', comment='#')

    # replace ncbi with genbank chromosome annotation
    for accession_number in gene_annotation.seqname.unique():
        if accession_number in mapping:
            gene_annotation.loc[gene_annotation.seqname == accession_number, 'seqname'] = mapping[accession_number]
        else:
            print('No mapping for accession number: {}'.format(accession_number))
            gene_annotation = gene_annotation[gene_annotation.seqname != accession_number]
    
    gene_annotation.to_csv(output, sep='\t', header=False, index = False)
    os.replace(file_tmp, file_gene_gtf)


############################################

def download_genome_fasta(source, release_ensemble, mapping, dir_output):
    """
    ...
    Parameters
    ----------
        source: 
            
        release_ensemble: 

        mapping:
            
        dir_output: 
            
    Returns
    -------
        file_genome_fasta: 
            
    """
    if source == 'ensemble':
        ftp_link = 'ftp.ensembl.org'
        directory = 'pub/release-{}/fasta/homo_sapiens/dna'.format(release_ensemble)
        file_pattern = 'dna_rm.primary_assembly'
    elif source == 'ncbi':
        ftp_link = 'ftp.ncbi.nlm.nih.gov'
        directory = 'refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers'
        file_pattern = 'genomic.fna'
    else:
        raise ValueError('Error: unknown source "{}"'.format(source))
    
    file_output = ftp_download(ftp_link, directory, file_pattern, dir_output)
    file_genome_fasta = utils.decompress_gzip(os.path.join(dir_output, file_output))

    if source == 'ncbi':
        process_ncbi_genome_fasta(dir_output, file_genome_fasta, mapping)

    return file_genome_fasta


############################################

def process_ncbi_genome_fasta(dir_output, file_genome_fasta, mapping):
    """
    ...
    Parameters
    ----------
        dir_output: 
            
        file_genome_fasta: 

        mapping:
            
    Returns
    -------
        --- None --- 
            
    """   
    file_tmp = os.path.join(dir_output, 'temp.fna')
    output = open(file_tmp, 'w')

    for chromosome_sequnece in SeqIO.parse(file_genome_fasta,'fasta'):
        accession_number = chromosome_sequnece.id
        if accession_number in mapping:
            chromosome_sequnece.id = mapping[accession_number]
            chromosome_sequnece.name = mapping[accession_number]
            chromosome_sequnece.description = chromosome_sequnece.description.replace(accession_number, mapping[accession_number])
            SeqIO.write(chromosome_sequnece, output, 'fasta')
        else:
            print('No mapping for accession number: {}'.format(accession_number))
    
    os.replace(file_tmp, file_genome_fasta)


############################################




