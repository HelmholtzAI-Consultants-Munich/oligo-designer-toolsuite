
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


def ftp_download(ftp_link, directory, file_name, dir_output):
    """
    Download file from ftp server.
    Parameters
    ----------
        ftp_link: string
            Link to ftp server, e.g. 'ftp.ncbi.nlm.nih.gov'.
        directory: string
            Directory on ftp server, where the file is located, e.g. 'refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers'.
        file_name: string
            Name of file that should be downloaded from ftp server, e.g. 'GCF_000001405.39_GRCh38.p13_genomic.fna'.
        dir_output: string
            Path to directory for downloaded files.
    Returns
    -------
        file_output: string
            Path to downloaded file.
    """
    ftp = FTP(ftp_link)
    ftp.login() #login to ftp server
    ftp.cwd(directory) #move to directory

    allfiles = [] 
    ftp.retrlines('LIST ', allfiles.append) 

    for file in allfiles:
        if file_name in file: 
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

def download_chr_mapping(ftp, dir_output):
    """
    Download file with mapping of chromosome names between GenBank and Ref-Seq accession number from ftp server and create a mapping dictionary.
    Parameters
    ----------
        ftp: dict
            Dictionary with ftp parameters (ftp_link: link to ftp server; directory: directory on ftp server; file_name: name of file).
        dir_output: string
            Path to directory for downloaded files.
    Returns
    -------
        mapping: dict
            Dictionary with mapping of chromsome names from GenBank to Ref-Seq.
    """
    file_output = ftp_download(ftp['ftp_link'], ftp['directory'], ftp['file_name'], dir_output)
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

def download_gene_gtf(source, mapping, ftp, dir_output):
    """
    Download gene annotation in gtf file format from ftp server and unzip file. If gene annotation comes from ncbi, map chromosome annotation to Ref-Seq accession number.
    Parameters
    ----------
        source: string
            Source for ftp download, e.g. 'ncbi' or 'ensemble'.
        mapping: dict
            Chromosome mapping dictionary (GenBank to Ref-Seq).
        ftp: dict
            Dictionary with ftp parameters (ftp_link: link to ftp server; directory: directory on ftp server; file_name: name of file).
        dir_output: string
            Path to directory for downloaded files.
    Returns
    -------
        file_gene_gtf: string
            Path to downloaded file. 
    """
    file_output = ftp_download(ftp['ftp_link'], ftp['directory'], ftp['file_name'], dir_output)
    file_gene_gtf = utils.decompress_gzip(os.path.join(dir_output, file_output))

    if source == 'ncbi':
        process_ncbi_gene_gtf(dir_output, file_gene_gtf, mapping)
        
    return file_gene_gtf


############################################

def process_ncbi_gene_gtf(dir_output, file_gene_gtf, mapping):
    """
    Process gene annotation file downloaded from NCBI: map chromosome annotation to Ref-Seq.
    Parameters
    ----------
        dir_output: string
            Path to directory where the processed file should be saved.
        file_gene_gtf: string
            Path to gtf file with gene annotation.
        mapping: dict
            Chromosome mapping dictionary (GenBank to Ref-Seq).
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

def download_genome_fasta(source, mapping, ftp, dir_output):
    """
    Download genome sequence in fasta file format from ftp server and unzip file. If genome sequence comes from ncbi, map chromosome annotation to Ref-Seq accession number.
    Parameters
    ----------
        source: string
            Source for ftp download, e.g. 'ncbi' or 'ensemble'.
        mapping: dict
            Chromosome mapping dictionary (GenBank to Ref-Seq).
        ftp: dict
            Dictionary with ftp parameters (ftp_link: link to ftp server; directory: directory on ftp server; file_name: name of file).
        dir_output: string
            Path to directory for downloaded files.
    Returns
    -------
        file_genome_fasta: string
            Path to downloaded file. 
    """
    file_output = ftp_download(ftp['ftp_link'], ftp['directory'], ftp['file_name'], dir_output)
    file_genome_fasta = utils.decompress_gzip(os.path.join(dir_output, file_output))

    if source == 'ncbi':
        process_ncbi_genome_fasta(dir_output, file_genome_fasta, mapping)

    return file_genome_fasta


############################################

def process_ncbi_genome_fasta(dir_output, file_genome_fasta, mapping):
    """
    Process genome sequence file downloaded from NCBI: map chromosome annotation to Ref-Seq.
    Parameters
    ----------
        dir_output: string
            Path to directory where the processed file should be saved.
        file_genome_fasta: string
            Path to fasta file with genome sequence.
        mapping: dict
            Chromosome mapping dictionary (GenBank to Ref-Seq).
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




