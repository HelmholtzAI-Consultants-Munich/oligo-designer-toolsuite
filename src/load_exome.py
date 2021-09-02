
############################################
# imports
############################################

import os
import itertools
import pybedtools
from gtfparse import read_gtf
from Bio import SeqIO
from ftplib import FTP
import pandas as pd
import subprocess as sp


############################################
# pipeline
############################################

def ftp_download(dir_data, ftp_link, directory, file_pattern):
    # login and move to directory
    ftp = FTP(ftp_link)
    ftp.login() 
    ftp.cwd(directory)

    allfiles = [] 
    ftp.retrlines('LIST ', allfiles.append) 

    for file in allfiles:
        if file_pattern in file: 
            if " -> " in file: # if file is a symbolic link
                file_download = file.split(" -> ")[1]
                file_output = file_download.split('/')[-1]
                ftp.retrbinary('RETR ' + file_download, open(os.path.join(dir_data, file_output), 'wb').write)
            else:
                file_download = file.split(" ")[-1]
                file_output = file_download
                ftp.retrbinary('RETR ' + file_download, open(os.path.join(dir_data, file_output), 'wb').write)

    ftp.quit()

    return file_output


############################################

def download_gene_gtf(dir_data, source, release_ensemble):
    # ensemble: release_ensemble = 104
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
    
    file_output = ftp_download(dir_data, ftp_link, directory, file_pattern)

    file_gene_annotation = '{}{}'.format(dir_data, file_output)
    cmd = "gunzip {}".format(file_gene_annotation)
    sp.run(cmd, shell=True)
    file_gene_annotation = file_gene_annotation.split('.gz')[0]
        
    return file_gene_annotation


############################################

def download_genome_fasta(dir_data, source, release_ensemble, file_chr_mapping):
    # ensemble: release = 104
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
    
    file_output = ftp_download(dir_data, ftp_link, directory, file_pattern)

    file_genome_sequence = '{}{}'.format(dir_data, file_output)
    cmd = "gunzip {}".format(file_genome_sequence)
    sp.run(cmd, shell=True)
    file_genome_sequence = file_genome_sequence.split('.gz')[0]

    if source == 'ncbi':
        process_ncbi_genome(dir_data, file_genome_sequence, file_chr_mapping)

    return file_genome_sequence


############################################

def download_chr_mapping(dir_data):
    ftp_link = 'ftp.ncbi.nlm.nih.gov'
    directory = 'genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/'
    file_pattern = 'assembly_report.txt'

    file_output = ftp_download(dir_data, ftp_link, directory, file_pattern)

    file_mapping = '{}{}'.format(dir_data, file_output)

    with open(file_mapping) as output:
        *_comments, names = itertools.takewhile(lambda line: line.startswith('#'), output)
        names = names[1:].split()
    output.close()

    assembly_report = pd.read_table(file_mapping, names=names, sep="\t", comment='#')

    mapping_chromosome = assembly_report[assembly_report['Sequence-Role'] == 'assembled-molecule']
    mapping_chromosome = pd.Series(mapping_chromosome['Sequence-Name'].values, index=mapping_chromosome['RefSeq-Accn']).to_dict()

    mapping_scaffolds = assembly_report[assembly_report['Sequence-Role'] != 'assembled-molecule']
    mapping_scaffolds = pd.Series(mapping_scaffolds['GenBank-Accn'].values, index=mapping_scaffolds['RefSeq-Accn']).to_dict()

    mapping = dict(mapping_chromosome)
    mapping.update(mapping_scaffolds)

    return mapping


############################################

def process_ncbi_genome(dir_data, file_genome_sequence, mapping):
    
    file_output = '{}temp.fna'.format(dir_data)
    output = open(file_output, 'w')

    for chromosome_sequnece in SeqIO.parse(file_genome_sequence,'fasta'):

        accession_number = chromosome_sequnece.id
        if accession_number in mapping:
            chromosome_sequnece.id = mapping[accession_number]
            chromosome_sequnece.name = mapping[accession_number]
            chromosome_sequnece.description = chromosome_sequnece.description.replace(accession_number, mapping[accession_number])
            SeqIO.write(chromosome_sequnece, output, 'fasta')
        else:
            print('No mapping for accession number: {}'.format(accession_number))
    
    cmd = 'mv {} {}'.format(file_output, file_genome_sequence)
    sp.run(cmd, shell=True)


############################################

def process_ncbi_annotation(exon_annotation, mapping):

    exon_annotation['exon_id'] = exon_annotation['transcript_id'] + '_exon' + exon_annotation['exon_number']
    exon_annotation['source'] = 'Refseq'

    accession_numbers = exon_annotation.seqname.unique()
    for accession_number in accession_numbers:
        if accession_number in mapping:
            exon_annotation.seqname[exon_annotation.seqname == accession_number] = mapping[accession_number]
        else:
            print('No mapping for accession number: {}'.format(accession_number))
            exon_annotation = exon_annotation[exon_annotation.seqname != accession_number]
    return exon_annotation


############################################

def get_exons_annotation(file_gene_annotation, source, file_chr_mapping):

    file_exon_annotation = '{}.exons.gtf'.format(file_gene_annotation.split('.gtf')[0])
    
    # extract exons from gene annotation file
    gene_annotation = read_gtf(file_gene_annotation)
    exon_annotation = gene_annotation[gene_annotation["feature"] == "exon"]
    if source == 'ncbi':
        exon_annotation = process_ncbi_annotation(exon_annotation, file_chr_mapping)

    exon_annotation[['seqname','source','exon_id','start','end','score','strand','frame','gene_id']].to_csv(file_exon_annotation, sep='\t', header=False, index = False)
    
    return file_exon_annotation


############################################

def get_exons_fasta(file_exon_annotation, file_genome_sequence):
    
    file_exon_sequence = '{}.fna'.format(file_exon_annotation.split('.gtf')[0])
    
    # get sequence for exons
    exon_annotation = pybedtools.BedTool(file_exon_annotation)
    genome_sequence = pybedtools.BedTool(file_genome_sequence)

    exon_annotation = exon_annotation.sequence(fi=genome_sequence, s=True, name=True)
    exon_annotation.save_seqs(file_exon_sequence)
    
    cmd = 'rm {}.fai'.format(file_genome_sequence)
    sp.run(cmd, shell=True)

    return file_exon_sequence


############################################

