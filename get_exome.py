
############################################
# imports
############################################

import argparse
import os
import pybedtools
from gtfparse import read_gtf
from ftplib import FTP

############################################
# pipeline
############################################

def args():
    """
    Returns command line arguments. 
    Parameters
    ----------
        --- none ---

    Returns
    -------
        args_parser: namespace
            Namespace object with argument attributes.
    """
    
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('-d','--dirdata',help='directory for data storage',type=str,required=True)
    args_parser.add_argument('-r','--release',help='ensemble release',type=str,required=True)
    return args_parser.parse_args()



def download_genome_fasta(dir_data, release='104'):
    ftp = FTP('ftp.ensembl.org')
    ftp.login() 

    # move to genome directory
    directory = 'pub/release-{}/fasta/homo_sapiens/dna'.format(release)
    ftp.cwd(directory)

    allfiles = ftp.nlst()

    for file in allfiles:
        if 'dna_rm.primary_assembly' in file:
            file_download = file
            ftp.retrbinary('RETR ' + file_download, open(os.path.join(dir_data, file_download), 'wb').write)

    ftp.quit()

    file_genome_sequence = '{}{}'.format(dir_data, file_download)
    os.system("gunzip {}".format(file_genome_sequence))
    file_genome_sequence = file_genome_sequence.split('.gz')[0]

    return file_genome_sequence



def download_gene_gtf(dir_data, release='104'):
    ftp = FTP('ftp.ensembl.org')
    ftp.login() 

    # move to annotation directory
    directory = 'pub/release-{}/gtf/homo_sapiens/'.format(release)
    ftp.cwd(directory)

    allfiles = ftp.nlst()

    for file in allfiles:
        if '{}.gtf.gz'.format(release) in file:
            file_download = file
            ftp.retrbinary('RETR ' + file_download, open(os.path.join(dir_data, file_download), 'wb').write)

    ftp.quit()

    file_gene_annotation = '{}{}'.format(dir_data, file_download)
    os.system("gunzip {}".format(file_gene_annotation))
    file_gene_annotation = file_gene_annotation.split('.gz')[0]

    return file_gene_annotation


def get_exons_annotation(file_gene_annotation):

    file_exon_annotation = '{}.exons.gtf'.format(file_gene_annotation.split('.gtf')[0])
    
    # extract exons from gene annotation file
    gene_annotation = read_gtf(file_gene_annotation)
    
    exon_annotation = gene_annotation[gene_annotation["feature"] == "exon"]
    #exon_annotation = exon_annotation[exon_annotation["transcript_biotype"] == "protein_coding"]
    exon_annotation[['seqname','source','exon_id','start','end','score','strand','frame','gene_id']].to_csv(file_exon_annotation, sep='\t', header=False, index = False)
    
    return file_exon_annotation
        
    
def get_exons_fasta(file_exon_annotation, file_genome_sequence):
    
    file_exon_sequence = '{}.fa'.format(file_exon_annotation.split('.gtf')[0])
    
    # get sequence for exons
    exon_annotation = pybedtools.BedTool(file_exon_annotation)
    genome_sequence = pybedtools.BedTool(file_genome_sequence)

    exon_annotation = exon_annotation.sequence(fi=genome_sequence, s=True, name=True)
    exon_annotation.save_seqs(file_exon_sequence)
    
    return file_exon_sequence



def main():
    
    parameters = args()

    dir_data = parameters.dirdata
    release = parameters.release

    file_genome_sequence = download_genome_fasta(dir_data, release)
    print(file_genome_sequence)

    file_gene_annotation = download_genome_fasta(dir_data, release)
    print(file_gene_annotation)

    print('get exome annotation')
    file_exon_annotation = get_exons_annotation(file_gene_annotation)
    print('get exome fasta')
    file_exon_sequence = get_exons_fasta(file_exon_annotation, file_genome_sequence)



if __name__ == '__main__':
    
    main()