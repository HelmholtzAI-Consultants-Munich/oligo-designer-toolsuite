############################################
# imports
############################################

import os
import re
import yaml
import gzip
import shutil

import pybedtools

from ftplib import FTP
from Bio.SeqUtils import MeltingTemp as mt

import logging
logging.getLogger('padlock_probe_designer')

############################################
# helper functions
############################################

def get_config(config):
    '''Load config file in yaml format.

    :param config: Path to yaml config file.
    :type config: string
    :return: User-defined parameters, where keys are the parameter names and values are the paremeter values.
    :rtype: dict
    '''
    with open(config, 'r') as ymlfile:
        return yaml.safe_load(ymlfile)

    
def print_config(config):
    '''Prints formatted config parameters as <parameter_name>: <parameter_value> to logging file.

    :param config: User-defined parameters, where keys are the parameter names and values are the paremeter values.
    :type config: dict
    '''
    logging.info('#########Parameter settings#########')
    for item, value in config.items(): 
        logging.info("{}: {}".format(item, value))


def get_Tm_parameters(Tm_parameters,sequence='probe'):
    '''Convert config parameters for 'table' attributes of MeltingTemp function into MeltingTemp attributes. 
    
    And load specific parameters of the given `sequence` type.
    :param Tm_parameters: Dictionary with parameters for MeltingTemp function.
    :type Tm_parameters: dict
    :param sequence: Name of the sequence for which parameters are loaded
    :type sequence: string
    :return: Dictionary with parameters for MeltingTemp function.
    :rtype: dict
    '''
    Tm_params = Tm_parameters['shared'].copy()
    if Tm_parameters[sequence]:
        Tm_params.update(Tm_parameters[sequence])
    
    Tm_params['nn_table'] = getattr(mt, Tm_params['nn_table'])
    Tm_params['tmm_table'] = getattr(mt, Tm_params['tmm_table'])
    Tm_params['imm_table'] = getattr(mt, Tm_params['imm_table'])
    Tm_params['de_table'] = getattr(mt, Tm_params['de_table'])

    return Tm_params

def get_Tm_correction_parameters(Tm_correction_parameters,sequence='probe'):
    '''Load specific Tm correction parameters of the given `sequence` type
    :param Tm_correction_parameters: Dictionary with parameters for MeltingTemp function.
    :type Tm_correction_parameters: dict
    :param sequence: Name of the sequence for which parameters are loaded
    :type sequence: string
    :return: Dictionary with parameters for MeltingTemp function.
    :rtype: dict
    '''
    Tm_params = Tm_correction_parameters['shared'].copy()
    if Tm_correction_parameters[sequence]:
        Tm_params.update(Tm_correction_parameters[sequence])

    return Tm_params


def ftp_download(ftp_link, directory, file_name, dir_output):
    '''Download file from ftp server.

    :param ftp_link: Link to ftp server, e.g. 'ftp.ncbi.nlm.nih.gov'.
    :type ftp_link: string
    :param directory: Directory on ftp server, where the file is located, e.g. 'refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers'.
    :type directory: string
    :param file_name: Name of file that should be downloaded from ftp server, e.g. 'GCF_000001405.39_GRCh38.p13_genomic.fna'.
    :type file_name: string
    :param dir_output: Path to directory for downloaded files.
    :type dir_output: string
    :return: Path to downloaded file.
    :rtype: string
    '''
    
    ftp = FTP(ftp_link)
    ftp.login()  # login to ftp server
    ftp.cwd(directory)  # move to directory

    files = ftp.nlst()

    for file in files:
        if re.match(file_name, file):
            file_output = os.path.join(dir_output, file)
            ftp.retrbinary('RETR ' + file, open(file_output, 'wb').write)

    ftp.quit()

    return file_output


def decompress_gzip(file_gzip):
    '''Decompress zip files.

    :param file_gzip: Path to zipped file.
    :type file_gzip: string
    :return: Path to unzipped file.
    :rtype: string
    '''
    file_output = file_gzip.split('.gz')[0]
    with gzip.open(file_gzip, 'rb') as f_in:
        with open(file_output, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_gzip)
    return file_output


def get_ncbi_ftp_parameters(species, annotation_release, genome_assembly, dir_output):
    '''_summary_

    :param species: _description_
    :type species: _type_
    :param annotation_release: _description_
    :type annotation_release: _type_
    :param genome_assembly: _description_
    :type genome_assembly: _type_
    :param dir_output: _description_
    :type dir_output: _type_
    :return: _description_
    :rtype: _type_
    '''
    ftp_link = 'ftp.ncbi.nlm.nih.gov'
    if species == 'human':
        ftp_directory = 'refseq/H_sapiens/annotation/annotation_releases/'
    if species == 'mouse':
        ftp_directory = 'refseq/M_musculus/annotation_releases/'

    if annotation_release == 'current':
        ftp_directory = ftp_directory + 'current/'
    else:
        ftp_directory = ftp_directory + f'{annotation_release}/'
    
    file_readme = ftp_download(ftp_link, ftp_directory, 'README', dir_output)
    with open(file_readme, 'r') as handle:
        for line in handle:
            if line.startswith('ASSEMBLY NAME:'):
                assembly_name = line.strip().split('\t')[1]
            if line.startswith('ASSEMBLY ACCESSION:'):
                assembly_accession = line.strip().split('\t')[1]
                break
    os.remove(file_readme)
    ftp_directory = ftp_directory + f'{assembly_accession}_{assembly_name}'
    ftp_file_gtf = f'{assembly_accession}_{assembly_name}_genomic.gtf.gz'
    ftp_file_fasta = f'{assembly_accession}_{assembly_name}_genomic.fna.gz'
    ftp_file_chr_mapping = f'{assembly_accession}_{assembly_name}_assembly_report.txt'

    return ftp_link, ftp_directory, ftp_file_gtf, ftp_file_fasta, ftp_file_chr_mapping


def get_ensemble_ftp_parameters(species, annotation_release, genome_assembly, dir_output):
    '''_summary_

    :param species: _description_
    :type species: _type_
    :param annotation_release: _description_
    :type annotation_release: _type_
    :param genome_assembly: _description_
    :type genome_assembly: _type_
    :param dir_output: _description_
    :type dir_output: _type_
    :return: _description_
    :rtype: _type_
    '''
    ftp_link = 'ftp.ensembl.org'

    if species == 'human':
        species_id = 'homo_sapiens'
    if species == 'mouse':
        species_id = 'mus_musculus'

    if annotation_release == 'current':
        file_readme = ftp_download(ftp_link, 'pub/', 'current_README', dir_output)
        with open(file_readme, 'r') as handle:
            for line in handle:
                if line.startswith('Ensembl Release'):
                    annotation_release = line.strip().split(' ')[2]
        os.remove(file_readme)

    ftp_directory_gtf = f'pub/release-{annotation_release}/gtf/{species_id}/'
    ftp_directory_fasta = f'pub/release-{annotation_release}/fasta/{species_id}/dna/'
    ftp_file_gtf = f'{species_id.capitalize()}.{genome_assembly}.{annotation_release}.gtf'
    ftp_file_fasta = f'{species_id.capitalize()}.{genome_assembly}.dna_rm.primary_assembly.fa'

    return ftp_link, ftp_directory_gtf, ftp_directory_fasta, ftp_file_gtf, ftp_file_fasta


def get_fasta(file_bed, file_reference_fasta, file_fasta, split=False):
    '''Get sequence for regions annotated in input gft file using genome fasta file.

    :param file_bed: Path to bed file with annotated genomic regions.
    :type file_bed: string
    :param file_reference_fasta: Path to fasta file with reference sequence, e.g. transcriptome.
    :type file_reference_fasta: string
    :param file_fasta: Path to fasta file where retrieved sequences are written to.
    :type file_fasta: string
    :param split: Use -split option of bedtools getfasta, defaults to False
    :type split: bool
    '''
    annotation = pybedtools.BedTool(file_bed)
    genome_sequence = pybedtools.BedTool(file_reference_fasta)

    annotation = annotation.sequence(fi=genome_sequence, s=True, name=True, split=split)
    annotation.save_seqs(file_fasta)
    