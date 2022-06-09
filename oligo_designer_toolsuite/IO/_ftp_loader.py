
#import pybedtools

import os
import re
import gzip
import shutil

from ftplib import FTP

class BaseFtpLoader():
    '''_summary_
    '''
    def __init__(self, ftp_link, ftp_directory, file_name, dir_output) -> None: 
        self.ftp_link=ftp_link, 
        self.ftp_directory=ftp_directory 
        self.file_name=file_name
        self.dir_output=dir_output

    def ftp_download(self):

        '''
        Download file from ftp server.
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
        ###Code####
        ftp = FTP(self.ftp_link)
        ftp.login()  # login to ftp server
        ftp.cwd(self.ftp_directory)  # move to directory

        files = ftp.nlst()

        for file in files:
            if re.match(self.file_name, file):
                file_output = os.path.join(self.dir_output, file)
                ftp.retrbinary('RETR ' + file, open(file_output, 'wb').write)

        ftp.quit()

        return file_output
  
        
        pass 

    def decompress_gzip(self, file_gzip):
        '''
        Decompress zip files.
        :param file_gzip: Path to zipped file.
        :type file_gzip: string
        :return: Path to unzipped file.
        :rtype: string
        
        ###Code###
        file_output = file_gzip.split('.gz')[0]
        with gzip.open(file_gzip, 'rb') as f_in:
            with open(file_output, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(file_gzip)
        return file_output
        '''
        pass

    def download_gtf(self):
        '''
        :param dir_output: _description_
        :type dir_output: _type_
        :return: _description_
        :rtype: _type_
        '''

        ### TO DO: Generalize gtf download
    
        file_readme = self.ftp_download(self.ftp_link, self.ftp_directory, 'README', self.dir_output)
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
        ftp_file_chr_mapping = f'{assembly_accession}_{assembly_name}_assembly_report.txt'

        return ftp_link, ftp_directory, ftp_file_gtf, ftp_file_chr_mapping

        
        pass

    def download_fasta(self, dir_output):
        '''
        :param dir_output: _description_
        :type dir_output: _type_
        :return: _description_
        :rtype: _type_
        '''

        ### TO DO: Generalize gtf download
        
        ftp_link = self.generate_FTP_link()
        file_readme = self.ftp_download(ftp_link, ftp_directory, 'README', dir_output)
        with open(file_readme, 'r') as handle:
            for line in handle:
                if line.startswith('ASSEMBLY NAME:'):
                    assembly_name = line.strip().split('\t')[1]
                if line.startswith('ASSEMBLY ACCESSION:'):
                    assembly_accession = line.strip().split('\t')[1]
                    break
        os.remove(file_readme)
        ftp_directory = ftp_directory + f'{assembly_accession}_{assembly_name}'
        ftp_file_fasta = f'{assembly_accession}_{assembly_name}_genomic.fna.gz'

        return ftp_link, ftp_directory, ftp_file_fasta




class FtpLoaderEnsemble(BaseFtpLoader):

    def __init__(self, species, genome_assembly, annotation_release) -> None:
        super().__init__()
        self.species = species
        self.genome_assembly = genome_assembly
        self.annotation_release = annotation_release
        
    def generate_FTP_link():
        return 'ftp.ensembl.org'


class FTPLoaderNCBI(BaseFtpLoader):

    def __init__(self, species, genome_assembly, annotation_release) -> None:
        '''_summary_
        :param species: _description_
        :type species: _type_
        :param annotation_release: _description_
        :type annotation_release: _type_
        :param genome_assembly: _description_
        :type genome_assembly: _type_
        '''
        super().__init__()
        self.species = species
        self.genome_assembly = genome_assembly
        self.annotation_release = annotation_release

        
        if species == 'human':
            ftp_directory = 'refseq/H_sapiens/annotation/annotation_releases/'
        if species == 'mouse':
            ftp_directory = 'refseq/M_musculus/annotation_releases/'

        if annotation_release == 'current':
            ftp_directory = ftp_directory + 'current/'
        else:
            ftp_directory = ftp_directory + f'{annotation_release}/'

        self.ftp_directory=ftp_directory
        
        
    def generate_FTP_link(self):
        return 'ftp.ncbi.nlm.nih.gov'
        

    def _download_mapping_chr_names(self):
        '''Download file with mapping of chromosome names between GenBank and Ref-Seq accession number 
        from ftp server and create a mapping dictionary.
        :return: Dictionary with mapping of chromsome names from GenBank to Ref-Seq.
        :rtype: dict
        '''

        ##Use _download_chr_mapping in annotation_loader

        pass

        


    def _map_chr_names_gene_gtf():

        ### process_ncbi_gene_gtf in annotation_loader
        pass

    def _map_chr_names_genome_fasta():

        ### process_ncbi_genome_fasta in annotation_loader
        pass

    def download_gtf(self):

        ### download_gene_gtf in annotation_loader

        pass

    def download_fasta(self):

        ### download_genome_fasta in annotation_loader
   
        pass


