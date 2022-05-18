from ctypes import FormatError
import os
import shutil
import warnings
from pathlib import Path

import random

import oligo_designer_toolsuite.IO._ftp_loader as ftp_loader
import oligo_designer_toolsuite.IO._data_parser as data_parser


class BaseDB():
    '''_summary_
    '''
    def __init__(self) -> None:

        self.genome_assembly = 'GRCh38'
        self.annotation_release_ncbi = 'latest'
        self.annotation_release_ensemble = '104'


    def read_DB(self, file_DB):
        '''Read database from file, check if databse has correct format (fasta file).

        :param file_DB: _description_
        :type file_DB: _type_
        '''
        if os.path.exists(file_DB):
            if data_parser.check_fasta_format(file_DB):
                self.file_DB = file_DB
            else:
                raise ValueError('Database has incorrect format!')
        else:
            raise ValueError('Database file does not exist!')
        pass


    def download_ncbi_annotation(self, genome_assembly, annotation_release, dir_output):
        if genome_assembly is None:
            warnings.warn(f'No genome assembly defined. Using default assembly {self.genome_assembly}!')
            genome_assembly = self.genome_assembly

        if annotation_release is None:
            warnings.warn(f'No ncbi annotation release defined. Using default release {self.annotation_release_ncbi}!')
            annotation_release = self.annotation_release_ncbi

        file_annotation = ftp_loader.download_gene_gtf_ncbi(dir_output, genome_assembly, annotation_release)
        file_sequence = ftp_loader.download_genome_fasta_ncbi(dir_output, genome_assembly, annotation_release)

        return file_annotation, file_sequence


    def download_ensemble_annotation(self, genome_assembly, annotation_release, dir_output):
        if genome_assembly is None:
            warnings.warn(f'No genome assembly defined. Using default assembly {self.genome_assembly}!')
            genome_assembly = self.genome_assembly

        if annotation_release is None:
            warnings.warn(f'No ensemble annotation release defined. Using default release {self.annotation_release_ensemble}!')
            annotation_release = self.annotation_release_ensemble

        file_annotation = ftp_loader.download_gene_gtf_ensemble(dir_output, genome_assembly, annotation_release)
        file_sequence = ftp_loader.download_genome_fasta_ensemble(dir_output, genome_assembly, annotation_release)

        return file_annotation, file_sequence


    def _get_gene_sequence(self, file_sequence, file_annotation, dir_output):
        '''Generate fasta file with DNA sequence corresponding to whole gene annotation, including exons, introns, 5'UTR

        :param file_sequence: _description_
        :type file_sequence: _type_
        :param file_annotation: _description_
        :type file_annotation: _type_
        '''
        pass

    def _get_gene_transcript(self, file_sequence, file_annotation, blockSize, dir_output):
        pass
        

    def _get_gene_CDS(self, file_sequence, file_annotation, blockSize, dir_output):
        pass




class ReferenceDB(BaseDB):
    '''_summary_
    '''
    def __init__(self) -> None:
        super().__init__()

        self.file_DB = None

        self.blockSize = 200 # maximum oligo length of "long oligos" is 180


    def create_DB_from_ncbi_annotation(self, genome_assembly=None, annotation_release=None, genome=False, gene_transcript=True, dir_output='./annotation'):

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_DB = os.path.join(dir_output, f'reference_DB_{genome_assembly}_NCBI_release{annotation_release}_genome{genome}_gene_transcript{gene_transcript}')

        file_annotation, file_sequence = self.download_ncbi_annotation(genome_assembly, annotation_release, dir_output)
        files_fasta = self._get_fasta_files(self, genome, gene_transcript, file_sequence, file_annotation, self.blockSize, dir_output)
        data_parser.merge_fasta(files_fasta, self.file_DB)


    def create_DB_from_ensemble_annotation(self, genome_assembly=None, annotation_release=None, genome=False, gene_transcript=True, dir_output='./annotation'):

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_DB = os.path.join(dir_output, f'reference_DB_{genome_assembly}_NCBI_release{annotation_release}_genome{genome}_gene_transcript{gene_transcript}')

        file_annotation, file_sequence = self.download_ensemble_annotation(genome_assembly, annotation_release, dir_output)
        files_fasta = self._get_fasta_files(self, genome, gene_transcript, file_sequence, file_annotation, self.blockSize, dir_output)
        data_parser.merge_fasta(files_fasta, self.file_DB)


    def create_DB_from_custom_annotation(self, file_annotation, file_sequence, genome=False, gene_transcript=True, dir_output='./annotations'):

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_DB = os.path.join(dir_output, f'reference_DB_custom_genome{genome}_gene_transcript{gene_transcript}')

        if not data_parser.check_gtf_format(file_annotation):
            raise ValueError('Annotation File has incorrect format!')
            os._exit(0)
        if not data_parser.check_fasta_format(file_sequence):
            raise ValueError('Sequence File has incorrect format!')
            os._exit(0)

        files_fasta = self._get_fasta_files(self, genome, gene_transcript, file_sequence, file_annotation, self.blockSize, dir_output)
        data_parser.merge_fasta(files_fasta, self.file_DB)
        

    def _get_fasta_files(self, genome, gene_transcript, file_sequence, file_annotation, blockSize, dir_output):

        files_fasta = []

        if genome:
            file_genome = os.path.join(dir_output,'genome.fna')
            shutil.copyfile(file_sequence, file_genome)
            files_fasta.append(file_genome)

        if gene_transcript:
            file_gene_transcript = self._get_gene_transcript(file_sequence, file_annotation, blockSize, dir_output)
            files_fasta.append(file_gene_transcript)

        return files_fasta



class OligoDB(BaseDB):
    '''_summary_
    '''
    def __init__(self, probe_length_min, probe_length_max) -> None:
        super().__init__()

        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max
        self.blockSize = probe_length_max

        self.file_DB = None
        self.DB = None

        self.annotation_source = None

        
    def generate_DB_from_ncbi_annotation(self, genes = None, region='gene_transcript', genome_assembly=None, annotation_release=None, dir_output='./annotation'):

        self.annotation_source = 'ncbi'

        file_annotation, file_sequence = self.download_ncbi_annotation(genome_assembly, annotation_release, dir_output)
        file_region = self._create_target_regions(genes, region, file_annotation, file_sequence, dir_output)
        self.DB = self._generate_oligos(file_region)
        

    def generate_DB_from_ensemble_annotation(self, genes = None, region='gene_transcript', genome_assembly=None, annotation_release=None, dir_output='./annotation'):

        self.annotation_source = 'ensemble'

        file_annotation, file_sequence = self.download_ensemble_annotation(genome_assembly, annotation_release, dir_output)
        file_region = self._create_target_regions(genes, region, file_annotation, file_sequence, dir_output)
        self.DB = self._generate_oligos(file_region)


    def generate_DB_from_custom_annotation(self, file_annotation, file_sequence, genes = None, region='gene_transcript', dir_output='./annotation'):

        self.annotation_source = 'custom'
        
        if not data_parser.check_gtf_format(file_annotation):
            raise ValueError('Annotation File has incorrect format!')
            os._exit(0)
        if not data_parser.check_fasta_format(file_sequence):
            raise ValueError('Sequence File has incorrect format!')
            os._exit(0)

        file_region = self._create_target_regions(genes, region, file_annotation, file_sequence, dir_output)
        self.DB = self._generate_oligos(file_region)


    def write_DB(self, dir_output):
        self.file_DB = os.path.join(dir_output, f'oligo_DB')
        # -> output DB as fasta file
        pass


    def mask_prohibited_sequences(self):
        pass


    def mask_repeats(self):
        pass


    def filter_xyz(self):
        #wrap different filter from 'oligo_filter'
        pass


    def _create_target_regions(self, genes, region, file_annotation, file_sequence, dir_output):
        if genes is None:
            genes = self._get_gene_list_from_annotation(file_annotation)

        if region == 'gene_sequence':
            file_region = self._get_gene_sequence(file_sequence, file_annotation, dir_output)

        if region == 'gene_transcript':
            file_region = self._get_gene_transcript(file_sequence, file_annotation, self.blockSize, dir_output)
            
        if region == 'gene_CDS':
            file_region = self._get_gene_CDS(file_sequence, file_annotation, self.blockSize, dir_output)

        return file_region


    def _get_gene_list_from_annotation(self, file_annotation):
        gene_annotation = data_parser.read_gtf(file_annotation)
        genes = gene_annotation.loc[gene_annotation['feature'] == 'gene']
        genes = list(genes['gene_id'].unique())
        random.shuffle(genes)
            
        return genes

    def _generate_oligos(self, file_region):
        pass