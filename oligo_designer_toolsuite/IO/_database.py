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

        self.species = 'human'
        self.genome_assembly = 'GRCh38'
        self.annotation_release = 'current'


    def read_DB(self, file_DB):
        if os.path.exists(file_DB):
            if data_parser.check_fasta_format(file_DB):
                self.file_DB = file_DB
            else:
                raise ValueError('Database has incorrect format!')
        else:
            raise ValueError('Database file does not exist!')
        pass


    def download_annotation(self, species, genome_assembly, annotation_source, annotation_release, dir_output):
        if species is None:
            warnings.warn(f'No species defined. Using default species {self.species}!')
            species = self.species

        if genome_assembly is None:
            warnings.warn(f'No genome assembly defined. Using default assembly {self.genome_assembly}!')
            genome_assembly = self.genome_assembly

        if annotation_release is None:
            warnings.warn(f'No ncbi annotation release defined. Using default release {self.annotation_release}!')
            annotation_release = self.annotation_release

        if annotation_source == 'NCBI':
            FTP = ftp_loader.FTPLoaderNCBI(species, genome_assembly, annotation_release)
        elif annotation_source == 'ensemble'
            FTP = ftp_loader.FtpLoaderEnsemble(species, genome_assembly, annotation_release)
        
        file_annotation = FTP.download_gtf(dir_output)
        file_sequence = FTP.download_fasta(dir_output)

        return file_annotation, file_sequence

    def _get_gene_transcript(self, file_sequence, file_annotation, blockSize, dir_output):
        pass

    '''
    def _get_gene(self, file_sequence, file_annotation, dir_output):
        pass

    def _get_gene_CDS(self, file_sequence, file_annotation, blockSize, dir_output):
        pass
    '''



class ReferenceDB(BaseDB):
    '''_summary_
    '''
    def __init__(self, blockSize = 200) -> None:
        super().__init__()

        self.file_DB = None
        self.blockSize = blockSize # maximum oligo length of "long oligos" is 180


    def create_DB_from_ncbi_annotation(self, species=None, genome_assembly=None, annotation_release=None, genome=False, gene_transcript=True, dir_output='./annotation'):

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        
        file_annotation, file_sequence = self.download_annotation(species, genome_assembly, 'NCBI', annotation_release, dir_output)
        files_fasta = self._get_fasta_files(self, genome, gene_transcript, file_sequence, file_annotation, self.blockSize, dir_output)

        self.file_DB = os.path.join(dir_output, f'reference_DB_{self.species}_{self.genome_assembly}_NCBI_release{self.annotation_release}_genome{genome}_gene_transcript{gene_transcript}')
        data_parser.merge_fasta(files_fasta, self.file_DB)


    def create_DB_from_ensemble_annotation(self, species=None, genome_assembly=None, annotation_release=None, genome=False, gene_transcript=True, dir_output='./annotation'):

        Path(dir_output).mkdir(parents=True, exist_ok=True)

        file_annotation, file_sequence = self.download_annotation(species, genome_assembly, 'ensemble', annotation_release, dir_output)
        files_fasta = self._get_fasta_files(self, genome, gene_transcript, file_sequence, file_annotation, self.blockSize, dir_output)

        self.file_DB = os.path.join(dir_output, f'reference_DB_{self.species}_{self.genome_assembly}_ensemble_release{self.annotation_release}_genome{genome}_gene_transcript{gene_transcript}')
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

        
    def generate_DB_from_ncbi_annotation(self, genes = None, region='gene_transcript', species=None, genome_assembly=None, annotation_release=None, dir_output='./annotation'):

        file_annotation, file_sequence = self.download_annotation(species, genome_assembly, 'NCBI', annotation_release, dir_output)
        file_region = self._create_target_regions(genes, region, file_annotation, file_sequence, dir_output)
        self.DB = self._generate_oligos(file_region)
        

    def generate_DB_from_ensemble_annotation(self, genes = None, region='gene_transcript', species=None, genome_assembly=None, annotation_release=None, dir_output='./annotation'):

        file_annotation, file_sequence = self.download_annotation(species, genome_assembly, 'ensemble', annotation_release, dir_output)
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
        #self.file_DB = ...
        # -> output DB as fasta file
        pass


    def _create_target_regions(self, genes, region, file_annotation, file_sequence, dir_output):
        if genes is None:
            genes = self._get_gene_list_from_annotation(file_annotation)

        if region == 'gene':
            file_region = self._get_gene(file_sequence, file_annotation, dir_output)

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

    '''
    def mask_prohibited_sequences(self):
        pass


    def mask_repeats(self):
        pass


    def filter_xyz(self):
        #wrap different filter from 'oligo_filter'
        pass
    '''