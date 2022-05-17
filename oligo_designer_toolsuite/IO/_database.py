import os
import shutil

import oligo_designer_toolsuite.IO._ftp_loader as ftp_loader


class BaseDB():
    '''_summary_
    '''
    def __init__(self, species) -> None:
        self.species = species # species of sequence, e.g. human, mouse
        self.blockSize = None
        self.file_DB = None
        pass


    def read_DB(self, file_DB):
        # check if file exists
        self.file_DB = file_DB
        pass


    def create_DB_from_NCBI_annotation(self, transcriptome=True, genome=True, genome_assembly='GRCh38', annotation_release='latest', dir_output='./annotations'):
        file_DB = os.path.join(dir_output, f'reference_DB_{self.species}_NCBI_{genome_assembly}_{annotation_release}_transcriptome{transcriptome}_genome{genome}')
        # check if DB exists
        # if file_DB exists:
        #   self.file_DB = file_DB
        #   return
        file_annotation = ftp_loader.download_gene_gtf_NCBI(dir_output)
        file_sequence = ftp_loader.download_genome_fasta_NCBI(dir_output)

        if transcriptome and not genome:
            file_transcriptome = self._get_transcriptome(file_sequence, file_annotation)
            shutil.move(file_transcriptome, file_DB)
        
        elif not transcriptome and genome:
            shutil.copyfile(file_sequence, file_DB)
        
        elif transcriptome and genome:
            file_transcriptome = self._get_transcriptome(file_sequence, file_annotation)
            with open(file_DB, 'w') as handle:
                shutil.copyfileobj(open(file_transcriptome,'rb'), handle)
                shutil.copyfileobj(open(file_sequence,'rb'), handle)
            shutil.rmtree(file_transcriptome)
            
        else:
            raise ValueError('<genomic> and <transcript> set to False. Could not create Reference Database.')

        self.file_DB = file_DB

    def create_DB_from_ensemble_annotation(self, transcriptome=True, genome=True, genome_assembly='GRCh38', annotation_release='104', dir_output='./annotations'):
        file_DB = os.path.join(dir_output, f'reference_DB_{self.species}_ensemble_{genome_assembly}_{annotation_release}_transcriptome{transcriptome}_genome{genome}')
        # check if DB exists
        # if file_DB exists:
        #   self.file_DB = file_DB
        #   return
        file_annotation = ftp_loader.download_gene_gtf_ensemble(dir_output)
        file_sequence = ftp_loader.download_genome_fasta_ensemble(dir_output)
        
        if transcriptome and not genome:
            file_transcriptome = self._get_transcriptome(file_sequence, file_annotation)
            shutil.move(file_transcriptome, file_DB)
        
        elif not transcriptome and genome:
            shutil.copyfile(file_sequence, file_DB)
        
        elif transcriptome and genome:
            file_transcriptome = self._get_transcriptome(file_sequence, file_annotation)
            with open(file_DB, 'w') as handle:
                shutil.copyfileobj(open(file_transcriptome,'rb'), handle)
                shutil.copyfileobj(open(file_sequence,'rb'), handle)
            shutil.rmtree(file_transcriptome)
            
        else:
            raise ValueError('<genomic> and <transcript> set to False. Could not create Reference Database.')

        self.file_DB = file_DB

    def create_DB_from_custom_annotation(self, transcriptome=True, genome=True, file_annotation=None, file_sequence=None, file_transcriptome=None, dir_output='./annotations'):
        file_DB = os.path.join(dir_output, f'reference_DB_{self.species}_custom_transcriptome{transcriptome}_genome{genome}')
        # check if DB exists
        # if file_DB exists:
        #   self.file_DB = file_DB
        #   return
        if transcriptome and not genome:
            if file_transcriptome:
                shutil.copyfile(file_transcriptome, file_DB)
                print(f'Take transcriptome annotation from: {file_transcriptome}')
            elif file_annotation and file_sequence:
                file_transcriptome = self._get_transcriptome(file_sequence, file_annotation)
                shutil.move(file_transcriptome, file_DB)
            else:
                raise ValueError('No custom files provided for transcriptome!')
            
        elif not transcriptome and genome:
            if file_sequence:
                shutil.copyfile(file_sequence, file_DB)
            else:
                raise ValueError('No custom files provided for genome!')

        elif transcriptome and genome:
            if file_transcriptome:
                print(f'Take transcriptome annotation from: {file_transcriptome}')
            elif file_annotation and file_sequence:
                file_transcriptome = self._get_transcriptome(file_sequence, file_annotation)
            else:
                raise ValueError('No custom files provided for transcriptome!')
            
            if file_sequence == None:
                raise ValueError('No custom files provided for genome!')
            
            with open(file_DB, 'w') as handle:
                shutil.copyfileobj(open(file_transcriptome,'rb'), handle)
                shutil.copyfileobj(open(file_sequence,'rb'), handle)
            
            if file_transcriptome == None:
                shutil.rmtree(file_transcriptome)
            
        else:
            raise ValueError('<genomic> and <transcript> set to False. Could not create Reference Database.')

        self.file_DB = file_DB

    def _get_transcriptome(self):
        pass




class ReferenceDB(BaseDB):
    '''_summary_
    '''
    def __init__(self, species) -> None:
        super().__init__(species)
        self.blockSize = 200 # maximum length of "long oligos" is 180bp





class OligoDB(BaseDB):
    '''_summary_
    '''
    def __init__(self, species) -> None:
        super().__init__(species)
        self.blockSize = 10 #change to probe length   

    def mask_prohibited_sequences():
        pass

    def mask_repeats():
        pass

    def filter_xyz():
        #wrap different filter from 'oligo_filter'
        pass
