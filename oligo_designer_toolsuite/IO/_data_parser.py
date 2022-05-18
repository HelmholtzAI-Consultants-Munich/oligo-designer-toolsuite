
import os
import shutil


def read_gtf():
    pass

def read_gff():
    pass

def get_sequence_from_annotation():
    pass

def check_gtf_format():
    pass

def check_fasta_format():
    pass

def merge_fasta(files_fasta, file_merged_fasta):

    if files_fasta == []:
        raise ValueError('No fasta files provided for merge.')

    with open(file_merged_fasta, 'w') as handle_DB:
        for file in files_fasta:
            if os.path.exists(file):
                shutil.copyfileobj(open(file,'rb'), handle_DB)
                shutil.rmtree(file)
            else:
                raise ValueError(f'Fasta file {file} does not exist!')