import pandas as pd
from gtfparse import read_gtf
import pybedtools


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
    
    file_genome_sequence = './data/Homo_sapiens.GRCh38.104.dna_rm.primary_assembly.fa'
    file_gene_annotation = './data/Homo_sapiens.GRCh38.104.gtf'
    print('get exome annotation')
    file_exon_annotation = get_exons_annotation(file_gene_annotation)
    print('get exome fasta')
    file_exon_sequence = get_exons_fasta(file_exon_annotation, file_genome_sequence)



if __name__ == '__main__':
    
    main()