############################################
# imports
############################################

import os
import time
import pickle
import argparse
import pandas as pd

from gtfparse import read_gtf
import Levenshtein as Lev
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt

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
    args_parser.add_argument('-fe','--file_exon_sequence',help='exon sequence file (fasta file)',type=str,required=True)
    args_parser.add_argument('-fg','--file_gene_annotation',help='gene annotation file (gtf file)',type=str,required=True)
    args_parser.add_argument('-o','--dir_output',help='output directory',type=str,required=True)
    args_parser.add_argument('-len','--probe_length',help='length of probes',type=int,required=False, default='45')
    args_parser.add_argument('-Pmin','--min_probes_per_gene',help='minimum number of probes per gene',type=int,required=False, default='4')
    args_parser.add_argument('-GCmin','--GC_content_min',help='minimum GC content of probes',type=float,required=False, default='40')
    args_parser.add_argument('-GCmax','--GC_content_max',help='maximum GC content of probes',type=float,required=False, default='60')
    args_parser.add_argument('-Tmin','--Tm_min',help='minimum melting temperature of probes',type=float,required=False, default='68')
    args_parser.add_argument('-TCmax','--Tm_max',help='maximum melting temperature of probes',type=float,required=False, default='75')
    args_parser.add_argument('-c','--coverage',help='minimum coverage between probes and target sequence (blast parameter "coverage") given in percent (%), ranging from 0 tp 100.',type=float,required=False, default='50')
    args_parser.add_argument('-pe','--percent_identity',help='Maximum similarity between probes and target sequences (blast parameter "percent identity") given in percent (%), ranging from 0 tp 100.',type=float,required=False, default='80')

    return args_parser.parse_args()



def get_mapping_gene_exon(dir_output, file_gene_annotation):

    file_mapping_exon_to_gene = '{}mapping_exon_to_gene.pkl'.format(dir_output)
    file_mapping_gene_to_exon = '{}mapping_gene_to_exon.pkl'.format(dir_output)

    if not os.path.isfile(file_mapping_exon_to_gene) or not os.path.isfile(file_mapping_gene_to_exon):

        gene_annotation = read_gtf(file_gene_annotation)
        exon_annotation = gene_annotation[gene_annotation["feature"] == "exon"]

        mapping_exon_to_gene = pd.Series(exon_annotation.gene_id.values,index=exon_annotation.exon_id).to_dict()
        mapping_gene_to_exon = pd.Series(exon_annotation.exon_id.values,index=exon_annotation.gene_id).to_dict()

        pickle.dump(mapping_exon_to_gene, open(file_mapping_exon_to_gene,'wb'))
        pickle.dump(mapping_gene_to_exon, open(file_mapping_gene_to_exon,'wb'))

    return file_mapping_exon_to_gene, file_mapping_gene_to_exon



def get_mapping_gene_probes(dir_output, file_exon_sequence, mapping_exon_to_gene, probe_length, GC_content_min, GC_content_max, Tm_min, Tm_max):

    file_mapping_gene_to_probes = '{}mapping_gene_to_probes.pkl'.format(dir_output)

    if not os.path.isfile(file_mapping_gene_to_probes):
        count_exons = 0
        count_probes = 0

        mapping_gene_to_probes = {}

        for exon in SeqIO.parse(file_exon_sequence, "fasta"):

            count_exons += 1
            exon_id = exon.id.split('::')[0]
            gene_id = mapping_exon_to_gene[exon_id]

            exon_sequence = exon.seq
            if len(exon_sequence) > probe_length:
                probes_of_exon = set([exon_sequence[i:i+probe_length] for i in range(len(exon_sequence)-(probe_length-1)) if ('N' not in exon_sequence[i:i+probe_length])])
                count_probes += len(probes_of_exon)
                
                for probe in probes_of_exon:
                    gc_content = GC(probe)
                    if (GC_content_min < gc_content < GC_content_max):
                        Tm = mt.Tm_NN(probe)
                        if (Tm_min < Tm < Tm_max):
                            if gene_id in mapping_gene_to_probes:
                                if probe in mapping_gene_to_probes[gene_id]:
                                    mapping_gene_to_probes[gene_id][probe].append(exon_id)
                                else:
                                    mapping_gene_to_probes[gene_id][probe] = [exon_id]
                            else:
                                mapping_gene_to_probes[gene_id] = {probe: [exon_id]}
                        
            #if count_exons % 100000 == 0:
            #    print('Processed exons: {}'.format(count_exons))
            #    print('Probes found in total: {}'.format(count_probes))
            #    break

        print('Total number of exons processed: {} with {} probes in total.'.format(count_exons, count_probes))
        pickle.dump(mapping_gene_to_probes, open(file_mapping_gene_to_probes,'wb'))

    return file_mapping_gene_to_probes



def get_unique_probes(dir_output, mapping_gene_to_probes, mapping_gene_to_probes_unique):

    file_mapping_gene_to_probes_unique = '{}mapping_gene_to_probes_unique.pkl'.format(dir_output)

    if not os.path.isfile(file_mapping_gene_to_probes_unique):

        count_genes = 0
        count_probes = 0

        probes_unique = {}
        probes_non_unique = {}

        for gene_id, mapping_probe_to_exon in mapping_gene_to_probes.items():
            count_genes +=1

            for probe, exon_id in mapping_probe_to_exon.items():
                count_probes +=1

                if probe in probes_non_unique:
                    probes_non_unique[probe] += 1
                    mapping_gene_to_probes_unique[gene_id].pop(probe)

                elif probe in probes_unique:
                    mapping_gene_to_probes_unique[gene_id].pop(probe)
                    mapping_gene_to_probes_unique[probes_unique[probe][1]].pop(probe)
                    probes_unique.pop(probe)
                    probes_non_unique[probe] = 2

                else:
                    probes_unique[probe] = [exon_id, gene_id]

                #if count_probes % 100000 == 0:
                #    print('Processed genes: {}'.format(count_genes))
                #    print('Processed probes: {}'.format(count_probes))
                #    print('Unique probes: {}'.format(len(probes_unique)))

        print('Total number of genes processed: {}, total number of probes processed: {} with {} unique probes in total.'.format(count_genes, count_probes, len(probes_unique)))    
        pickle.dump(mapping_gene_to_probes_unique, open(file_mapping_gene_to_probes_unique,'wb')) 
    
    return file_mapping_gene_to_probes_unique



def filter_probes_with_hamming_distance(dir_output, mapping_gene_to_probes_unique, min_probes_per_gene, probe_length, percent_identity):
  
    file_mapping_gene_to_probes_hamming = '{}mapping_gene_to_probes_hamming.pkl'.format(dir_output)

    if not os.path.isfile(file_mapping_gene_to_probes_hamming):

        count_genes = 0
        count_probes = 0
        count_selected_probes = 0

        mapping_gene_to_probes_hamming = {}
        threshold = probe_length / 100 * (100-percent_identity)

        print('Filter out probes that map to target sequences with < {} missmatches.'.format(threshold))

        for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_unique.items():

            count_genes +=1

            if len(mapping_probe_to_exon) > min_probes_per_gene:
                
                for probe, exons in mapping_probe_to_exon.items():
                    count_probes += 1
                    foundMatch = False

                    for gene_id2, mapping_probe_to_exon2 in mapping_gene_to_probes_unique.items():

                            if gene_id != gene_id2:
                                
                                for probe2, exons2 in mapping_probe_to_exon2.items():
                                    hamming_dist = Lev.hamming(str(probe),str(probe2))
                                    if  hamming_dist < threshold:
                                            foundMatch = True
                                            break
                            
                            if foundMatch:
                                break

                    if not foundMatch:   
                        count_selected_probes +=1

                        if gene_id not in mapping_gene_to_probes_hamming:
                            mapping_gene_to_probes_hamming[gene_id] = []

                        mapping_gene_to_probes_hamming[gene_id].append(probe)

                    #if count_probes % 100000 == 0:
                    #    print('Processed genes: {}'.format(count_genes))
                    #    print('Processed probes: {}'.format(count_probes))
                    #    print('Probes passed hamming distance filter: {}'.format(count_selected_probes))
                
        print('Total number of genes processed: {}, total number of probes processed: {} with {} probes passed hamming distance filter.'.format(count_genes, count_probes, count_selected_probes))    
        pickle.dump(mapping_gene_to_probes_hamming, open(file_mapping_gene_to_probes_hamming,'wb')) 
    
    return file_mapping_gene_to_probes_hamming



def blastn_filter(gene_id, mapping_probe_to_exon, mapping_gene_to_probes_unique, min_probes_per_gene, mapping_gene_to_probes_blastn):

    if len(mapping_probe_to_exon) > min_probes_per_gene:
                
        for probe, exons in mapping_probe_to_exon.items():
            
            foundMatch = False

            for gene_id2, mapping_probe_to_exon2 in mapping_gene_to_probes_unique.items():

                    if gene_id != gene_id2:
                        
                        for probe2, exons2 in mapping_probe_to_exon2.items():
                            cmd = 'blastn -query {} -db {} -outfmt 10 -out blast_output.txt -word_size 7 -strand plus'
                            output1 = sp.check_output(cmd, shell=True)
                            #hamming_dist = Lev.hamming(str(probe),str(probe2))
                            #if  hamming_dist < probe_length / 100 * percent_identity:
                            #        foundMatch = True
                            #        break
                    
                    if foundMatch:
                        break

            if not foundMatch:   

                if gene_id not in mapping_gene_to_probes_blastn:
                    mapping_gene_to_probes_blastn[gene_id] = []

                mapping_gene_to_probes_blastn[gene_id].append(probe)

            


def filter_probes_with_blastn(dir_output, mapping_gene_to_probes_unique, min_probes_per_gene, probe_length, percent_identity):
  
    file_mapping_gene_to_probes_blastn = '{}mapping_gene_to_probes_blastn.pkl'.format(dir_output)

    if not os.path.isfile(file_mapping_gene_to_probes_blastn):

        count_genes = 0
        count_probes = 0
        count_selected_probes = 0

        mapping_gene_to_probes_blastn = {}

        for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_unique.items():
            count_genes +=1
            blastn_filter(gene_id, mapping_probe_to_exon, mapping_gene_to_probes_unique, min_probes_per_gene, mapping_gene_to_probes_blastn)
            
                
        print('Total number of genes processed: {}, total number of probes processed: {} with {} probes passed hamming distance filter.'.format(count_genes, count_probes, count_selected_probes))    
        pickle.dump(mapping_gene_to_probes_blastn, open(file_mapping_gene_to_probes_blastn,'wb')) 
    
    return file_mapping_gene_to_probes_blastn



def main():
    
    parameters = args()

    file_exon_sequence = parameters.file_exon_sequence
    file_gene_annotation = parameters.file_gene_annotation
    dir_output = parameters.dir_output
    probe_length = parameters.probe_length
    min_probes_per_gene = parameters.min_probes_per_gene
    GC_content_min = parameters.GC_content_min
    GC_content_max = parameters.GC_content_max
    Tm_min = parameters.Tm_min
    Tm_max = parameters.Tm_max
    Tm_max = parameters.Tm_max
    coverage = parameters.coverage
    percent_identity = parameters.percent_identity

    print('\nParameter setting:')
    print('Exon sequence file: {}'.format(file_exon_sequence))
    print('Gene annotation file: {}'.format(file_gene_annotation))
    print('Output directory: {}'.format(dir_output))
    print('Length of probes: {}'.format(probe_length))
    print('Minimum number of probes per gene: {}'.format(min_probes_per_gene))
    print('Minimum GC content of probes: {}'.format(GC_content_min))
    print('Maximum GC content of probes: {}'.format(GC_content_max))
    print('Minimum melting temperature of probes: {}'.format(Tm_min))
    print('Maximum melting temperature of probes: {}'.format(Tm_max))
    print('Minimum coverage between probes and target sequence (blast parameter "coverage"): {} %'.format(coverage))
    print('Maximum similarity between probes and target sequences (blast parameter "percent identity"): {} % \n'.format(percent_identity))

    t = time.time()
    file_mapping_exon_to_gene, file_mapping_gene_to_exon = get_mapping_gene_exon(dir_output, file_gene_annotation)
    print('Time to process gene-exon mapping dicts: {} s'.format(time.time() - t))

    mapping_exon_to_gene = pickle.load(open(file_mapping_exon_to_gene,'rb'))
    #mapping_gene_to_exon = pickle.load(open(file_mapping_gene_to_exon,'rb'))

    t = time.time()
    file_mapping_gene_to_probes = get_mapping_gene_probes(dir_output, file_exon_sequence, mapping_exon_to_gene, probe_length, GC_content_min, GC_content_max, Tm_min, Tm_max)
    print('Time to process gene-probe mapping dict: {} s'.format(time.time() - t))

    mapping_gene_to_probes = pickle.load(open(file_mapping_gene_to_probes,'rb'))
    mapping_gene_to_probes_unique = pickle.load(open(file_mapping_gene_to_probes,'rb'))

    t = time.time()
    file_mapping_gene_to_probes_unique = get_unique_probes(dir_output, mapping_gene_to_probes, mapping_gene_to_probes_unique)
    print('Time to process unique probes: {} s'.format(time.time() - t))

    mapping_gene_to_probes_unique = pickle.load(open(file_mapping_gene_to_probes_unique,'rb'))
    removed_genes = []
    for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_unique.items():
        if len(mapping_probe_to_exon) < min_probes_per_gene:
            removed_genes.append(gene_id)

    print('{} genes were filtered out because they had < {} probes per gene.'.format(len(removed_genes), min_probes_per_gene))

    t = time.time()
    file_mapping_gene_to_probes_hamming = filter_probes_with_hamming_distance(dir_output, mapping_gene_to_probes_unique, min_probes_per_gene, probe_length, percent_identity)
    print('Time to process hamming distance filter: {} s'.format(time.time() - t))
    #mapping_gene_to_probes_hamming = pickle.load(open(file_mapping_gene_to_probes_hamming,'rb'))



if __name__ == '__main__':
    
    main()