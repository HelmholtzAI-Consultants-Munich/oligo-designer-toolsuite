############################################
# imports
############################################

import os
import time
import pickle
import subprocess as sp
import pandas as pd


import Levenshtein as Lev

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

############################################
# filter probes with blastn or hamming distance 
############################################


def filter_probes_with_blastn(dir_output, file_exon_sequence, mapping_gene_to_probes_unique, mapping_exon_to_gene, min_probes_per_gene, probe_length, word_size, coverage, percent_identity, num_threads):
  
    #print_run = 500
    count_genes = 0
    count_probes = 0
    count_probes_last_run = 0
    count_selected_probes = 0

    # create temporary directory
    dir_output_blast = '{}probe_queries_combined/'.format(dir_output)
    cmd = 'mkdir {}'.format(dir_output_blast)
    sp.run(cmd, shell=True)

    # create blast database
    cmd = NcbimakeblastdbCommandline(input_file=file_exon_sequence, dbtype='nucl')
    out, err = cmd()

    t = time.time()
    mapping_gene_to_probes_blastn = {}

    probes = []
    for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_unique.items():
        count_genes += 1
        #if count_genes > print_run:
        #    break
        
        # exclude genes with less than required number of probes
        if len(mapping_probe_to_exon) > min_probes_per_gene:
            
            # write all probes to fasta file
            for probe, exons in mapping_probe_to_exon.items():
                count_probes += 1
                sequence = SeqRecord(probe, '{}_{}'.format(gene_id, probe), '', '')
                probes.append(sequence)

        else:
            mapping_gene_to_probes_blastn[gene_id] = None

    file_fasta_probes_per_gene = '{}probes.fna'.format(dir_output_blast)
                
    output = open(file_fasta_probes_per_gene, 'w')
    SeqIO.write(probes, output, 'fasta')
    output.close()


    # run BlastN
    file_blast_output = '{}blast.txt'.format(dir_output_blast)
    cmd = NcbiblastnCommandline(query=file_fasta_probes_per_gene,db=file_exon_sequence, outfmt=10, out=file_blast_output,word_size=word_size, strand='plus',num_threads=num_threads)
    out, err = cmd()

    # read and process results of BlastN
    blast_results = pd.read_csv(file_blast_output,header=None, usecols=[0,1,2,3,6,7]) # don't load columns: 'missmatches','gap_opens','start_t','end_t', 'evalue','bit_score'
    blast_results.columns = ['query','target_exon','percent_identity','alignment_length','query_start','query_end']
    
    blast_results['query_gene_id'] = blast_results['query'].str.split('_').str[0]
    blast_results['query_probe'] = blast_results['query'].str.split('_').str[1]
    blast_results['target_exon'] = blast_results['target_exon'].str.split('::').str[0]
    
    
    for gene_id in blast_results['query_gene_id'].unique():
        probes_with_match = []
        blast_results_gene = blast_results[blast_results['query_gene_id'] == gene_id]
        blast_results_gene.reset_index(inplace=True, drop=True)
        for idx in blast_results_gene.index:
            probe = blast_results_gene['query_probe'][idx]
            gene_id_target = mapping_exon_to_gene[blast_results_gene['target_exon'][idx]]

            # only consider hit in other genes
            if gene_id_target != gene_id:
                # filter out probes that have more than 50% sequence coverage with target sequence, 80% identity with target sequence, and target sequence covers ligation site +- 5
                if int(blast_results_gene['alignment_length'][idx]) > (probe_length * coverage) and float(blast_results_gene['percent_identity'][idx]) > percent_identity and int(blast_results_gene['query_start'][idx]) < (probe_length //2 - 4) and int(blast_results_gene['query_end'][idx]) > (probe_length //2 + 5):
                    probes_with_match.append(probe)

        # remove all probes with matches
        mapping_probe_to_exon = mapping_gene_to_probes_unique[gene_id]
        for probe in set(probes_with_match):
            mapping_probe_to_exon.pop(probe)
        mapping_gene_to_probes_blastn[gene_id] = mapping_probe_to_exon
        count_selected_probes += len(mapping_probe_to_exon)
    
            
    count_probes_single_run = count_probes - count_probes_last_run
    time_elapsed = time.time() - t

    print('In total, {} genes with {} probes were processed and {} probes passed the blastn filter in total.'.format(count_genes, count_probes, count_selected_probes))
    print('In one run, {} genes with {} probes were processed, which took {} s and {} s per probe'.format(count_genes,count_probes_single_run,time_elapsed,time_elapsed/count_probes_single_run))

    count_probes_last_run = count_probes
    t = time.time()

    # remove blast files
    cmd = 'rm -rf {}'.format(dir_output_blast)
    sp.run(cmd, shell=True)

    return mapping_gene_to_probes_blastn


############################################

'''
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


'''




