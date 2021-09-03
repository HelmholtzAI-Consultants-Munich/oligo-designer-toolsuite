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

############################################
# filter probes with blastn or hamming distance 
############################################


def filter_probes_with_blastn(dir_output, file_exon_sequence, mapping_gene_to_probes_unique, mapping_exon_to_gene, min_probes_per_gene, probe_length, coverage, percent_identity, num_threads):
  
    file_mapping_gene_to_probes_blastn = '{}mapping_gene_to_probes_blastn.pkl'.format(dir_output)

    if not os.path.isfile(file_mapping_gene_to_probes_blastn):
        print_run = 10
        count_genes = 0
        count_probes = 0
        count_probes_last_run = 0
        count_selected_probes = 0
        
        dir_output_blast = '{}probe_queries_single/'.format(dir_output)
        cmd = 'mkdir {}'.format(dir_output_blast)
        sp.run(cmd, shell=True)

        cmd = 'makeblastdb -in {} -dbtype nucl'.format(file_exon_sequence)
        sp.run(cmd, shell=True)

        t = time.time()
        mapping_gene_to_probes_blastn = {}

        for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_unique.items():
            if count_genes > 1000:
                return
            count_genes += 1
            
            # exclude genes with less than required number of probes
            if len(mapping_probe_to_exon) > min_probes_per_gene:
                
                # write all probes to fasta file
                probes = []
               
                for probe, exons in mapping_probe_to_exon.items():
                    count_probes += 1
                    sequence = SeqRecord(probe, '{}'.format(probe), '', '')
                    probes.append(sequence)

                file_fasta_probes_per_gene = '{}probes_{}.fna'.format(dir_output_blast, gene_id)
                
                output = open(file_fasta_probes_per_gene, 'w')
                SeqIO.write(probes, output, 'fasta')
                output.close()

                # run BlastN
                file_blast_output = '{}blast_{}.txt'.format(dir_output_blast, gene_id)
                cmd = 'blastn -query {} -db {} -outfmt 10 -out {} -word_size 7 -strand plus -num_threads {} 2>/dev/null'.format(file_fasta_probes_per_gene, file_exon_sequence, file_blast_output, num_threads)
                sp.run(cmd, shell=True)

                # read and process results of BlastN
                blast_results = pd.read_csv(file_blast_output,header=None, usecols=[0,1,2,3,6,7]) # don't load columns: 'missmatches','gap_opens','start_t','end_t', 'evalue','bit_score'
                blast_results.columns = ['probe_query','exon_target','percent_identity','alignment_length','start_q','end_q']
                blast_results['exon_target'] = blast_results['exon_target'].str.split('::').str[0]
                
                probes_with_match = []
                for idx in blast_results.index:
                    probe = blast_results['probe_query'][idx]
                    gene_id_target = mapping_exon_to_gene[blast_results['exon_target'][idx]]

                    # only consider hit in other genes
                    if gene_id_target != gene_id:
                        # filter out probes that have more than 50% sequence coverage with target sequence, 80% identity with target sequence, and target sequence covers ligation site +- 5
                        if int(blast_results['alignment_length'][idx]) > (probe_length * coverage) and float(blast_results['percent_identity'][idx]) > percent_identity and int(blast_results['start_q'][idx]) < (probe_length //2 - 4) and int(blast_results['end_q'][idx]) > (probe_length //2 + 5):
                            probes_with_match.append(probe)

                # remove all probes with matches
                for probe in set(probes_with_match):
                    mapping_probe_to_exon.pop(probe)
                count_selected_probes += len(mapping_probe_to_exon)
                mapping_gene_to_probes_blastn[gene_id] = mapping_probe_to_exon

            
            if count_genes % print_run == 0:
                count_probes_single_run = count_probes - count_probes_last_run
                time_elapsed = time.time() - t

                print('In total, {} genes with {} probes were processed and {} probes passed the blastn filter in total.'.format(count_genes, count_probes, count_selected_probes))
                print('In one run, {} genes with {} probes were processed, which took {} s and {} s per probe'.format(print_run,count_probes_single_run,time_elapsed,time_elapsed/count_probes_single_run))
               
                count_probes_last_run = count_probes
                t = time.time()

    print('Total number of genes processed: {}, total number of probes processed: {} with {} probes passed hamming distance filter.'.format(count_genes, count_probes, count_selected_probes))    
    pickle.dump(mapping_gene_to_probes_blastn, open(file_mapping_gene_to_probes_blastn,'wb')) 

    # remove blast files
    cmd = 'rm -rf {}'.format(dir_output_blast)
    sp.run(cmd, shell=True)

    return file_mapping_gene_to_probes_blastn



############################################

def filter_probes_with_blastn_2(dir_output, file_exon_sequence, mapping_gene_to_probes_unique, mapping_exon_to_gene, min_probes_per_gene, probe_length, coverage, percent_identity, num_threads):
  
    file_mapping_gene_to_probes_blastn = '{}mapping_gene_to_probes_blastn.pkl'.format(dir_output)

    if not os.path.isfile(file_mapping_gene_to_probes_blastn):
        count_genes = 0
        count_probes = 0
        count_selected_probes = 0

        dir_output_blast = '{}probe_queries_joint/'.format(dir_output)
        cmd = 'mkdir {}'.format(dir_output_blast)
        sp.run(cmd, shell=True)

        cmd = 'makeblastdb -in {} -dbtype nucl'.format(file_exon_sequence)
        sp.run(cmd, shell=True)

        mapping_gene_to_probes_blastn = {}

        probes = []
        for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_unique.items():
            if count_genes > 10:
                break
            count_genes += 1
            
            # exclude genes with less than required number of probes
            if len(mapping_probe_to_exon) > min_probes_per_gene:
                
                # write all probes to fasta file
                
                for probe, exons in mapping_probe_to_exon.items():
                    count_probes += 1
                    sequence = SeqRecord(probe, '{}_{}'.format(gene_id, probe), '', '')
                    probes.append(sequence)

        file_fasta_probes_per_gene = '{}probes.fna'.format(dir_output_blast)
                
        output = open(file_fasta_probes_per_gene, 'w')
        SeqIO.write(probes, output, 'fasta')
        output.close()

        # run BlastN
        file_blast_output = '{}blast.txt'.format(dir_output_blast)
        cmd = 'blastn -query {} -db {} -outfmt 10 -out {} -word_size 7 -strand plus -num_threads {} 2>/dev/null'.format(file_fasta_probes_per_gene, file_exon_sequence, file_blast_output, num_threads)
        print(cmd)
        sp.run(cmd, shell=True)

        # read and process results of BlastN
        blast_results = pd.read_csv(file_blast_output,header=None, usecols=[0,1,2,3,6,7]) # don't load columns: 'missmatches','gap_opens','start_t','end_t', 'evalue','bit_score'
        blast_results.columns = ['query','exon_target','percent_identity','alignment_length','start_q','end_q']
        
        blast_results['gene_id'] = blast_results['query'].str.split('_').str[0]
        blast_results['probe'] = blast_results['query'].str.split('_').str[1]
        blast_results['exon_target'] = blast_results['exon_target'].str.split('::').str[0]

        print(blast_results)
        
        for gene_id in blast_results['gene_id'].unique():
            print(gene_id)
            probes_with_match = []
            blast_results_gene = blast_results[blast_results['gene_id'] == gene_id]
            blast_results_gene.reset_index(inplace=True, drop=True)
            for idx in blast_results_gene.index:
                probe = blast_results_gene['probe'][idx]
                gene_id_target = mapping_exon_to_gene[blast_results_gene['exon_target'][idx]]

                # only consider hit in other genes
                if gene_id_target != gene_id:
                    # filter out probes that have more than 50% sequence coverage with target sequence, 80% identity with target sequence, and target sequence covers ligation site +- 5
                    if int(blast_results_gene['alignment_length'][idx]) > (probe_length * coverage) and float(blast_results_gene['percent_identity'][idx]) > percent_identity and int(blast_results_gene['start_q'][idx]) < (probe_length //2 - 4) and int(blast_results_gene['end_q'][idx]) > (probe_length //2 + 5):
                        probes_with_match.append(probe)

            # remove all probes with matches
            #print('Number of probes before blast filter: {}'.format(len(mapping_probe_to_exon)))
            mapping_probe_to_exon = mapping_gene_to_probes_unique[gene_id]
            for probe in set(probes_with_match):
                mapping_probe_to_exon.pop(probe)
            #print('Number of probes after blast filter: {}'.format(len(mapping_probe_to_exon)))
            count_selected_probes += len(mapping_probe_to_exon)
            mapping_gene_to_probes_blastn[gene_id] = mapping_probe_to_exon

    if count_genes % 10 == 0:
        print('Processed genes: {}'.format(count_genes))
        print('Processed probes: {}'.format(count_probes))
        print('Probes passed blastn filter: {}'.format(count_selected_probes))
                

    print('Total number of genes processed: {}, total number of probes processed: {} with {} probes passed hamming distance filter.'.format(count_genes, count_probes, count_selected_probes))    
    pickle.dump(mapping_gene_to_probes_blastn, open(file_mapping_gene_to_probes_blastn,'wb')) 

    # remove blast files
    cmd = 'rm -rf {}'.format(dir_output_blast)
    sp.run(cmd, shell=True)

    return file_mapping_gene_to_probes_blastn


############################################

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




