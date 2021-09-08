############################################
# imports
############################################

import os
import pickle
import logging
import pandas as pd

from gtfparse import read_gtf

from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt


############################################
# probe processing
############################################

def get_mapping_gene_exon(dir_output, file_gene_annotation, get_annotation):

    file_mapping_exon_to_gene = '{}mapping_exon_to_gene.pkl'.format(dir_output)

    if get_annotation and not os.path.isfile(file_mapping_exon_to_gene): 

        gene_annotation = read_gtf(file_gene_annotation)
        exon_annotation = gene_annotation[gene_annotation["feature"] == "exon"]

        if not 'exon_id' in exon_annotation.columns:
            exon_annotation['exon_id'] = exon_annotation['transcript_id'] + '_exon' + exon_annotation['exon_number']
            exon_annotation['source'] = 'Refseq'

        mapping_exon_to_gene = pd.Series(exon_annotation.gene_id.values,index=exon_annotation.exon_id).to_dict()
        pickle.dump(mapping_exon_to_gene, open(file_mapping_exon_to_gene,'wb'))
    
    else:
        mapping_exon_to_gene = pickle.load(open(file_mapping_exon_to_gene,'rb'))

    return mapping_exon_to_gene


############################################

def get_mapping_gene_probes(dir_output, file_exon_sequence, mapping_exon_to_gene, probe_length, GC_content_min, GC_content_max, Tm_min, Tm_max, get_annotation):

    file_mapping_gene_to_probes = '{}mapping_gene_to_probes.pkl'.format(dir_output)

    if get_annotation and not os.path.isfile(file_mapping_gene_to_probes):
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
                                    mapping_gene_to_probes[gene_id][probe]['exon_id'].append(exon_id)
                                else:
                                    mapping_gene_to_probes[gene_id][probe] = {'exon_id': [exon_id], 'GC': gc_content, 'Tm': Tm}
                            else:
                                mapping_gene_to_probes[gene_id] = {probe: {'exon_id': [exon_id], 'GC': gc_content, 'Tm': Tm}}
                        
            #if count_exons % 1000 == 0:
            #    print('Processed exons: {}'.format(count_exons))
            #    print('Probes found in total: {}'.format(count_probes))
            #    break

        pickle.dump(mapping_gene_to_probes, open(file_mapping_gene_to_probes,'wb'))

        logging.info('In total {} exons processed with {} probes.'.format(count_exons, count_probes))

    else:
        mapping_gene_to_probes = pickle.load(open(file_mapping_gene_to_probes,'rb'))

    return mapping_gene_to_probes


############################################

def get_unique_probes(dir_output, mapping_gene_to_probes, mapping_gene_to_probes_unique, get_annotation):

    file_mapping_gene_to_probes_unique = '{}mapping_gene_to_probes_unique.pkl'.format(dir_output)

    if get_annotation and not os.path.isfile(file_mapping_gene_to_probes_unique):

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
               
        pickle.dump(mapping_gene_to_probes_unique, open(file_mapping_gene_to_probes_unique,'wb')) 

        logging.info('Total number of genes processed: {}, total number of probes processed: {} with {} unique probes in total.'.format(count_genes, count_probes, len(probes_unique)))
    
    else:
        mapping_gene_to_probes_unique = pickle.load(open(file_mapping_gene_to_probes_unique,'rb'))

    return mapping_gene_to_probes_unique

    ############################################