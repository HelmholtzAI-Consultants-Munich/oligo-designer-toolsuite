############################################
# imports
############################################

import os
import multiprocessing

from pyfaidx import Fasta

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt

import src.utils as utils


############################################
# functions
############################################

'''
## pool is an alternative to process, hwoever, it did not lead better performance on CPU server

def get_probes_pool(number_batchs, batch_size, exon_annotation, genes, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_probes):
    
    # create index file
    Fasta(file_genome_fasta)

    Tm_parameters['nn_table'] = getattr(mt, Tm_parameters['nn_table'])
    Tm_parameters['tmm_table'] = getattr(mt, Tm_parameters['tmm_table'])
    Tm_parameters['imm_table'] = getattr(mt, Tm_parameters['imm_table'])
    Tm_parameters['de_table'] = getattr(mt, Tm_parameters['de_table'])

    print('Multiprocessing with pool')
    number_cpus = multiprocessing.cpu_count() - 2
    pool = multiprocessing.Pool(number_cpus)
    print('Number of cpus: {}'.format(number_cpus))

    jobs = []
    for batch_id in range(number_batchs): 
        genes_batch = genes[(batch_size * batch_id):(min(batch_size * (batch_id+1), len(genes)+1))]
        
        p = pool.apply_async(get_probes_per_batch, args=(batch_id, genes_batch, exon_annotation, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_probes, ))
        jobs.append(p)

    utils.processes_wait(jobs)
    pool.close()
    
    os.remove('{}.fai'.format(file_genome_fasta))
'''

############################################

def get_probes(number_batchs, batch_size, exon_annotation, genes, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_probes):
    
    # create index file
    Fasta(file_genome_fasta)

    Tm_parameters['nn_table'] = getattr(mt, Tm_parameters['nn_table'])
    Tm_parameters['tmm_table'] = getattr(mt, Tm_parameters['tmm_table'])
    Tm_parameters['imm_table'] = getattr(mt, Tm_parameters['imm_table'])
    Tm_parameters['de_table'] = getattr(mt, Tm_parameters['de_table'])

    jobs = []
    for batch_id in range(number_batchs): 
        genes_batch = genes[(batch_size * batch_id):(min(batch_size * (batch_id+1), len(genes)+1))]
    
        proc = multiprocessing.Process(target=get_probes_per_batch, args=(batch_id, genes_batch, exon_annotation, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_probes, ))
        jobs.append(proc)
        proc.start()

    for job in jobs:
        job.join()

    os.remove('{}.fai'.format(file_genome_fasta))


############################################

def get_probes_per_batch(batch_id, genes_batch, exon_annotation, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_probes):
    file_exon_gtf_batch = os.path.join(dir_output_probes, 'exons_batch{}.gtf'.format(batch_id))
    file_exon_fasta_batch = os.path.join(dir_output_probes, 'exons_batch{}.fna'.format(batch_id))
    file_probe_fasta_batch = os.path.join(dir_output_probes, 'probes_batch{}.fna'.format(batch_id))
    
    _get_exome_fasta(genes_batch, exon_annotation, file_exon_gtf_batch, file_exon_fasta_batch, file_genome_fasta)
    gene_probes = _GC_Tm_filter(file_exon_fasta_batch, genes_batch, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max)
    _get_probes_fasta(gene_probes, file_probe_fasta_batch)


############################################

def _get_exome_fasta(genes_batch, exon_annotation, file_exon_gtf_batch, file_exon_fasta_batch, file_genome_fasta):

    exon_annotation_genes = exon_annotation.loc[exon_annotation['gene_id'].isin(genes_batch)].copy()
    exon_annotation_genes['exon'] = exon_annotation_genes['gene_id'] + '_tid' + exon_annotation_genes['transcript_id'] + '_eid' + exon_annotation_genes['exon_id']
    exon_annotation_genes = exon_annotation_genes.sort_values(by=['gene_id', 'transcript_id', 'exon_id'])
    exon_annotation_genes[['seqname','source','exon','start','end','score','strand','frame']].to_csv(file_exon_gtf_batch, sep='\t', header=False, index = False)
    
    # get sequence for exons
    utils.get_fasta(file_exon_gtf_batch, file_genome_fasta, file_exon_fasta_batch)


############################################

def _GC_Tm_filter(file_exon_fasta_batch, genes_batch, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max):

    gene_probes = {key: {} for key in genes_batch}
    gene_total_probes = {key: 0 for key in genes_batch}

    for exon in SeqIO.parse(file_exon_fasta_batch, "fasta"):

        identifier = exon.id.split('::')[0]
        gene_id = identifier.split('_tid')[0]
        transcript_id = identifier.split('_tid')[1].split('_eid')[0]
        exon_id = identifier.split('_eid')[1]

        coordinates = exon.id.split('::')[1]
        chrom = coordinates.split(':')[0]
        start = int(coordinates.split(':')[1].split('-')[0])
        strand = coordinates.split('(')[1].split(')')[0]

        sequence = exon.seq

        if len(sequence) > probe_length:
            total_probes = gene_total_probes[gene_id]
            number_probes = len(sequence)-(probe_length-1)
            probes_sequence = [sequence[i:i+probe_length] for i in range(number_probes)]
            probes_start = [start + i for i in range(number_probes)]
            probes_end = [start + i + probe_length for i in range(number_probes)]

            probes_id = [total_probes + i for i in range(number_probes)]
            gene_total_probes[gene_id] = total_probes + number_probes

            for i in range(len(probes_sequence)):
                probe_sequence = probes_sequence[i]

                if 'N' not in probe_sequence:
                    gc_content = round(GC(probe_sequence),2)

                    if (GC_content_min < gc_content < GC_content_max):
                        Tm = round(mt.Tm_NN(probe_sequence, **Tm_parameters),2)

                        if (Tm_min < Tm < Tm_max):
                            probe_id = probes_id[i]
                            probe_start = probes_start[i]
                            probe_end = probes_end[i]

                            tmp = gene_probes[gene_id]
                            if probe_sequence in tmp:
                                tmp[probe_sequence]['transcript_id'].append(transcript_id)
                                tmp[probe_sequence]['exon_id'].append(exon_id)
                                tmp[probe_sequence]['probe_id'].append(probe_id)
                                tmp[probe_sequence]['start'].append(probe_start)
                                tmp[probe_sequence]['end'].append(probe_end)
                            else:
                                tmp[probe_sequence] = {'transcript_id': [transcript_id], 'exon_id': [exon_id], 'probe_id': [probe_id], 'chr': chrom, 'start': [probe_start], 'end': [probe_end], 'strand': strand, 'gc': gc_content, 'Tm': Tm}
                            gene_probes[gene_id] = tmp
                            
    #print('Number of possible probes in batch {}: {}'.format(file_exon_fasta_batch.split('.fna')[0].split('_batch')[1], sum(gene_total_probes.values())))
    #print('Number of selected probes in batch {}: {}'.format(file_exon_fasta_batch.split('.fna')[0].split('_batch')[1], sum({key: len(value) for key, value in gene_probes.items()}.values())))
    return gene_probes

############################################

def _get_probes_fasta(gene_probes, file_probe_fasta_batch):
    output = []
    for gene_id, probes in gene_probes.items():
        for probe_sequence, probe_attributes in probes.items():
            header = '{}_tid{}_eid{}_pid{}_seq{}_chr{}_start{}_end{}_strand{}_gc{}_tm{}'.format(gene_id, ';'.join(probe_attributes['transcript_id']), 
                        ';'.join(probe_attributes['exon_id']), ';'.join(str(p) for p in probe_attributes['probe_id']), probe_sequence, probe_attributes['chr'], 
                        ';'.join(str(s) for s in probe_attributes['start']), ';'.join(str(e) for e in probe_attributes['end']), probe_attributes['strand'], probe_attributes['gc'], probe_attributes['Tm'])
            sequence = SeqRecord(probe_sequence, header, '', '')
            output.append(sequence)
         
    with open(file_probe_fasta_batch, 'w') as handle:
        SeqIO.write(output, handle, 'fasta')
        








