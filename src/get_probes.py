############################################
# imports
############################################

import os
import multiprocessing

from pyfaidx import Fasta

from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt

import src.utils as utils


############################################
# functions
############################################

def get_probes(number_batchs, batch_size, exon_annotation, genes, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_annotations):
    """
    Get the fasta sequence of all possible probes with user-defined length for all input genes. 
    Generated probes are filtered by GC content, melting temperature (Tm) and undefined nucleotides ('N') in their sequence. 
    This process can be executed in a parallele fashion on a user-defined number of threads. 
    Parameters
    ----------
        number_batchs: int
            Number of threads for multiprocessing.
        batch_size: int
            Number of genes processed in one batch.
        exon_annotation: pandas.DataFrame
            Dataframe with exon annotation.
        genes: list
            List of genes for which probes should be designed.
        probe_length: int
            Length of designed probe.
        GC_content_min: float
            Minimal GC content of a probe.
        GC_content_max: float
            Maximum GC content of a probe.
        Tm_parameters: dict
            Dictionary with parameters for 'MeltingTemp' fucntion that computed the melting temperature for a sequence.
        Tm_min: float
            Minimal melting temperature of a probe.
        Tm_max: float
            Maximum melting temperature of a probe.
        file_genome_fasta: string
            Path to fasta file with genome sequence.
        dir_output_annotations: string
            Path to directory for annotation files.  
    Returns
    -------
        --- none ---
    """
    # create index file
    Fasta(file_genome_fasta)

    jobs = []
    for batch_id in range(number_batchs): 
        genes_batch = genes[(batch_size * batch_id):(min(batch_size * (batch_id+1), len(genes)+1))]
    
        proc = multiprocessing.Process(target=get_probes_per_batch, args=(batch_id, genes_batch, exon_annotation, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_annotations, ))
        jobs.append(proc)
        proc.start()

    print(jobs)

    for job in jobs:
        job.join()

    os.remove('{}.fai'.format(file_genome_fasta))


############################################

def get_probes_per_batch(batch_id, genes_batch, exon_annotation, probe_length, GC_content_min, GC_content_max, Tm_parameters, Tm_min, Tm_max, file_genome_fasta, dir_output_annotations):
    """
    Get the fasta sequence of all possible probes for all genes in the batch.
    Parameters
    ----------
        batch_id: int
            Batch ID.
        genes_batch: list
            List of genes for which probes should be designed.
        exon_annotation: pandas.DataFrame
            Dataframe with exon annotation.
        probe_length: int
            Length of designed probe.
        GC_content_min: float
            Minimal GC content of a probe.
        GC_content_max: float
            Maximum GC content of a probe.
        Tm_parameters: dict
            Dictionary with parameters for 'MeltingTemp' fucntion that computed the melting temperature for a sequence.
        Tm_min: float
            Minimal melting temperature of a probe.
        Tm_max: float
            Maximum melting temperature of a probe.
        file_genome_fasta: string
            Path to fasta file with genome sequence.
        dir_output_annotations: string
            Path to directory for annotation files.   
    Returns
    -------
        --- none ---  
    """
    file_exon_gtf_batch = os.path.join(dir_output_annotations, 'exons_batch{}.gtf'.format(batch_id))
    file_exon_fasta_batch = os.path.join(dir_output_annotations, 'exons_batch{}.fna'.format(batch_id))
    file_probe_info_batch = os.path.join(dir_output_annotations, 'probes_info_batch{}.txt'.format(batch_id))
    file_probe_sequence_batch = os.path.join(dir_output_annotations, 'probes_sequence_batch{}.txt'.format(batch_id))
    
    _get_exome_fasta(genes_batch, exon_annotation, file_exon_gtf_batch, file_exon_fasta_batch, file_genome_fasta)
    gene_probes = _get_probes_info(file_exon_fasta_batch, genes_batch, probe_length, Tm_parameters)
    _write_probes_info(gene_probes, file_probe_info_batch, file_probe_sequence_batch)


############################################

def _get_exome_fasta(genes_batch, exon_annotation, file_exon_gtf_batch, file_exon_fasta_batch, file_genome_fasta):
    """
    Extract exons for the current batch and write exon regions to gtf file. 
    Get sequence for annotated exon regions (gtf file) from genome sequence (fasta file) and write exon sequences to fasta file.
    Parameters
    ----------
        genes_batch: list
            List of genes for which probes should be designed.
        exon_annotation: pandas.DataFrame
            Dataframe with exon annotation.
        file_exon_gtf_batch: string
            Path to gtf exon annotation output file.
        file_exon_fasta_batch: string
            Path to fasta exon sequence output file.
        file_genome_fasta: string
            Path to fasta file with genome sequence.
    Returns
    -------
        --- none ---       
    """
    exon_annotation_genes = exon_annotation.loc[exon_annotation['gene_id'].isin(genes_batch)].copy()
    exon_annotation_genes['exon'] = exon_annotation_genes['gene_id'] + '_tid' + exon_annotation_genes['transcript_id'] + '_eid' + exon_annotation_genes['exon_id']
    exon_annotation_genes = exon_annotation_genes.sort_values(by=['gene_id', 'transcript_id', 'exon_id'])
    exon_annotation_genes[['seqname','source','exon','start','end','score','strand','frame']].to_csv(file_exon_gtf_batch, sep='\t', header=False, index = False)
    
    # get sequence for exons
    utils.get_fasta(file_exon_gtf_batch, file_genome_fasta, file_exon_fasta_batch)
        

############################################

def _get_probes_info(file_exon_fasta_batch, genes_batch, probe_length, Tm_parameters):
    gene_probes = {key: {} for key in genes_batch}

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
            number_probes = len(sequence)-(probe_length-1)
            probes_sequence = [sequence[i:i+probe_length] for i in range(number_probes)]
            probes_start = [start + i for i in range(number_probes)]
            probes_end = [start + i + probe_length for i in range(number_probes)]

            for i in range(len(probes_sequence)):
                probe_sequence = probes_sequence[i]
                
                if 'N' not in probe_sequence:

                    tmp = gene_probes[gene_id]

                    probe_start = probes_start[i]
                    probe_end = probes_end[i]

                    if probe_sequence in tmp:
                        tmp[probe_sequence]['transcript_id'].append(transcript_id)
                        tmp[probe_sequence]['exon_id'].append(exon_id)
                        tmp[probe_sequence]['start'].append(probe_start)
                        tmp[probe_sequence]['end'].append(probe_end)
                    else:
                        gc_content = round(GC(probe_sequence),2)
                        Tm = round(mt.Tm_NN(probe_sequence, **Tm_parameters),2)
                        tmp[probe_sequence] = {'transcript_id': [transcript_id], 'exon_id': [exon_id], 'chr': chrom, 'start': [probe_start], 'end': [probe_end], 'strand': strand, 'gc': gc_content, 'Tm': Tm}
                    gene_probes[gene_id] = tmp
                            
    return gene_probes


############################################

def _write_probes_info(gene_probes, file_probe_info_batch, file_probe_sequence_batch):

    with open(file_probe_info_batch, 'w') as handle1:
        handle1.write('gene_id\ttranscript_id\texon_id\tprobe_sequence\tchromosome\tstart\tend\tstrand\tGC_content\tmelting_temperature\n')

        with open(file_probe_sequence_batch, 'w') as handle2:
            
            for gene_id, probes in gene_probes.items():
                for probe_sequence, probe_attributes in probes.items():

                    output1 = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_id, ';'.join(probe_attributes['transcript_id']), ';'.join(probe_attributes['exon_id']), 
                            probe_sequence, probe_attributes['chr'], ';'.join(str(s) for s in probe_attributes['start']), ';'.join(str(e) for e in probe_attributes['end']), probe_attributes['strand'], 
                            probe_attributes['gc'], probe_attributes['Tm'])
                    handle1.write(output1)

                    output2 = '{}\n'.format(probe_sequence)
                    handle2.write(output2)