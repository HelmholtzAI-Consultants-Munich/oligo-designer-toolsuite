############################################
# imports
############################################

import os
import pandas as pd
import multiprocessing

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

import src.utils as utils

############################################
# functions 
###########################################

def run_blast_search(number_batchs, word_size, num_threads_blast, file_gene_gtf, file_genome_fasta, dir_output_annotations, dir_output_blast):   
    """
    Run BlastN alignment tool to find regions of local similarity between sequences, where sequences are probes and transcripts.
    BlastN identifies the transcript regions where probes match with a certain coverage and similarity.
    Parameters
    ----------
        number_batchs: int
            Number of threads for multiprocessing.
        word_size: int
            Word size for the blastn seed (exact match to target).
        num_threads_blast: int
            Number of threads for blastN.
        file_gene_gtf: string
            Path to gtf file with gene annotation.
        file_genome_fasta: string
            Path to fasta file with genome sequence.
        dir_output_annotations: string
            Path to directory for annotation files. 
        dir_output_blast: string
            Path to output directory for blastN aligment search results.
    Returns
    -------
        --- none ---
    """
    file_transcriptome_gtf = _get_transcriptome_gtf(file_gene_gtf, dir_output_annotations)
    file_transcriptome_fasta = os.path.join(dir_output_annotations, 'transcriptome.fna')
    utils.get_fasta(file_transcriptome_gtf, file_genome_fasta, file_transcriptome_fasta)

    # create blast database
    cmd = NcbimakeblastdbCommandline(input_file=file_transcriptome_fasta, dbtype='nucl')
    out, err = cmd()
    print('Blast database created.')

    # run blast with multi process
    jobs = []
    for batch_id in range(number_batchs):
        proc = multiprocessing.Process(target=_run_blast, args=(batch_id, word_size, num_threads_blast, file_transcriptome_fasta, dir_output_annotations, dir_output_blast, ))
        jobs.append(proc)
        proc.start()

    print(jobs)

    for job in jobs:
        job.join()   


############################################

def filter_probes_by_blast_results(number_batchs, probe_length, percent_identity, coverage, min_probes_per_gene, dir_output_blast, dir_output_results):
    """
    Process the results from BlastN alignment search and filter probes based on the results. 
    Parameters
    ----------
        number_batchs: int
            Number of threads for multiprocessing.
        probe_length: int
            Length of designed probe.
        percent_identity: float
            Maximum similarity between probes and target sequences (blast parameter "percent identity"), ranging from 0 tp 1 (no missmatch).
        coverage: float
            Minimum coverage between probes and target sequence (blast parameter "coverage"), ranging from 0 tp 1 (full coverage).
        min_probes_per_gene
            Minimum number of probes per gene.
        dir_output_blast: string
            Path to output directory for blastN aligment search results.
        dir_output_probes: string
            Path to output directory for probe design results. 
    Returns
    -------
        --- none ---
    """   
    file_removed_genes = os.path.join(dir_output_results, 'genes_with_insufficient_probes.txt')

    # process blast results
    jobs = []
    for batch_id in range(number_batchs):
        proc = multiprocessing.Process(target=_process_blast_results, args=(batch_id, probe_length, percent_identity, coverage, min_probes_per_gene, file_removed_genes, dir_output_blast, dir_output_results, ))
        jobs.append(proc)
        proc.start()

    for job in jobs:
        job.join()


############################################

def _get_transcriptome_gtf(file_gene_gtf, dir_output_annotations):
    """
    Get annotation (gtf file) of transcripts. 
    Parameters
    ----------
        file_gene_gtf: string
            Path to gtf file with gene annotation.
        dir_output_annotations: string
            Path to directory for annotation files. 
    Returns
    -------
        file_transcriptome_gtf: string
            Path to gtf file with transcriptome annotation.
    """
    file_transcriptome_gtf = os.path.join(dir_output_annotations, 'transcriptome.gtf')

    transcriptome_annotation = utils.load_transcriptome_annotation(file_gene_gtf)
    transcriptome_annotation['transcript'] = transcriptome_annotation['gene_id'] + '_tid' + transcriptome_annotation['transcript_id']
    transcriptome_annotation[['seqname','source','transcript','start','end','score','strand','frame']].to_csv(file_transcriptome_gtf, sep='\t', header=False, index = False)
   
    return file_transcriptome_gtf


############################################

def _run_blast(batch_id, word_size, num_threads_blast, file_transcriptome_fasta, dir_output_annotations, dir_output_blast):
    """
    Run BlastN alignment search for all probes of one batch. 
    Parameters
    ----------
        batch_id: int
            Batch ID.
        word_size: int
            Word size for the blastn seed (exact match to target).
        num_threads_blast: int
            Number of threads for blastN.
        file_transcriptome_fasta: tring
            Path to fasta file with transcriptome sequences.
        dir_output_annotations: string
            Path to directory for annotation files. 
        dir_output_blast: string
            Path to output directory for blastN aligment search results.
    Returns
    -------
        --- none ---
    """
    file_probe_fasta_batch = os.path.join(dir_output_annotations, 'probes_batch{}.fna'.format(batch_id))
    file_blast_batch = os.path.join(dir_output_blast, 'blast_batch{}.txt'.format(batch_id))

    cmd = NcbiblastnCommandline(query=file_probe_fasta_batch,db=file_transcriptome_fasta, outfmt="10 qseqid sseqid pident length qstart qend", out=file_blast_batch, word_size=word_size, strand='plus', num_threads=num_threads_blast)
    out, err = cmd()


############################################

def _process_blast_results(batch_id, probe_length, percent_identity, coverage, min_probes_per_gene, file_removed_genes, dir_output_blast, dir_output_probes):
    """
    Process the output of the BlastN alignment search. 
    Parameters
    ----------
        batch_id: int
            Batch ID.
        probe_length: int
            Length of designed probe.
        percent_identity: float
            Maximum similarity between probes and target sequences (blast parameter "percent identity"), ranging from 0 tp 1 (no missmatch).
        coverage: float
            Minimum coverage between probes and target sequence (blast parameter "coverage"), ranging from 0 tp 1 (full coverage).
        min_probes_per_gene
            Minimum number of probes per gene.
        file_removed_genes: string
            Path to text file with genes that were removed due to low number of probes.
        dir_output_blast: string
            Path to output directory for blastN aligment search results.
        dir_output_probes: string
            Path to output directory for probe design results. 
    Returns
    -------
        --- none ---
    """
    blast_results = _read_blast_output(batch_id, dir_output_blast)
    _filter_probes(blast_results, probe_length, percent_identity, coverage, min_probes_per_gene, file_removed_genes, dir_output_probes)


############################################

def _read_blast_output(batch_id, dir_output_blast):
    """
    Load the output of the BlastN alignment search into a DataFrame and process the results.
    Parameters
    ----------
        batch_id: int
            Batch ID.
        dir_output_blast: string
            Path to output directory for blastN aligment search results.
    Returns
    -------
        blast_results: pandas.DataFrame
            DataFrame with processed blast alignment search results.  
    """
    file_blast_batch = os.path.join(dir_output_blast, 'blast_batch{}.txt'.format(batch_id))

    blast_results = pd.read_csv(file_blast_batch, header=None)
    blast_results.columns = ['query','target_exon','percent_identity','alignment_length','query_start','query_end']

    blast_results['query_gene_id'] = blast_results['query'].str.split('_tid').str[0]
    blast_results['target_exon'] = blast_results['target_exon'].str.split('::').str[0]

    return blast_results


############################################

def _filter_probes(blast_results, probe_length, percent_identity, coverage, min_probes_per_gene, file_removed_genes, dir_output_probes):
    """
    Use the results of the BlastN alignement search to remove probes with high similarity and coverage based on user defined thresholds.
    Parameters
    ----------
        blast_results: pandas.DataFrame
            DataFrame with processed blast alignment search results.
        probe_length: int
            Length of designed probe.
        percent_identity: float
            Maximum similarity between probes and target sequences (blast parameter "percent identity"), ranging from 0 tp 1 (no missmatch).
        coverage: float
            Minimum coverage between probes and target sequence (blast parameter "coverage"), ranging from 0 tp 1 (full coverage).
        min_probes_per_gene
            Minimum number of probes per gene.
        file_removed_genes: string
            Path to text file with genes that were removed due to low number of probes.
        dir_output_probes: string
            Path to output directory for probe design results. 
    Returns
    -------
        --- none ---
    """
    removed_genes = [] 

    for gene_id in blast_results['query_gene_id'].unique():
        probes_wo_match = _get_probes_wo_match(gene_id, blast_results, probe_length, percent_identity, coverage)

        if len(probes_wo_match) < min_probes_per_gene:
            removed_genes.append(gene_id)
        else:
            _write_output(gene_id, probes_wo_match, dir_output_probes)

    with open(file_removed_genes, 'a') as output:
        for gene_id in removed_genes:
            output.write('{}\n'.format(gene_id))


############################################

def _get_probes_wo_match(gene_id, blast_results, probe_length, percent_identity, coverage):
    """
    Get a list of suitable probes for a gene. 
    Parameters
    ----------
        gene_id: string
            Gene ID of processed gene.
        blast_results: pandas.DataFrame
            DataFrame with processed blast alignment search results.
        probe_length: int
            Length of designed probe.
        percent_identity: float
            Maximum similarity between probes and target sequences (blast parameter "percent identity"), ranging from 0 tp 1 (no missmatch).
        coverage: float
            Minimum coverage between probes and target sequence (blast parameter "coverage"), ranging from 0 tp 1 (full coverage).
    Returns
    -------
        probes_wo_match: list
            List of suitable probes that don't have matches in the transcriptome.
    """
    blast_results_gene = blast_results[blast_results['query_gene_id'] == gene_id]
    blast_results_gene.reset_index(inplace=True, drop=True)

    probes_wo_match = {}
    probes_with_match = []

    for idx in blast_results_gene.index:
        header = blast_results_gene['query'][idx]
        target_exon = blast_results_gene['target_exon'][idx]
        gene_id_target = target_exon.split('_tid')[0]

        if header not in probes_with_match:
            probes_wo_match[header] = ''
            if gene_id_target != gene_id:
                if float(blast_results_gene['percent_identity'][idx]) > percent_identity and int(blast_results_gene['alignment_length'][idx]) > (probe_length * coverage) and int(blast_results_gene['query_start'][idx]) < (probe_length //2 - 4) and int(blast_results_gene['query_end'][idx]) > (probe_length //2 + 5):
                    probes_with_match.append(header)
                    probes_wo_match.pop(header)
    
    probes_wo_match = list(set(probes_wo_match.keys()))

    return probes_wo_match


############################################

def _write_output(gene_id, probes_wo_match, dir_output_probes):
    """
    Write results of probe design pipeline to file and create one file with suitable probes per gene. 
    Parameters
    ----------
        gene_id: string
            Gene ID of processed gene.
        probes_wo_match: list
            List of suitable probes that don't have matches in the transcriptome.
        dir_output_probes: string
            Path to output directory for probe design results. 
    Returns
    -------
        --- none ---    
    """
    file_output = '{}probes_{}'.format(dir_output_probes, gene_id)

    with open(file_output, 'w') as handle:
        handle.write('probe_sequence\tchromosome\tstart\tend\tstrand\tgene_id\ttranscript_id\texon_id\tprobe_id\tGC_content\tmelting_temperature\n')

        for header in probes_wo_match:
            transcript_id = header.split('_tid')[1].split('_eid')[0]
            exon_id = header.split('_eid')[1].split('_pid')[0]
            probe_id = header.split('_pid')[1].split('_seq')[0]
            probe_seq = header.split('_seq')[1].split('_chr')[0]
            chrom = header.split('_chr')[1].split('_start')[0]
            start = header.split('_start')[1].split('_end')[0]
            end = header.split('_end')[1].split('_strand')[0]
            strand = header.split('_strand')[1].split('_gc')[0]
            gc = header.split('_gc')[1].split('_tm')[0]
            Tm = header.split('_tm')[1]

            if len(set(start.split(';'))) > 1:
                frame_multi_s = pd.DataFrame({'transcript_id':transcript_id.split(';'), 'exon_id':exon_id.split(';'), 'probe_id':probe_id.split(';'), 'start': start.split(';'), 'end':end.split(';')})
                for s in set(start.split(';')):
                    frame_s = frame_multi_s.loc[frame_multi_s['start'] == s]

                    output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(probe_seq, chrom, frame_s.start.unique()[0], frame_s.end.unique()[0], strand, gene_id, ';'.join(frame_s.transcript_id), ';'.join(frame_s.exon_id), ';'.join(frame_s.probe_id), gc, Tm)
                    handle.write(output)
            else:
                output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(probe_seq, chrom, start.split(';')[0], end.split(';')[0], strand, gene_id, transcript_id, exon_id, probe_id, gc, Tm)
                handle.write(output)








        

