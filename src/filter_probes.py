############################################
# imports
############################################

import os
import pandas as pd
import multiprocessing

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

import src.utils as utils

############################################
# functions 
###########################################

def filter_probes_by_GC_Tm_exactmatch(number_batchs, GC_content_min, GC_content_max, Tm_min, Tm_max, dir_output_annotations): 
    duplicated_sequences = _get_duplicated_sequences(number_batchs, dir_output_annotations)
    
    # run filter with multiprocess
    jobs = []
    for batch_id in range(number_batchs):
        proc = multiprocessing.Process(target=_filter_probes_GC_Tm_exactmatch, args=(batch_id, duplicated_sequences, GC_content_min, GC_content_max, Tm_min, Tm_max, dir_output_annotations, ))
        jobs.append(proc)
        proc.start()

    print(jobs)

    for job in jobs:
        job.join()  


############################################

def _get_duplicated_sequences(number_batchs, dir_output_annotations):
    sequences = []
    
    for batch_id in range(number_batchs):
        file_probe_sequence_batch = os.path.join(dir_output_annotations, 'probes_sequence_batch{}.txt'.format(batch_id))

        with open(file_probe_sequence_batch, "r") as handle:
            sequences_batch = [line.rstrip() for line in handle]
        sequences.extend(sequences_batch)
        os.remove(file_probe_sequence_batch)
    
    seen_sequences = []
    duplicated_sequences = []
    
    for sequence in sequences:
        if sequence not in seen_sequences:
            seen_sequences.append(sequence)
        else:
            duplicated_sequences.append(sequence)
    
    duplicated_sequences = list(set(duplicated_sequences))
    print(duplicated_sequences)
    return duplicated_sequences


############################################

def _filter_probes_GC_Tm_exactmatch(batch_id, duplicated_sequences, GC_content_min, GC_content_max, Tm_min, Tm_max, dir_output_annotations):
    file_probe_info_batch = os.path.join(dir_output_annotations, 'probes_info_batch{}.txt'.format(batch_id))
    file_probe_fasta_batch = os.path.join(dir_output_annotations, 'probes_batch{}.fna'.format(batch_id))

    probes_info = pd.read_csv(file_probe_info_batch, sep='\t')
    os.remove(file_probe_info_batch)

    probes_info_filtered = []

    for gene in probes_info.gene_id.unique():
        probe_info_gene = probes_info.loc[probes_info.gene_id == gene]
        probe_info_gene.reset_index(inplace=True, drop=True)
        
        for row in probe_info_gene.index:
            sequence = probe_info_gene.iloc[row, probe_info_gene.columns.get_loc('probe_sequence')]
            
            if sequence not in duplicated_sequences:
                gc_content = probe_info_gene.iloc[row, probe_info_gene.columns.get_loc('GC_content')]
                
                if (GC_content_min < gc_content < GC_content_max):
                    Tm = probe_info_gene.iloc[row, probe_info_gene.columns.get_loc('melting_temperature')]
                    
                    if (Tm_min < Tm < Tm_max):
                        probe = probe_info_gene.iloc[row].to_dict()
                        probe['probe_id'] = '{}_pid{}'.format(gene, row)
                        probes_info_filtered.append(probe)

            else:
                print('Duplicated probe: {}'.format(sequence))
    
    probes_info_filtered = pd.DataFrame(probes_info_filtered)
    _write_probes(probes_info_filtered, file_probe_info_batch, file_probe_fasta_batch)


############################################

def _write_probes(probes_info_filtered, file_probe_info_batch, file_probe_fasta_batch):
    # save info table
    probes_info_filtered[['probe_id', 'probe_sequence', 'gene_id', 'transcript_id', 'exon_id', 'chromosome', 'start', 'end', 'strand', 'GC_content', 'melting_temperature']].to_csv(file_probe_info_batch, sep='\t', index=False)

    # save sequence of probes in fasta format
    output = []
    for row in probes_info_filtered.index:
        header = probes_info_filtered.iloc[row, probes_info_filtered.columns.get_loc('probe_id')]
        sequence = probes_info_filtered.iloc[row, probes_info_filtered.columns.get_loc('probe_sequence')]
        output.append(SeqRecord(Seq(sequence), header, '', ''))

    with open(file_probe_fasta_batch, 'w') as handle:
        SeqIO.write(output, handle, 'fasta')
         

############################################

def run_blast_search(number_batchs, word_size, percent_identity, num_threads_blast, file_gene_gtf, file_genome_fasta, dir_output_annotations, dir_output_blast):   
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
        proc = multiprocessing.Process(target=_run_blast, args=(batch_id, word_size, percent_identity, num_threads_blast, file_transcriptome_fasta, dir_output_annotations, dir_output_blast, ))
        jobs.append(proc)
        proc.start()

    print(jobs)

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
    transcriptome_annotation[['seqname','source','gene_id','start','end','score','strand','frame']].to_csv(file_transcriptome_gtf, sep='\t', header=False, index = False)
   
    return file_transcriptome_gtf


############################################

def _run_blast(batch_id, word_size, percent_identity, num_threads_blast, file_transcriptome_fasta, dir_output_annotations, dir_output_blast):
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

    cmd = NcbiblastnCommandline(query=file_probe_fasta_batch,db=file_transcriptome_fasta, outfmt="10 qseqid sseqid length qstart qend", out=file_blast_batch, word_size=word_size, perc_identity=percent_identity, num_threads=num_threads_blast)
    out, err = cmd()


############################################

def filter_probes_by_blast_results(number_batchs, probe_length, coverage, ligation_site, min_probes_per_gene, dir_output_annotations, dir_output_blast, dir_output_results):
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
        proc = multiprocessing.Process(target=_process_blast_results, args=(batch_id, probe_length, coverage, ligation_site, min_probes_per_gene, file_removed_genes, dir_output_annotations, dir_output_blast, dir_output_results, ))
        jobs.append(proc)
        proc.start()

    for job in jobs:
        job.join()


############################################

def _process_blast_results(batch_id, probe_length, coverage, ligation_site, min_probes_per_gene, file_removed_genes, dir_output_annotations, dir_output_blast, dir_output_probes):
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
        dir_output_annotations: string
            Path to directory for annotation files. 
        dir_output_blast: string
            Path to output directory for blastN aligment search results.
        dir_output_probes: string
            Path to output directory for probe design results. 
    Returns
    -------
        --- none ---
    """
    probes_info = _load_probes_info(batch_id, dir_output_annotations)
    blast_results = _read_blast_output(batch_id, dir_output_blast)
    _filter_probes_blast(probes_info, blast_results, probe_length, coverage, ligation_site, min_probes_per_gene, file_removed_genes, dir_output_probes)


############################################

def _load_probes_info(batch_id, dir_output_annotations):
    file_probe_info_batch = os.path.join(dir_output_annotations, 'probes_info_batch{}.txt'.format(batch_id))
    probes_info = pd.read_csv(file_probe_info_batch, sep='\t')
    return probes_info


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
    blast_results.columns = ['query','target','alignment_length','query_start','query_end']

    blast_results['query_gene_id'] = blast_results['query'].str.split('_pid').str[0]
    blast_results['target_gene_id'] = blast_results['target'].str.split('::').str[0]

    return blast_results


############################################

def _filter_probes_blast(probes_info, blast_results, probe_length, coverage, ligation_site, min_probes_per_gene, file_removed_genes, dir_output_probes):
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
        probes_wo_match = _get_probes_wo_match(gene_id, blast_results, probe_length, coverage, ligation_site)

        if len(probes_wo_match) < min_probes_per_gene:
            removed_genes.append(gene_id)
        else:
            _write_output(probes_info, gene_id, probes_wo_match, dir_output_probes)

    with open(file_removed_genes, 'a') as output:
        for gene_id in removed_genes:
            output.write('{}\n'.format(gene_id))


############################################

def _get_probes_wo_match(gene_id, blast_results, probe_length, coverage, ligation_site):
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
        gene_id_target = blast_results_gene['target_gene_id'][idx]

        if header not in probes_with_match:
            probes_wo_match[header] = ''
            if gene_id_target != gene_id:
                if ligation_site > 0:
                    if int(blast_results_gene['alignment_length'][idx]) > (probe_length * coverage) and int(blast_results_gene['query_start'][idx]) < (probe_length //2 - (ligation_site-1)) and int(blast_results_gene['query_end'][idx]) > (probe_length //2 + ligation_site):
                        probes_with_match.append(header)
                        probes_wo_match.pop(header)
                else:
                    if int(blast_results_gene['alignment_length'][idx]) > (probe_length * coverage):
                        probes_with_match.append(header)
                        probes_wo_match.pop(header)
    
    probes_wo_match = list(set(probes_wo_match.keys()))

    return probes_wo_match


############################################

def _write_output(probes_info, gene_id, probes_wo_match, dir_output_probes):
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
    valid_probes = probes_info[probes_info['probe_id'].isin(probes_wo_match)]
    file_output = '{}probes_{}'.format(dir_output_probes, gene_id)

    valid_probes.to_csv(file_output, sep='\t', index=False)
    






        

