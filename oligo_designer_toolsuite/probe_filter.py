############################################
# imports
############################################

import os
import re
import shutil
import multiprocessing
import iteration_utilities

import pandas as pd

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

############################################
# probe filter class
############################################


class ProbeFilter:
    '''This class is used to filter the set of all possible guides per gene based on:
    1. exact matches,
    2. results of blast alignment tool.
    For each gene the set of filtered probes is saved and a list of genes where no probes were left after filtering.

    :param config: User-defined parameters, where keys are the parameter names and values are the paremeter values.
    :type config: dict
    :param logging: Logger object to store important information.
    :type logging: logging.basicConfig
    :param dir_output: User-defined output directory.
    :type dir_output: string
    '''

    def __init__(self, config, logging, dir_output, file_transcriptome_fasta, genes, dir_annotations = None):
        """Constructor method
        """
         # set logger
        self.logging = logging

        # set directory
        if dir_annotations == None:
            self.dir_annotations = os.path.join(dir_output, 'annotations')
            Path(self.dir_annotations).mkdir(parents=True, exist_ok=True)
        else: 
            self.dir_annotations = dir_annotations
        self.dir_blast = os.path.join(dir_output, 'blast')
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_probes = os.path.join(dir_output, 'probes')
        Path(self.dir_probes).mkdir(parents=True, exist_ok=True)

        self.file_removed_genes = os.path.join(dir_output, 'genes_with_insufficient_probes.txt')

        # set parameters
        self.number_batchs = config['number_batchs']
        self.number_subbatches = config['number_subbatches']
        self.num_threads_blast = config['num_threads']
        self.word_size = config['word_size'] 
        self.percent_identity = config['percent_identity']
        self.probe_length_min = config['probe_length_min']
        self.probe_length_max = config['probe_length_max']
        self.coverage = config['coverage']
        self.ligation_site = config['ligation_site']
        self.min_probes_per_gene = config['min_probes_per_gene']
        self.file_transcriptome_fasta = file_transcriptome_fasta
        self.removed_genes = genes

        # initialize additional paremeters
        self.duplicated_sequences = None


    def filter_probes_by_exactmatch(self):
        '''Probes with exact matches are filtered out.
        This process can be executed in a parallele fashion on a user-defined number of threads (defined in config file).
        '''

        def _get_duplicated_sequences():
            '''Get a list of probe sequences that have a exact match within the pool of all 
            possible probe sequences for the list of input genes.

            :return: List of probe sequences with exact matches in the pool of probes.
            :rtype: list
            '''
            sequences = []
            
            for batch_id in range(self.number_batchs):
                file_probe_sequence_batch = os.path.join(self.dir_annotations, 'probes_sequence_batch{}.txt'.format(batch_id))

                with open(file_probe_sequence_batch, 'r') as handle:
                    sequences_batch = [line.rstrip() for line in handle]
                sequences.extend(sequences_batch)
                os.remove(file_probe_sequence_batch)
            
            duplicated_sequences = list(iteration_utilities.unique_everseen(iteration_utilities.duplicates(sequences)))
            #self.logging.info('Number of duplicated probe sequences: {}'.format(len(duplicated_sequences)))
            
            return duplicated_sequences

        def _filter_probes_exactmatch(batch_id):
            '''Remove sequences with exact matches within the pool of all possible probe sequences for the list of input genes.

            :param batch_id: Batch ID.
            :type batch_id: int
            '''
            file_probe_info_batch = os.path.join(self.dir_annotations, 'probes_info_batch{}.txt'.format(batch_id))
            file_probe_fasta_batch = os.path.join(self.dir_annotations, 'probes_sequence_batch{}.fna'.format(batch_id))

            probes_info = pd.read_csv(file_probe_info_batch, sep='\t',
                                    dtype={'gene_id': str, 'transcript_id': str, 'exon_id': str, 'probe_sequence': str, 
                                           'chromosome': str, 'start': str, 'end': str, 'strand': str, 
                                           'GC_content': float, 'melting_temperature': float, 'melt_temp_arm1':float,
                                           'melt_temp_arm2': float, 'melt_temp_dif_arms': float, 'length': int, 
                                           'ligation_site': int})            
            os.remove(file_probe_info_batch)

            probes_info_filetred = probes_info[~ probes_info.probe_sequence.isin(self.duplicated_sequences)]
            probes_info_filetred.reset_index(inplace=True, drop=True)
            # probe_ids = ['{}_pid{}'.format(probes_info_filetred.gene_id[probe_id], probe_id) for probe_id in range(len(probes_info_filetred.index))]
            probe_ids = ['{}_pid{}'.format(g_id, i) for i, g_id in enumerate(probes_info_filetred['gene_id'])]
            probes_info_filetred.insert(0, 'probe_id', probe_ids)
            _write_probes(probes_info_filetred, file_probe_info_batch, file_probe_fasta_batch)


        def _write_probes(probes_info_filtered, file_probe_info_batch, file_probe_fasta_batch):
            '''Save filtered probe information in tsv file. Save probe sequences as fasta file.

            :param probes_info_filtered: Dataframe with probe information, filtered based on sequence properties.
            :type probes_info_filtered: pandas.DataFrame
            :param file_probe_info_batch: Path to tsv file with probe infos.
            :type file_probe_info_batch: string
            :param file_probe_fasta_batch: Path to fast file with probe sequences.
            :type file_probe_fasta_batch: string
            '''
            # save info table
            probes_info_filtered[['probe_id', 'probe_sequence', 'gene_id', 'transcript_id', 'exon_id', 
                                  'chromosome', 'start', 'end', 'strand', 
                                  'GC_content', 'melting_temperature', 'melt_temp_arm1', 'melt_temp_arm2', 'melt_temp_dif_arms', 
                                  'length', 'ligation_site']].to_csv(file_probe_info_batch, sep='\t', index=False)

            # save sequence of probes in fasta format
            genes = probes_info_filtered.gene_id.unique()
            subbatch_size = int(len(genes)/self.number_subbatches) + (len(genes) % self.number_subbatches > 0)
            for subbatch_id in range(self.number_subbatches):
                file_probe_fasta_subbatch = file_probe_fasta_batch.replace('.fna', '_{}.fna'.format(subbatch_id))
                
                genes_subbatch = genes[(subbatch_id*subbatch_size):((subbatch_id+1)*subbatch_size)]
                probes_info_filtered_subbatch = probes_info_filtered.loc[probes_info_filtered['gene_id'].isin(genes_subbatch)].copy()
                probes_info_filtered_subbatch.reset_index(inplace=True, drop=True)

                output = []
                for row in probes_info_filtered_subbatch.index:
                    header = probes_info_filtered_subbatch.iloc[row, probes_info_filtered_subbatch.columns.get_loc('probe_id')]
                    sequence = Seq(probes_info_filtered_subbatch.iloc[row, probes_info_filtered_subbatch.columns.get_loc('probe_sequence')])
                    output.append(SeqRecord(sequence, header, '', ''))

                with open(file_probe_fasta_subbatch, 'w') as handle:
                    SeqIO.write(output, handle, 'fasta')

        # get list of exact matches in probes pool
        self.duplicated_sequences = _get_duplicated_sequences()
                
        # run filter with multiprocess
        jobs = []
        for batch_id in range(self.number_batchs):
            proc = multiprocessing.Process(target=_filter_probes_exactmatch, args=(batch_id, ))
            jobs.append(proc)
            proc.start()

        #print('\n {} \n'.format(jobs))

        for job in jobs:
            job.join()  


    def run_blast_search(self):   
        '''Run BlastN alignment tool to find regions of local similarity between sequences, where sequences are probes and transcripts.
        BlastN identifies the transcript regions where probes match with a certain coverage and similarity.
        '''

        def _run_blast(batch_id):
            '''Run BlastN alignment search for all probes of one batch. 

            :param batch_id: Batch ID.
            :type batch_id: int
            '''
            for subbatch_id in range(self.number_subbatches):
                file_probe_fasta_batch = os.path.join(self.dir_annotations, 'probes_sequence_batch{}_{}.fna'.format(batch_id, subbatch_id))
                file_blast_batch = os.path.join(self.dir_blast, 'blast_batch{}_{}.txt'.format(batch_id, subbatch_id))

                cmd = NcbiblastnCommandline(query=file_probe_fasta_batch, db=self.file_transcriptome_fasta, outfmt="10 qseqid sseqid length qstart qend qlen", out=file_blast_batch,
                                            strand='plus', word_size=self.word_size, perc_identity=self.percent_identity, num_threads=self.num_threads_blast) 
                out, err = cmd()

        # create blast database
        cmd = NcbimakeblastdbCommandline(input_file=self.file_transcriptome_fasta, dbtype='nucl')
        out, err = cmd()

        # run blast with multi process
        jobs = []
        for batch_id in range(self.number_batchs):
            proc = multiprocessing.Process(target=_run_blast, args=(batch_id, ))
            jobs.append(proc)
            proc.start()

        #print('\n {} \n'.format(jobs))

        for job in jobs:
            job.join()   

    
    def filter_probes_by_blast_results(self):
        '''Process the results from BlastN alignment search and filter probes based on the results. 
        '''

        def _process_blast_results(batch_id):
            '''Process the output of the BlastN alignment search. 

            :param batch_id: Batch ID.
            :type batch_id: int
            '''
            probes_info = _load_probes_info(batch_id)

            num_probes_wo_match = 0
            for subbatch_id in range(self.number_subbatches):
                blast_results = _read_blast_output(batch_id, subbatch_id)
                num_probes_wo_match += _filter_probes_blast(probes_info, blast_results)


        def _load_probes_info(batch_id):
            '''Load filtered probe infomration from tsv file. 

            :param batch_id: Batch ID.
            :type batch_id: int
            :return: Dataframe with probe information, filtered based on sequence properties.
            :rtype: pandas.DataFrame
            '''
            file_probe_info_batch = os.path.join(self.dir_annotations, 'probes_info_batch{}.txt'.format(batch_id))
            probes_info = pd.read_csv(file_probe_info_batch, sep='\t',  
                                      dtype={'gene_id': str, 'transcript_id': str, 'exon_id': str, 'probe_sequence': str, 
                                             'chromosome': str, 'start': str, 'end': str, 'strand': str, 
                                             'GC_content': float, 'melting_temperature': float, 'melt_temp_arm1':float,
                                             'melt_temp_arm2': float, 'melt_temp_dif_arms': float, 'length': int, 
                                             'ligation_site': int})
            return probes_info

        def _read_blast_output(batch_id, subbatch_id):
            '''Load the output of the BlastN alignment search into a DataFrame and process the results.

            :param batch_id: Batch ID.
            :type batch_id: int
            :return: DataFrame with processed blast alignment search results.  
            :rtype: pandas.DataFrame
            '''
            file_blast_batch = os.path.join(self.dir_blast, 'blast_batch{}_{}.txt'.format(batch_id, subbatch_id))
            blast_results = pd.read_csv(file_blast_batch, header=None, sep=',', low_memory=False,
                                        names=['query','target','alignment_length','query_start','query_end','query_length'], engine='c', 
                                        dtype={'query': str, 'target': str, 'alignment_length': int, 'query_start': int, 'query_end': int, 'query_length': int})

            blast_results['query_gene_id'] = blast_results['query'].str.split('_pid').str[0]
            blast_results['target_gene_id'] = blast_results['target'].str.split('::').str[0]

            return blast_results

        def _filter_probes_blast(probes_info, blast_results):
            '''Use the results of the BlastN alignement search to remove probes with high similarity, 
            probe coverage and ligation site coverage based on user defined thresholds.

            :param probes_info: Dataframe with probe information, filtered based on sequence properties.
            :type probes_info: pandas.DataFrame
            :param blast_results: DataFrame with processed blast alignment search results. 
            :type blast_results: pandas.DataFrame
            '''
            # blast_results_matches = blast_results[~(blast_results[['query_gene_id','target_gene_id']].nunique(axis=1) == 1)]
            blast_results_matches = blast_results[blast_results['query_gene_id'] != blast_results['target_gene_id']]
            blast_results_matches_filtered = []
            
            for probe_length in range(self.probe_length_min, self.probe_length_max + 1):

               min_alignment_length = probe_length * self.coverage / 100
               blast_results_matches_probe_length = blast_results_matches[blast_results_matches.query_length == probe_length]
               blast_results_matches_probe_length = blast_results_matches_probe_length[blast_results_matches_probe_length.alignment_length > min_alignment_length]

               if self.ligation_site > 0:
                   ligation_site_start = probe_length // 2 - (self.ligation_site - 1)
                   ligation_site_end = probe_length // 2 + self.ligation_site
                   blast_results_matches_probe_length = blast_results_matches_probe_length[blast_results_matches_probe_length.query_start < ligation_site_start]
                   blast_results_matches_probe_length = blast_results_matches_probe_length[blast_results_matches_probe_length.query_end > ligation_site_end]
                
               blast_results_matches_filtered.append(blast_results_matches_probe_length) 
               
            blast_results_matches_filtered = pd.concat(blast_results_matches_filtered)

            probes_with_match = blast_results_matches_filtered['query'].unique()
            probes_wo_match = blast_results[~blast_results['query'].isin(probes_with_match)]

            for gene_id in blast_results['query_gene_id'].unique():
               probes_wo_match_gene = probes_wo_match[probes_wo_match.query_gene_id == gene_id]
               probes_wo_match_gene = probes_wo_match_gene['query'].unique()

               if len(probes_wo_match_gene) > 0: #gene has to have at least one probe
                   _write_output(probes_info, gene_id, probes_wo_match_gene)

            return len(probes_wo_match['query'].unique())   

        def _write_output(probes_info, gene_id, probes_wo_match):
            '''Write results of probe design pipeline to file and create one file with suitable probes per gene. 

            :param probes_info: Dataframe with probe information, filtered based on sequence properties.
            :type probes_info: pandas.DataFrame
            :param gene_id: Gene ID of processed gene.
            :type gene_id: string
            :param probes_wo_match: List of suitable probes that don't have matches in the transcriptome.
            :type probes_wo_match: list
            '''
            file_output = os.path.join(self.dir_probes, 'probes_{}.txt'.format(gene_id))
            valid_probes = probes_info[probes_info['probe_id'].isin(probes_wo_match)]
            valid_probes.to_csv(file_output, sep='\t', index=False)

        def _write_removed_genes():
            '''Write list of genes for which not enough probes could be designed for.
            '''
            # create file where removed genes are saved
            _, _, probe_files = next(os.walk(self.dir_probes))
            for probe_file in probe_files:
                gene_id = probe_file[len('probes_'):-len('.txt')]
                if gene_id in self.removed_genes:
                    self.removed_genes.remove(gene_id)
                
            with open(self.file_removed_genes, 'w') as output:
                for gene_id in self.removed_genes:
                    output.write('{}\n'.format(gene_id))

        jobs = []
        for batch_id in range(self.number_batchs):
            proc = multiprocessing.Process(target=_process_blast_results, args=(batch_id, ))
            jobs.append(proc)
            proc.start()

        #print('\n {} \n'.format(jobs))
        
        for job in jobs:
            job.join()
        
        _write_removed_genes()

        # remove intermediate files
        shutil.rmtree(self.dir_blast)
        for file in os.listdir(self.dir_annotations):
            if re.search('probes_*', file):
                os.remove(os.path.join(self.dir_annotations, file))


    

