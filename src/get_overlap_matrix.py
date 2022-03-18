############################################
# imports
############################################

import os
import argparse
import fnmatch
import pandas as pd
from functools import reduce
from pathlib import Path

import multiprocessing

############################################
# functions
############################################

def _get_files(dirs_in):
    
    files = []
    number_dirs = 0
    for dir in dirs_in:
        files.append(pd.DataFrame([[file.split('probes_')[1].split('.')[0], file] for file in _list_files_in_dir(dir, r'probes_*')], columns=['gene', 'file_dir{}'.format(number_dirs)]))
        number_dirs += 1
    files = reduce(lambda  left,right: pd.merge(left, right, on=['gene'], how='outer'), files)
    files = files.loc[:,~files.columns.duplicated()]

    return files


############################################

def _list_files_in_dir(dir, pattern):
    for root, dirs, files, rootfd in os.fwalk(dir, follow_symlinks=True):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


############################################

def _get_overlap_matrix(files, dir_out):

    jobs = []
    for idx in files.index:
        files_gene = files.iloc[idx]
        for index, value in files_gene.items():
            if index == 'gene':
                gene = value
                probes = []
            elif not pd.isna(value):
                probes.append(pd.read_csv(value, sep='\t'))
        probes = pd.concat(probes, axis=0, ignore_index=True)
        probes.to_csv(os.path.join(dir_out, 'probes_{}.txt'.format(gene)), sep='\t', index=False)
    
        proc = multiprocessing.Process(target=_compute_overlap_matrix, args=(gene, probes, dir_out, ))
        jobs.append(proc)
        proc.start()

    for job in jobs:
        job.join()
        

############################################

def _compute_overlap_matrix(gene, probes, dir_out):
    matrix = pd.DataFrame(0, columns=probes.probe_id, index=probes.probe_id)
    for i in probes.index:
        probe1_starts =  [int(s) for s in str(probes.loc[i,'start']).split(";")]
        probe1_ends =  [int(s) for s in str(probes.loc[i,'end']).split(";")]
        probe1_intervals = [[start,end] for start,end in zip(probe1_starts, probe1_ends)]
        pid1 = probes.loc[i,'probe_id']
        for j in probes.index:
            probe2_starts =  [int(s) for s in str(probes.loc[j,'start']).split(";")]
            probe2_ends =  [int(s) for s in str(probes.loc[j,'end']).split(";")]
            probe2_intervals = [[start,end] for start,end in zip(probe2_starts, probe2_ends)]
            pid2 = probes.loc[j,'probe_id']
            if _get_overlap(probe1_intervals, probe2_intervals):
                matrix.loc[pid1,pid2] = 1
                matrix.loc[pid2,pid1] = 1
            else:
                matrix.loc[pid1,pid2] = 0
                matrix.loc[pid2,pid1] = 0
            if j > i:
                break
    matrix.to_csv(os.path.join(dir_out, 'overlap_matrix_{}.txt'.format(gene)), sep='\t')


############################################

def _get_overlap(seq1_intervals, seq2_intervals):
    seqs_overlap = False
    for a in seq1_intervals:
        for b in seq2_intervals:
            overlap = min(a[1], b[1]) - max(a[0], b[0])
            seqs_overlap |= overlap > -1
    return seqs_overlap


############################################

def get_overlap_matrix(dir_in,dir_out):
    """Generate overlap matrix for probes of each gene in directory dir_in
    
    Note: the other functions here work with multiple `dirs_in` which is not needed when we have all probes files at one
    place. Could be changed for the other functions in the future.
    """

    Path(dir_out).mkdir(parents=True, exist_ok=True)    
    
    files = _get_files([dir_in])
    _get_overlap_matrix(files, dir_out)

