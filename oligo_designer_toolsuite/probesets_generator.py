############################################
# imports
############################################

import os
import logging
import shutil
import fnmatch
import multiprocessing

import numpy as np
import pandas as pd
import networkx as nx

from pathlib import Path
from functools import reduce

############################################
# probe set generator class
############################################

class ProbesetsGenerator:
    '''This class is used to generate ranked, non-overlapping probe sets.
    '''

    def __init__(self, config, dir_output, dir_probes = None):
        """Constructor method
        """
         # set logger
        self.logging = logging.getLogger('padlock_probe_designer')

        # set directory
        if dir_probes == None:
            self.dir_probes = os.path.join(dir_output, 'probes')
            Path(self.dir_probes).mkdir(parents=True, exist_ok=True)
        else:
            self.dir_probes = dir_probes

        self.dir_overlapmatrix = os.path.join(dir_output, 'overlapmatrix')
        Path(self.dir_overlapmatrix).mkdir(parents=True, exist_ok=True)

        self.dir_probesets = os.path.join(dir_output, 'probesets')
        Path(self.dir_probesets).mkdir(parents=True, exist_ok=True)

        self.file_removed_genes = os.path.join(dir_output, 'genes_with_insufficient_probes.txt')

        # set parameters
        self.n_probes_per_gene = config['n_probes_per_gene']
        self.Tm_params = [config['Tm_opt'], config['Tm_min'], config['Tm_max']]
        self.GC_params = [config['GC_content_opt'], config['GC_content_min'], config['GC_content_max']]
        self.min_probes_per_gene = config['min_probes_per_gene']


    def get_overlap_matrix(self):
        """Generate overlap matrix for probes of each gene in directory dir_in
        """  
        def _list_files_in_dir(dir, pattern):
            for root, dirs, files, rootfd in os.fwalk(dir, follow_symlinks=True):
                for basename in files:
                    if fnmatch.fnmatch(basename, pattern):
                        filename = os.path.join(root, basename)
                        yield filename


        def _get_files():
            files = []
            number_dirs = 0
            #Note: the other functions here work with multiple `dirs_in` which is not needed when we have all probes files at one
            #place. Could be changed for the other functions in the future.
            for dir in [self.dir_probes]:
                files.append(pd.DataFrame([[file.split('probes_')[1].split('.')[0], file] for file in _list_files_in_dir(dir, r'probes_*')], columns=['gene', 'file_dir{}'.format(number_dirs)]))
                number_dirs += 1
            files = reduce(lambda  left,right: pd.merge(left, right, on=['gene'], how='outer'), files)
            files = files.loc[:,~files.columns.duplicated()]

            return files


        def _get_overlap(seq1_intervals, seq2_intervals):
            seqs_overlap = False
            for a in seq1_intervals:
                for b in seq2_intervals:
                    overlap = min(a[1], b[1]) - max(a[0], b[0])
                    seqs_overlap |= overlap > -1
            return seqs_overlap


        def _compute_overlap_matrix(gene, probes):
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
            matrix.to_csv(os.path.join(self.dir_overlapmatrix, 'overlap_matrix_{}.txt'.format(gene)), sep='\t')


        def _get_overlap_matrix(files):

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
                #probes.to_csv(os.path.join(self.dir_output_overlapmatrix, 'probes_{}.txt'.format(gene)), sep='\t', index=False)
            
                proc = multiprocessing.Process(target=_compute_overlap_matrix, args=(gene, probes, ))
                jobs.append(proc)
                proc.start()

            for job in jobs:
                job.join()

        files = _get_files()
        _get_overlap_matrix(files)


    def get_probe_sets(self, n_sets=5000):
        """Get ranked list of non overlapping probe sets for each gene in directory `dir_in_probes`.
        
        n_sets: int
            Maximum number of ranked returned sets.
            
        Saves table "ranked_probesets_<gene>.txt" with ranked sets of probe ids
        """
        # TODO: Maybe we want to add the weighting between GC content and Tm in the scoring function as parameters in the config
        #       and additionally we might want to optimise also based on the blast result parameters (max overlap/coverage) and
        #       accordingl add another weighting factor
        def _get_probe_score_function(TM_w=1,GC_w=1):
            """Get function for scoring individual probes
            
            Note: initially thought about scoring sets and not single probes, if there is a reason for scoring probe sets (e.g. 
            if you're searching for probes with most similar melting temperature etc.) then we might want to adapt this, but for
            now we're scoring probes and choose the set with the lowest max score.
            
            Tm_params: list [Tm_opt,Tm_min,Tm_max] (melting temperature optimum and interval)
            GC_params: list [GC_opt,GC_min,GC_max] (GC contant optimum and interval)
            Tm_w: weight of Tm
            GC_w: wight of GC
            
            """
            Tm_dif_max = self.Tm_params[2] - self.Tm_params[0]
            Tm_dif_min = self.Tm_params[0] - self.Tm_params[1]
            if Tm_dif_max == Tm_dif_min:
                Tm_error = lambda T: abs(T - self.Tm_params[0])/Tm_dif_max
            else:
                def Tm_error(T):
                    T_dif = T - self.Tm_params[0]
                    return abs(T_dif)/Tm_dif_max * (T_dif>0) + abs(T_dif)/Tm_dif_min * (T_dif<0)
                
            GC_dif_max = self.GC_params[2] - self.GC_params[0]
            GC_dif_min = self.GC_params[0] - self.GC_params[1]
            if GC_dif_max == GC_dif_min:
                GC_error = lambda gc: abs(gc - self.GC_params[0])/GC_dif_max
            else:
                def GC_error(gc):
                    gc_dif = gc - self.GC_params[0]
                    return abs(gc_dif)/GC_dif_max * (gc_dif>0) + abs(gc_dif)/GC_dif_min * (gc_dif<0)
            
            def get_probe_score(Tm, GC):
                """Score probe based on melting temperature and GC content.
                Only positive scores possible. Score 0 is the best score.
                """
                return TM_w * Tm_error(Tm) + GC_w * GC_error(GC)
            
            return get_probe_score            


        def _select_n_probes_by_heuristic(probes, adj_mat, n=5, n_trials=10000):
            """Get a good probe set via a heuristic search
            
            This heurisitc only works if we assign single scores to each gene (which we generally assume for now).
            And it's especially simple by optimising for minimal max score of a gene set (instead of e.g. the sum).
            Last sentence refers to the filtering of probes in `get_nonoverlapping_sets`: If we optimize for the minimal max 
            score we can just throw out those probes that have a higher score than the worst probe of the heuristic set.
            
            Heuristic:
            Choose the gene with the lowest score, choose the next gene not overlapping with the first and the remaining
            best score, ... go on till a set of n genes is selected. To increase the probably of finding a good set the 
            approach is repeated, starting with the 2nd, the 3rd, ... the n_trials'th best gene. 
            
            probes: pd.DataFrame
                index contains probe ids, columns "score" contains scores of probes
            adj_mat: pd.DataFrame
                Overlap matrix of probes
            n: int
                Number of probes to select
            n_trials: int
                Number of trails to find a good set via the heuristic approach
            
            Returns
            -------
            list
                Probe ids
            float
                Highest (i.e. worst) score of selected probes
            
            """
            max_score = probes["score"].max() * 1.1
            probes_sorted = probes.sort_values("score")
            probe_ids_sorted = probes_sorted.index.tolist()
            
            inv_mat_sorted = adj_mat.loc[probe_ids_sorted,probe_ids_sorted].values
            inv_mat_sorted = np.abs(inv_mat_sorted - 1)
            
            for first_idx, probe_id in enumerate(probe_ids_sorted[:min(len(probe_ids_sorted),n_trials)]):
                set_idxs = np.array([first_idx])
                for _ in range(n-1):
                    # find first probe in sorted array that is not overlapping with any selected probe
                    no_overlap = np.all(inv_mat_sorted[set_idxs],axis=0)
                    if np.any(no_overlap):
                        set_idxs = np.append(set_idxs,np.where(no_overlap)[0][0])
                    else:
                        break
                if len(set_idxs) == n:
                    score = np.max(probes_sorted["score"].values[set_idxs])
                    if score < max_score:
                        max_score = score
                        best_set = set_idxs
                        
            return np.array(probe_ids_sorted)[best_set].tolist(), max_score


        def _get_nonoverlapping_sets(probes, adj_mat, n_sets=100):
            """Get ranked list of non overlapping probe sets for one gene
            
            config: config file which contains parameters for best set criteria (Tm and GC)
            adj_mat: pd.DataFrame overlap matrix of probes with probe_ids in index and columns.
            probes: pd.DataFrame with melting temperature of probes
            n_sets: maximal number of sets the will be returned.
            
            Returns
            -------
            pd.DataFrame 
                It contains up to `n_sets` probesets as rows. The number of probes per probeset equals
                config["n_probes_per_gene"] or the next highest number that could be found. The probesets are
                ranked as follows: 
                Each probe has a score, the maximum score (i.e. the worst) defines the overall score of the
                probeset. Additionally sets with the same score are further ranked by the sum over scores of 
                probes. --> The output contains the according columns ["score","score_sum"].
            """
            n = self.n_probes_per_gene
            # Get scoring function for probes
            get_probe_score = _get_probe_score_function(TM_w=1, GC_w=1)
            
            # Score probes
            probes["score"] = probes[["melting_temperature","GC_content"]].apply(
                lambda TGC: get_probe_score(TGC[0], TGC[1]), axis=1
            )
            
            # Represent overlap matrix as graph
            G = nx.convert_matrix.from_numpy_matrix(adj_mat.values)
            GC = nx.algorithms.operators.unary.complement(G)
            
            # Initialize variable
            heuristic_set = None
            
            # First check if there are no cliques with n probes
            n_ = 0
            for clique in nx.algorithms.clique.find_cliques(GC):
                n_ = max(len(clique),n_)
                if n_ >= n:
                    break
            if n_ < n:
                n = n_
            else:
                # Run a heuristic search first to reduce the number of probes
                heuristic_set, max_error = _select_n_probes_by_heuristic(probes, adj_mat, n=n, n_trials=10000)
                probes_below_err = (probes["score"] <= max_error)
                probes = probes.loc[probes_below_err]
                G = nx.convert_matrix.from_numpy_matrix(adj_mat.loc[probes_below_err,probes_below_err].values)
                GC = nx.algorithms.operators.unary.complement(G)
            
            # Search for best set
            probesets = []
            count = 0
            # Note: Search could be further optimised by iteratively throwing out probes with worse scores then current best set
            for clique in nx.algorithms.clique.find_cliques(GC):
                count += 1
                # Limit the number of combinations we iterate through
                if count > 100000:
                    break
                if len(clique) >= n:
                    # Get probe_ids of clique
                    probe_ids = probes.iloc[clique].index.tolist()
                    best_n_probes = probes.loc[probe_ids,"score"].sort_values(ascending = True).head(n)#.index.tolist()
                    # Calculate performance of probeset
                    probeset_error = best_n_probes.max()
                    probesets.append(best_n_probes.index.tolist() + [probeset_error,best_n_probes.sum()])
            probesets = pd.DataFrame(columns=[f"probe_{i}" for i in range(n)] + ["score","score_sum"], data=probesets)
            probesets = probesets.sort_values(['score','score_sum'], ascending = [True,True])
            # Eventually add heuristic set if it's better than the best set found (can only happen for very high nrs of probes)
            if heuristic_set:
                error_sum = probes.loc[heuristic_set,"score"].sum()
                heuristic_is_better = (max_error < probesets["score"].min())
                h_is_better_wrt_score2 = (max_error == probesets["score"].min()) and (error_sum < probesets["score_sum"].min())
                if heuristic_is_better or h_is_better_wrt_score2:
                    probesets.loc[len(probesets)] = heuristic_set + [max_error,error_sum]
            if n_sets:
                probesets = probesets.head(n_sets)
            probesets = probesets.reset_index(drop=True)
            
            return probesets


        probes_files = [f for f in os.listdir(self.dir_probes) if f.startswith('probes_')]
        genes = [f.split('_')[1].split(".")[0] for f in probes_files]
        overlap_files = [f'overlap_matrix_{gene_id}.txt' for gene_id in genes]
        
        for gene_id, probe_f, overlap_f in zip(genes,probes_files, overlap_files):
            probes = pd.read_csv(os.path.join(self.dir_probes, probe_f),index_col=0,sep='\t')
            adj_mat = pd.read_csv(os.path.join(self.dir_overlapmatrix, overlap_f),index_col=0,sep='\t')
            table_probesets = _get_nonoverlapping_sets(probes, adj_mat, n_sets=n_sets)
            probes_per_gene = len([col for col in table_probesets.columns if col.startswith('probe')])
            if probes_per_gene < self.min_probes_per_gene:
                with open(self.file_removed_genes, 'a') as output:
                    output.write(f'{gene_id}\t{probes_per_gene}\n')
            else:
                table_probesets.to_csv(os.path.join(self.dir_probesets, f'ranked_probesets_{gene_id}.txt'), sep='\t')

        shutil.rmtree(self.dir_overlapmatrix)
