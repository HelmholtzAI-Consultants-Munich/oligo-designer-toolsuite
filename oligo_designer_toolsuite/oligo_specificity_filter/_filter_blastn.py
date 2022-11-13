import os
from pathlib import Path

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from joblib import Parallel, delayed

from . import SpecificityFilterBase


class Blastn(SpecificityFilterBase):
    def __init__(
        self,
        dir_specificity,
        word_size,
        percent_identity,
        coverage,
        probe_length_min,
        probe_length_max,
    ):
        """This class filters probes based on the blast alignment tool

        :param file_reference_DB: path to fasta file containing all probes
        :type file_reference_DB: str
        :param word_size: word size for the blastn seed (exact match to target)
        :type word_size: int
        :param percent_identity: maximum similarity between probes and target sequences, ranging from 0 to 100% (no missmatch)
        :type percent_identity: int
        :param coverage: minimum coverage between probes and target sequence, ranging from 0 to 100% (full coverage)
        :type coverage: int
        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int
        """
        super().__init__(dir_specificity)

        self.word_size = word_size
        self.percent_identity = percent_identity
        self.coverage = coverage
        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max

        self.dir_blast = os.path.join(self.dir_specificity, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, oligo_DB, file_reference_DB, n_jobs):
        """Apply blastN filter to all batches in parallel"""

        # create blast database
        cmd = NcbimakeblastdbCommandline(
            input_file=file_reference_DB,
            dbtype="nucl",
            out=os.path.join(self.dir_blast, "Blast_DB"),
        )
        out, err = cmd()

        # run the balst search
        genes = list(oligo_DB.keys())
        filtered_oligo_DBs = Parallel(n_jobs=n_jobs)(
            delayed(self._run_blast)(oligo_DB[gene], gene) for gene in genes
        )

        # reconstruct the oligos db and return it
        for gene, filtered_oligo_DB in zip(genes, filtered_oligo_DBs):
            oligo_DB[gene] = filtered_oligo_DB
        # remove the files
        for file in os.listdir(self.dir_blast):
            os.remove(os.path.join(self.dir_blast, file))
        for file in os.listdir(self.dir_fasta):
            os.remove(os.path.join(self.dir_fasta, file))
        return oligo_DB

    def _run_blast(self, gene_DB, gene):
        """Run BlastN alignment tool to find regions of local similarity between sequences, where sequences are probes and transcripts.
        BlastN identifies the transcript regions where probes match with a certain coverage and similarity.
        """
        # run the blast search and write the results
        file_probe_fasta_gene = self._create_fasta_file(gene_DB, self.dir_fasta, gene)
        file_blast_gene = os.path.join(self.dir_blast, f"blast_{gene}.txt")
        cmd = NcbiblastnCommandline(
            query=file_probe_fasta_gene,
            db=os.path.join(self.dir_blast, "Blast_DB"),
            outfmt="10 qseqid sseqid length qstart qend qlen",
            out=file_blast_gene,
            strand="plus",
            word_size=self.word_size,
            perc_identity=self.percent_identity,
            num_threads=1,  # ????
        )
        out, err = cmd()

        # read the reuslts of the blast seatch
        blast_results = self._read_blast_output(file_blast_gene)
        # filter the DB based on the blast results
        matching_probes = self._find_matching_probes(blast_results)
        filtered_gene_DB = self._filter_matching_probes(gene_DB, matching_probes)
        return filtered_gene_DB

    def _read_blast_output(self, file_blast_gene):
        """Load the output of the BlastN alignment search into a DataFrame and process the results.
        :return: DataFrame with processed blast alignment search results.
        :rtype: pandas.DataFrame
        """

        blast_results = pd.read_csv(
            file_blast_gene,
            header=None,
            sep=",",
            low_memory=False,
            names=[
                "query",
                "target",
                "alignment_length",
                "query_start",
                "query_end",
                "query_length",
            ],
            engine="c",
            dtype={
                "query": str,
                "target": str,
                "alignment_length": int,
                "query_start": int,
                "query_end": int,
                "query_length": int,
            },
        )
        # return the real matches, that is the ones not belonging to the same gene of the query probe
        blast_results["query_gene_id"] = blast_results["query"].str.split("_").str[0]
        blast_results["target_gene_id"] = blast_results["target"].str.split("::").str[0]
        return blast_results

    def _find_matching_probes(self, blast_results):
        """Use the results of the BlastN alignement search to remove probes with high similarity,
        probe coverage and ligation site coverage based on user defined thresholds.
        :param probes_info: Dataframe with probe information, filtered based on sequence properties.
        :type probes_info: pandas.DataFrame
        :param blast_results: DataFrame with processed blast alignment search results.
        :type blast_results: pandas.DataFrame
        """

        # TODO: implement another option since we don't necessary have te lagation site
        blast_matches = blast_results[
            blast_results["query_gene_id"] != blast_results["target_gene_id"]
        ]
        values = blast_matches["query_length"] * self.coverage / 100
        blast_matches.insert(len(blast_matches.columns), "min_alignment_length", values)

        blast_matches_filtered = blast_matches.loc[
            blast_matches.alignment_length > blast_matches.min_alignment_length
        ]
        probes_with_match = blast_matches_filtered["query"].unique()

        return probes_with_match

    def _filter_matching_probes(self, gene_DB, matching_probes):
        probe_ids = list(gene_DB.keys())
        for probe_id in probe_ids:
            if probe_id in matching_probes:
                del gene_DB[probe_id]
        return gene_DB
