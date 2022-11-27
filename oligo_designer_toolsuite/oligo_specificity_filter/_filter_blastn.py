import os
import re
from pathlib import Path

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from joblib import Parallel, delayed

from . import SpecificityFilterBase


class Blastn(SpecificityFilterBase):
    """This class filters probes based on the blast alignment tool. All the probes which have a match with a percentage identity higher
    than the one given in input are filtered out.
    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    :param word_size: word size for the blastn seed (exact match to target)
    :type word_size: int
    :param percent_identity: maximum similarity between probes and target sequences, ranging from 0 to 100% (no missmatch)
    :type percent_identity: int
    :param coverage: minimum coverage between probes and target sequence, ranging from 0 to 100% (full coverage)
    :type coverage: int
    """

    def __init__(
        self,
        dir_specificity,
        word_size,
        percent_identity,
        coverage,
    ):
        """Constructor."""
        super().__init__(dir_specificity)

        self.word_size = word_size
        self.percent_identity = percent_identity
        self.coverage = coverage

        self.dir_blast = os.path.join(self.dir_specificity, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, oligo_DB, file_reference_DB, n_jobs):
        """Apply the blastn filter in parallel on the given ``oligo_DB``. Each jobs filters a single gene, and  at the same time are generated at most ``n_job`` jobs.
        The filtered database is returned.

        :param oligo_DB: database containing the probes and their features
        :type oligo_DB: dict
        :param file_reference_DB: path to the file that will be used as reference for the alignement
        :type file_reference_DB: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: probe info of user-specified genes
        :rtype : dict
        """

        # create blast database
        database_exists = False
        database_name = os.path.basename(file_reference_DB)
        # Check if blast database exists
        for file in os.listdir(self.dir_blast):
            if re.search(f"^{database_name}.*", file):
                break

        if not database_exists:
            cmd = NcbimakeblastdbCommandline(
                input_file=file_reference_DB,
                dbtype="nucl",
                out=os.path.join(self.dir_blast, database_name),
            )
            out, err = cmd()

        # run the balst search
        genes = list(oligo_DB.keys())
        filtered_oligo_DBs = Parallel(n_jobs=n_jobs)(
            delayed(self._run_blast)(oligo_DB[gene], gene, database_name)
            for gene in genes
        )

        # reconstruct the oligos db and return it
        for gene, filtered_oligo_DB in zip(genes, filtered_oligo_DBs):
            oligo_DB[gene] = filtered_oligo_DB

        return oligo_DB

    def _run_blast(self, gene_DB, gene, database_name):
        """Run BlastN alignment tool to find regions of local similarity between sequences, where sequences are probes and background sequences (e.g. transcript, genome, etc.).
        BlastN identifies the transcript regions where probes match with a certain coverage and similarity.

        :param gene_DB: database containing the probes form one gene
        :type gene_DB: dict
        :param gene: id of thet gene processed
        :type gene: str
        """
        # run the blast search and write the results
        file_probe_fasta_gene = self._create_fasta_file(gene_DB, self.dir_fasta, gene)
        file_blast_gene = os.path.join(self.dir_blast, f"blast_{gene}.txt")
        cmd = NcbiblastnCommandline(
            query=file_probe_fasta_gene,
            db=os.path.join(self.dir_blast, database_name),
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
        # remove temporary files
        os.remove(os.path.join(self.dir_blast, file_blast_gene))
        os.remove(os.path.join(self.dir_fasta, file_probe_fasta_gene))
        return filtered_gene_DB

    def _read_blast_output(self, file_blast_gene):
        """Load the output of the BlastN alignment search into a DataFrame and process the results."""

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

        :param blast_results: DataFrame with processed blast alignment search results.
        :type blast_results: pandas.DataFrame
        """

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
        """Filer out form the database the sequences with a match.

        :param gene_DB: dictionary with all the probes belonging to the current gene
        :type gene_DB: dict
        :param matching_probes: list of the probes with a match
        :type matching_probes: list
        :return: gene_DB withou the matching probes
        :rtype: dict
        """
        probe_ids = list(gene_DB.keys())
        for probe_id in probe_ids:
            if probe_id in matching_probes:
                del gene_DB[probe_id]
        return gene_DB
