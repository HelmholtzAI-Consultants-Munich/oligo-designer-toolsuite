import os
import re
from pathlib import Path

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from joblib import Parallel, delayed

from . import SpecificityFilterBase


class Blastn(SpecificityFilterBase):
    """This class filters oligos based on the blast alignment tool. All the oligos which have a match with a percentage identity higher
    than the one given in input are filtered out.

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    :param word_size: word size for the blastn seed (exact match to target)
    :type word_size: int
    :param percent_identity: maximum similarity between oligos and target sequences, ranging from 0 to 100% (no missmatch)
    :type percent_identity: int
    :param coverage: minimum coverage between oligos and target sequence, ranging from 0 to 100% (full coverage)
    :type coverage: int
    :param strand: strand of the query sequence to search
    :type strand: str
    """

    def __init__(
        self,
        dir_specificity: str,
        word_size: int,
        percent_identity: float,
        coverage: float,
        strand: str,
    ):
        """Constructor."""
        super().__init__(dir_specificity)

        self.word_size = word_size
        self.percent_identity = percent_identity
        self.coverage = coverage
        self.strand = strand

        self.dir_blast = os.path.join(self.dir_specificity, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """Apply the blastn filter in parallel on the given ``database``. Each jobs filters a single region, and  at the same time are generated at most ``n_job`` jobs.
        The filtered database is returned.

        :param database: database containing the oligos and their features
        :type database: dict
        :param file_reference: path to the file that will be used as reference for the alignement
        :type file_reference: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: oligo info of user-specified regions
        :rtype: dict
        """

        # create blast database
        database_exists = False
        database_name = os.path.basename(file_reference)
        # Check if blast database exists
        for file in os.listdir(self.dir_blast):
            if re.search(f"^{database_name}.*", file):
                database_exists = True
                break

        if not database_exists:
            cmd = NcbimakeblastdbCommandline(
                input_file=file_reference,
                dbtype="nucl",
                out=os.path.join(self.dir_blast, database_name),
            )
            out, err = cmd()

        # run the balst search
        regions = list(database.keys())
        filtered_database_regions = Parallel(n_jobs=n_jobs)(
            delayed(self._run_blast)(database[region], region, database_name)
            for region in regions
        )

        # reconstruct the oligos db and return it
        for region, filtered_database_region in zip(regions, filtered_database_regions):
            database[region] = filtered_database_region

        return database

    def _run_blast(self, database_region, region, database_name):
        """Run BlastN alignment tool to find regions of local similarity between sequences, where sequences are oligos and background sequences (e.g. transcript, genome, etc.).
        BlastN identifies the transcript regions where oligos match with a certain coverage and similarity.

        :param database_region: database containing the oligos form one region
        :type database_region: dict
        :param region: id of the region processed
        :type region: str
        :param database_name: path to the blatn database
        :type database_name: str
        """
        # run the blast search and write the results
        file_oligo_fasta_gene = self._create_fasta_file(database_region, self.dir_fasta, region)
        file_blast_gene = os.path.join(self.dir_blast, f"blast_{region}.txt")
        cmd = NcbiblastnCommandline(
            query=file_oligo_fasta_gene,
            db=os.path.join(self.dir_blast, database_name),
            outfmt="6 qseqid sseqid length qstart qend qlen",
            out=file_blast_gene,
            strand=self.strand,
            word_size=self.word_size,
            perc_identity=self.percent_identity,
            num_threads=1,  # ????
        )
        out, err = cmd()

        # read the reuslts of the blast seatch
        blast_results = self._read_blast_output(file_blast_gene)
        # filter the DB based on the blast results
        matching_oligos = self._find_matching_oligos(blast_results)
        filtered_database_region = self._filter_matching_oligos(database_region, matching_oligos)
        # remove temporary files
        #os.remove(os.path.join(self.dir_blast, file_blast_gene))
        os.remove(os.path.join(file_blast_gene))

        # os.remove(os.path.join(self.dir_fasta, file_oligo_fasta_gene))
        os.remove(os.path.join(file_oligo_fasta_gene))
        return filtered_database_region

    def _read_blast_output(self, file_blast_gene):
        """Load the output of the BlastN alignment search into a DataFrame and process the results."""

        blast_results = pd.read_csv(
            file_blast_gene,
            header=None,
            sep="\t",
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
        # return the real matches, that is the ones not belonging to the same region of the query oligo
        blast_results["query_gene_id"] = blast_results["query"].str.split("_").str[0]
        blast_results["target_gene_id"] = blast_results["target"].str.split("::").str[0]
        return blast_results

    def _find_matching_oligos(self, blast_results):
        """Use the results of the BlastN alignement search to remove oligos with high similarity,
        oligo coverage and ligation site coverage based on user defined thresholds.

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

        oligos_with_match = blast_matches_filtered["query"].unique()

        return oligos_with_match
