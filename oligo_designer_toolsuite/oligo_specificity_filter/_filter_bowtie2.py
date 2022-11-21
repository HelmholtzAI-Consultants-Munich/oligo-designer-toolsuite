import os
import re
from pathlib import Path
from subprocess import Popen

import pandas as pd
from joblib import Parallel, delayed

from . import SpecificityFilterBase


class Bowtie2(SpecificityFilterBase):
    """This class filters probes based on the Bowtie 2 read alignment tool. It is recommended to use Bowtie 2 instead of Bowtie for reads longer than about 50 bp, as it gives better performance.
    The user can customize the filtering by specifying the min_score. The Bowtie 2 filter gives a alignmnet score to each probe. The higher the score, the more simlar the read sequence is to the reference sequence. The min_score parameter filters out probes with an aligment score greater than min_score

    Use ``conda install -c bioconda bowtie2 to install the Bowtie 2 package``

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    :param min_score: User defined threshhold for alignmnent score. If specified, the Bowtie 2 filter filters out probes with an aligment score greater than min_score. If None, min_score defaults to -0.6 + -0.6 * L, where L is the read length
    :type min_score: float
    """

    def __init__(
        self,
        dir_specificity,
        min_score=None,
    ):
        """Constructor."""
        super().__init__(dir_specificity)

        self.min_score = min_score

        self.dir_bowtie2 = os.path.join(self.dir_specificity, "bowtie2")
        Path(self.dir_bowtie2).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, oligo_DB, file_reference_DB, n_jobs):
        """Apply the bowtie2 filter in parallel on the given ``oligo_DB``. Each jobs filters a single gene, and  at the same time are generated at most ``n_job`` jobs.
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
        # Some bowtie initializations, change the names
        # Check if bowtie2 index exists
        index_exists = False
        index_name = os.path.basename(file_reference_DB)
        for file in os.listdir(self.dir_bowtie2):
            if re.search(f"^{index_name}.*", file):
                index_exists = True
                break

        # Create bowtie2 index if none exists
        if not index_exists:
            command1 = "bowtie2-build -f " + file_reference_DB + " " + index_name
            process = Popen(command1, shell=True, cwd=self.dir_bowtie2).wait()

        genes = list(oligo_DB.keys())
        filtered_oligo_DBs = Parallel(n_jobs=n_jobs)(
            delayed(self._run_bowtie2)(oligo_DB[gene], gene, index_name)
            for gene in genes
        )

        # reconstruct the oligos_DB
        for gene, filtered_oligo_DB in zip(genes, filtered_oligo_DBs):
            oligo_DB[gene] = filtered_oligo_DB
        # remove the files
        """for file in os.listdir(self.dir_bowtie2):
            os.remove(os.path.join(self.dir_bowtie2, file))
        for file in os.listdir(self.dir_fasta):
            os.remove(os.path.join(self.dir_fasta, file))"""
        return oligo_DB

    def _run_bowtie2(self, gene_DB, gene, index_name):
        """Run Bowtie 2 alignment tool to find regions of local similarity between sequences and reference.

        :param gene_DB: database containing the probes form one gene
        :type gene_DB: dict
        :param gene: id of thet gene processed
        :type gene: str
        """

        file_probe_fasta_gene = self._create_fasta_file(gene_DB, self.dir_fasta, gene)
        file_bowtie2_gene = os.path.join(
            self.dir_bowtie2,
            f"bowtie2_{gene}",
        )

        if self.min_score is not None:
            command = (
                "bowtie2 -f -a --score-min "
                + str(self.min_score)
                + " --no-hd --no-unal -x "
                + index_name
                + " -U "
                + file_probe_fasta_gene
                + " -S "
                + file_bowtie2_gene
            )
        else:
            command = (
                "bowtie2 -f -a "
                + " --no-hd --no-unal -x "
                + index_name
                + " -U "
                + file_probe_fasta_gene
                + " -S "
                + file_bowtie2_gene
            )
        process = Popen(command, shell=True, cwd=self.dir_bowtie2).wait()

        # read the results of the bowtie search
        bowtie2_results = self._read_bowtie2_output(file_bowtie2_gene)
        # filter the DB based on the bowtie results
        matching_probes = self._find_matching_probes(bowtie2_results)
        filtered_gene_DB = self._filter_matching_probes(gene_DB, matching_probes)
        os.remove(os.path.join(self.dir_bowtie2, file_bowtie2_gene))
        os.remove(os.path.join(self.dir_fasta, file_probe_fasta_gene))
        return filtered_gene_DB

    def _read_bowtie2_output(self, file_bowtie2_gene):
        """Load the output of the bowtie 2 alignment search into a DataFrame and process the results."""
        bowtie2_results = pd.read_csv(
            file_bowtie2_gene,
            header=None,
            sep="\t",
            low_memory=False,
            names=[
                "query",
                "flags",
                "reference",
                "ref_start",
                "mapping_quality",
                "CIGAR_alignment",
                "ref_seq_mate",
                "offset",
                "fragment_length",
                "sequence",
                "read_qualities",
                "alignment_score",
            ],
            engine="c",
            usecols=list(range(12)),
            dtype={
                "query": str,
                "flags": str,
                "reference": str,
                "ref_start": str,
                "mapping_quality": str,
                "CIGAR_alignment": str,
                "ref_seq_mate": str,
                "offset": str,
                "fragment_length": str,
                "sequence": str,
                "read_qualities": str,
                "alignment_score": str,
            },
        )

        bowtie2_results["query_gene_id"] = (
            bowtie2_results["query"].str.split("_").str[0]
        )
        bowtie2_results["reference_gene_id"] = (
            bowtie2_results["reference"].str.split("::").str[0]
        )
        return bowtie2_results

    def _find_matching_probes(self, bowtie2_results):
        """Use the results of the Bowtie 2 alignement search to remove probes with high similarity (i.e. low number of mismatches represented by a high alignment score) based on user-defined min_score.

        :param bowtie2_results: DataFrame with processed bowtie alignment search results.
        :type bowtie2_results: pandas.DataFrame
        """

        bowtie2_matches = bowtie2_results[
            bowtie2_results["query_gene_id"] != bowtie2_results["reference_gene_id"]
        ]
        probes_with_match = bowtie2_matches["query"].unique()
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
