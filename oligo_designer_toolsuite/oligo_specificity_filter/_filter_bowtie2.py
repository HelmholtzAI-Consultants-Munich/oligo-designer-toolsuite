############################################
# imports
############################################

import os
import re
from pathlib import Path
from subprocess import Popen

import pandas as pd
from joblib import Parallel, delayed

from . import Bowtie

############################################
# Oligo Bowtie2 Filter Classes
############################################


class Bowtie2(Bowtie):
    """This class filters oligos based on the Bowtie 2 read alignment tool. It is recommended to use Bowtie 2 instead of Bowtie for reads longer than about 50 bp, as it gives better performance.
    The user can customize the filtering by specifying the min_score. The Bowtie 2 filter gives an alignment score to each oligo. The higher the score, the more similar the read sequence is to the reference sequence.
    The min_score parameter filters out oligos with an alignment score greater than min_score

    Use ``conda install -c bioconda bowtie2`` to install the Bowtie 2 package

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    :param min_score: User defined threshold for alignment score. If specified, the Bowtie 2 filter filters out oligos with an alignment score greater than min_score. If None, min_score defaults to -0.6 + -0.6 * L, where L is the read length
    :type min_score: float
    """

    def __init__(
        self,
        dir_specificity: str,
        min_score: float = None,
    ):
        """Constructor."""
        super().__init__(dir_specificity)

        self.min_score = min_score

        self.dir_bowtie2 = os.path.join(self.dir_specificity, "bowtie2")
        Path(self.dir_bowtie2).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """Apply the bowtie2 filter in parallel on the given ``database``. Each jobs filters a single region, and  at the same time are generated at most ``n_job`` jobs.
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
        # Some bowtie initializations, change the names
        # Check if bowtie2 index exists
        index_exists = False
        index_name = os.path.basename(file_reference)
        for file in os.listdir(self.dir_bowtie2):
            if re.search(f"^{index_name}.*", file):
                index_exists = True
                break

        # Create bowtie2 index if none exists
        if not index_exists:
            command1 = (
                "bowtie2-build --quiet --threads "
                + str(n_jobs)
                + " -f "
                + file_reference
                + " "
                + index_name
            )
            process = Popen(command1, shell=True, cwd=self.dir_bowtie2).wait()

        regions = list(database.keys())
        filtered_database_regions = Parallel(n_jobs=n_jobs)(
            delayed(self._run_bowtie2)(database[region], region, index_name)
            for region in regions
        )

        # reconstruct the oligos_DB
        for region, filtered_database_region in zip(regions, filtered_database_regions):
            database[region] = filtered_database_region

        return database

    def _run_bowtie2(self, database_region, region, index_name):
        """Run Bowtie 2 alignment tool to find regions of local similarity between sequences and reference.

        :param database_region: database containing the oligos form one region
        :type database_region: dict
        :param region: id of the region processed
        :type region: str
        """

        file_oligo_fasta_gene = self._create_fasta_file(
            database_region, self.dir_fasta, region
        )
        file_bowtie2_gene = os.path.join(
            self.dir_bowtie2,
            f"bowtie2_{region}",
        )

        if self.min_score is not None:
            command = (
                "bowtie2 --quiet -f -a --score-min "
                + str(self.min_score)
                + " --no-hd --no-unal -x "
                + index_name
                + " -U "
                + file_oligo_fasta_gene
                + " -S "
                + file_bowtie2_gene
            )
        else:
            command = (
                "bowtie2 --quiet -f -a "
                + " --no-hd --no-unal -x "
                + index_name
                + " -U "
                + file_oligo_fasta_gene
                + " -S "
                + file_bowtie2_gene
            )
        process = Popen(command, shell=True, cwd=self.dir_bowtie2).wait()

        # read the results of the bowtie search
        bowtie2_results = self._read_bowtie2_output(file_bowtie2_gene)
        # filter the DB based on the bowtie results
        matching_oligos = self._find_matching_oligos(bowtie2_results)
        filtered_database_region = self._filter_matching_oligos(
            database_region, matching_oligos
        )
        # remove the temporary files
        os.remove(os.path.join(self.dir_bowtie2, file_bowtie2_gene))
        os.remove(os.path.join(self.dir_fasta, file_oligo_fasta_gene))
        return filtered_database_region

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
