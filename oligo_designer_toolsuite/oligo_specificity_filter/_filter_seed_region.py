############################################
# imports
############################################

import os
import re
import subprocess

from pathlib import Path
from joblib import Parallel, delayed

from . import Bowtie

############################################
# Oligo Seed Region Filter Classes
############################################

class BowtieSeedRegion(Bowtie):
    """This class filters oligos based on the Bowtie short read alignment tool on a specific sub-region of the oligo. The region taken in consideration is created according to the ``seed_region_creation`` class.
    The user can customize the filtering by specifying the num_mismatches, and all oligos with number mismatches lower or equal to num_mismatches inside the seed region are filtered out.

    Use conda install -c bioconda bowtie to install Bowtie package

    :param file_transcriptome_fasta: path to fasta file containing all oligos
    :type file_transcriptome_fasta: str
    :param seed_region_creation: Class to generate the region of the oligo where the mismatches are considered. Oligos that have less than min_mismatches in the seed region are filtered out
    :type seed_region_creation: SeedRegionCreationBase class
    :param num_mismatches: Threshhold value on the number of mismatches required for each oligo. Oligos where the number of mismatches are greater than or equal to this threshhold are considered valid. Possible values range from 0 to 3.
    :type num_mismatches: int
    """

    def __init__(
        self,
        dir_specificity: str,
        seed_region_creation,
        num_mismatches: int = 0,
        strand: str = None,
    ):
        """Constructor."""
        super().__init__(dir_specificity)

        if num_mismatches > 3:
            raise ValueError(
                "Choice of num_mismatches out of range for bowtie allignment tool. Please choose a value no greater than 3"
            )
        else:
            self.num_mismatches = num_mismatches

        self.seed_region_creation = seed_region_creation
        self.strand = strand

        self.dir_seed_region = os.path.join(
            self.dir_specificity, "bowtie"
        )  # ecplot some possible already computed indices
        Path(self.dir_seed_region).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, database: dict, file_reference_DB: str, n_jobs: int):
        """Apply the bowtie filter in parallel on the seed region of oligos in the given ``database``. Each jobs filters a single region, and  at the same time are generated at most ``n_job`` jobs.
        The filtered database is returned.

        :param database: database containing the oligos and their features
        :type database: dict
        :param file_reference_DB: path to the file that will be used as reference for the alignement
        :type file_reference_DB: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: oligo info of user-specified regions
        :rtype: dict"""
        # generater the seed region coordinates
        database = self.seed_region_creation.apply(database)

        # Some bowtie initializations, change the names
        index_exists = False
        index_name = os.path.basename(file_reference_DB)
        # Check if bowtie index exists
        for file in os.listdir(self.dir_seed_region):
            if re.search(f"^{index_name}.*", file):
                index_exists = True
                break
        # Create bowtie index if none exists
        if not index_exists:
            command1 = (
                "bowtie-build --quiet --threads "
                + str(n_jobs)
                + " -f "
                + file_reference_DB
                + " "
                + index_name
            )
            process = subprocess.Popen(
                command1, shell=True, cwd=self.dir_seed_region
            ).wait()

        oligo_DB_seed = self._extract_seed_regions(database)
        regions = list(oligo_DB_seed.keys())
        filtered_oligo_DBs = Parallel(n_jobs=n_jobs)(
            delayed(self._run_bowtie_seed_region)(
                oligo_DB_seed[region], database[region], region, index_name
            )
            for region in regions
        )

        # reconstruct the oligos_DB
        for region, filtered_oligo_DB in zip(regions, filtered_oligo_DBs):
            database[region] = filtered_oligo_DB

        return database

    def _run_bowtie_seed_region(
        self, database_region_seed, gene_DB, region, index_name
    ):
        """Run Bowtie alignment tool to find regions of local similarity between sequences, where sequences are oligos and transcripts.
        Bowtie identifies all allignments between the oligos and transcripts and returns the number of mismatches and mismatch position for each alignment.

        :return: DataFrame with processed bowtie alignment search results.
        :rtype: pandas.DataFrame
        """

        file_oligo_fasta_gene = self._create_fasta_file(
            database_region_seed, self.dir_fasta, region
        )
        file_bowtie_gene = os.path.join(
            self.dir_seed_region,
            f"bowtie_{region}.txt",
        )
        command = (
            "bowtie --quiet -x "
            + index_name
            + " -f -a -v "
            + str(self.num_mismatches)
            + " "
            + file_oligo_fasta_gene
            + " "
            + file_bowtie_gene
        )

        process = subprocess.Popen(
            command,
            shell=True,
            cwd=self.dir_seed_region,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).wait()

        # read the results of the bowtie search
        bowtie_results = self._read_bowtie_output(file_bowtie_gene)
        # filter the DB based on the bowtie results
        matching_oligos = self._find_matching_oligos(bowtie_results)
        filtered_database_region = self._filter_matching_oligos(
            gene_DB, matching_oligos
        )
        # remove the temporary files
        os.remove(os.path.join(self.dir_seed_region, file_bowtie_gene))
        os.remove(os.path.join(self.dir_fasta, file_oligo_fasta_gene))
        return filtered_database_region

    def _extract_seed_regions(self, database):
        """geneate a new oligos DB containing only the seed regions of the oligos."""
        oligo_DB_seed = {}
        for region in database.keys():
            oligo_DB_seed[region] = {}
            for oligo_id in database[region].keys():
                oligo_DB_seed[region][oligo_id] = {}
                start, end = (
                    database[region][oligo_id]["seed_region_start"],
                    database[region][oligo_id]["seed_region_end"],
                )
                seed_region_seq = database[region][oligo_id]["sequence"][
                    start : end + 1
                ]  # end must be included
                oligo_DB_seed[region][oligo_id]["sequence"] = seed_region_seq
        return oligo_DB_seed
