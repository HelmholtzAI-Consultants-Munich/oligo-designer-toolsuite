import inspect
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

# from typing_extensions import Literal # Python 3.7 or below
from typing import Literal

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.sequence_generator import (
    OligoSequenceGenerator,
)
from oligo_designer_toolsuite.pipelines._utils import log_parameters


class BaseOligoDesigner:
    """ """

    def __init__(
        self,
        dir_output: str = "output",
        log_name: str = "oligo_designer",
        write_removed_genes: bool = True,
        write_intermediate_steps: bool = True,
    ):
        """Constructor"""
        ##### store parameters #####
        self.dir_output = dir_output
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.write_removed_genes = write_removed_genes
        self.write_intermediate_steps = write_intermediate_steps
        self.log_name = log_name

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_{log_name}_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("log_name")
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.NOTSET,
            handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
        )
        logging.captureWarnings(True)

        ##### log parameters #####
        logging.info("Parameters Init:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        log_parameters(parameters)

    def _get_oligo_database_info(self, oligo_database: dict):
        """Count the number of oligos and genes in the database.

        :param oligo_database: Database with oligos.
        :type oligo_database: dict
        :return: Number of genes and oligos in the database.
        :rtype: int, int
        """
        genes = oligo_database.keys()
        num_genes = len(genes)
        num_oligos = 0
        for gene in genes:
            num_oligos += len(oligo_database[gene].keys())

        return num_genes, num_oligos

    def _get_oligo_length_min_max_from_database(self, oligo_database: dict):
        """Get minimum and maximum length of oligos stored in the oligo database.

        :param oligo_database: Database with oligos.
        :type oligo_database: dict
        :return: Min and max length of oligos
        :rtype: int, int
        """
        oligo_length_min = sys.maxsize
        oligo_length_max = 0

        for region in oligo_database.keys():
            for oligo in oligo_database[region].keys():
                length = oligo_database[region][oligo]["length"]
                if length < oligo_length_min:
                    oligo_length_min = length
                if length > oligo_length_max:
                    oligo_length_max = length

        return oligo_length_min, oligo_length_max

    

    def create_oligo_database(
        self,
        regions: list,
        oligo_length_min: int,
        oligo_length_max: int,
        files_fasta_oligo_database: list[str],
        isoform_consensus: Literal["intersection", "union"] = "union",
        min_oligos_per_region: int = 0,
        n_jobs: int = 1,
    ):
        ##### log parameters #####
        logging.info("Parameters Create Database:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        log_parameters(parameters)

        if isoform_consensus == "intersection":  # TODO: what does it mean??
            raise Exception(f"Isoform consensus: {isoform_consensus} not implemented yet.")

        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            filename_out=f"{self.log_name}_oligos",
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            region_ids=regions,
            n_jobs=n_jobs,
        )

        ##### creating the oligo database #####
        # oligo database
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_region,
            write_regions_with_insufficient_oligos=True,
            dir_output=self.dir_output,
        )
        # load the oligo sequences
        oligo_database.load_sequences_from_fasta(
            files_fasta=[oligo_fasta_file],
            sequence_type="oligo",
            region_ids=regions,
        )

        ##### loggig database information #####
        if self.write_removed_genes:
            logging.info(
                f"Genes with <= {min_oligos_per_region} oligos will be removed from the oligo database and their names will be stored in '{oligo_database.file_removed_regions}'."
            )

        num_genes, num_oligos = self._get_oligo_database_info(oligo_database.database)
        logging.info(
            f"Step - Generate oligos: the database contains {num_oligos} oligos from {num_genes} genes."
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="oligo_database_initial.txt")
        else:
            file_database = ""

        return oligo_database, file_database

    def load_oligo_database(
        self,
        file_database: str,
        min_oligos_per_region: int = 0,
    ):
        ##### log parameters #####
        logging.info("Parameters Load Database:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        log_parameters(parameters)

        ##### loading the oligo database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_region,
            dir_output=self.dir_output,
        )
        oligo_database.load_database(file_database)

        ##### loggig database information #####
        if self.write_removed_genes:
            logging.info(
                f"Genes with <= {min_oligos_per_region} oligos will be removed from the oligo database and their names will be stored in '{oligo_database.file_removed_regions}'."
            )

        num_regions, num_oligos = self._get_oligo_database_info(oligo_database.database)
        logging.info(
            f"Step - Generate oligos: the database contains {num_oligos} oligos from {num_regions} genes."
        )

        return oligo_database
    
    