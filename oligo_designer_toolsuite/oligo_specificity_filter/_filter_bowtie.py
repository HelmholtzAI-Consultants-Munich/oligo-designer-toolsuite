############################################
# Imports
############################################

import os
import subprocess
from typing import List, Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import AlignmentSpecificityFilter

from ..utils._checkers_and_helpers import check_if_list
from ..utils._sequence_processor import get_sequence_from_annotation

############################################
# Oligo Bowtie Filter Classes
############################################


class BowtieFilter(AlignmentSpecificityFilter):
    def __init__(
        self,
        search_parameters: dict = {},
        hit_parameters: dict = {},
        names_search_output: list = [
            "query",
            "strand",
            "reference",
            "reference_start",
            "query_sequence",
            "read_quality",
            "num_instances",
            "mismatch_positions",
        ],
        filter_name: str = "bowtie_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BowtieFilter class."""
        super().__init__(filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters  # currently not used
        self.names_search_output = names_search_output

    def _create_index(self, file_reference: str, n_jobs: int) -> str:

        ## Create bowtie index
        file_reference = os.path.abspath(file_reference)
        filename_reference_index = os.path.basename(file_reference)

        cmd = (
            "bowtie-build --offrate 4"
            + " --threads "
            + str(n_jobs)
            + " -f "
            + file_reference
            + " "
            + filename_reference_index
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output, stdout=subprocess.DEVNULL).wait()
        return filename_reference_index

    def _run_search(
        self,
        oligo_database: OligoDatabase,
        file_index: str,
        sequence_type: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ) -> pd.DataFrame:
        region_ids = check_if_list(region_ids)
        if region_ids:
            region_name = "_".join(region_ids)
        else:
            region_name = "all_regions"

        file_oligo_database = oligo_database.write_database_to_fasta(
            filename=f"oligo_database_bowtie_{region_name}",
            save_description=False,
            region_ids=region_ids,
            sequence_type=sequence_type,
        )
        file_bowtie_results = os.path.join(self.dir_output, f"bowtie_results_{region_name}.txt")

        cmd_parameters = ""
        for parameter, value in self.search_parameters.items():
            cmd_parameters += f" {parameter} {value}"

        cmd = (
            "bowtie"
            + " -x "
            + file_index
            + " -f"  # fasta file is input
            + " -a"  # report all alignments -> TODO: does this make sense or set e.g. -k 100
            + cmd_parameters
            + " "
            + file_oligo_database
            + " "
            + file_bowtie_results
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output, stdout=subprocess.DEVNULL).wait()

        # read the reuslts of the bowtie search
        bowtie_results = self._read_search_output(
            file_search_results=file_bowtie_results,
            names_search_output=self.names_search_output,
        )

        # remove temporary files
        os.remove(file_oligo_database)
        os.remove(file_bowtie_results)

        # return loaded results
        return bowtie_results

    def _find_hits(
        self,
        oligo_database: OligoDatabase,  # not used in this filter
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
        region_ids: Union[str, List[str]],  # not used in this filter
    ) -> pd.DataFrame:
        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        return search_results

    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str) -> list:
        required_fields = [
            "query",
            "strand",
            "reference",
            "reference_start",
            "query_sequence",
        ]
        if not all(field in self.names_search_output for field in required_fields):
            raise ValueError(
                f"Some of the required fields {required_fields} are missing in the search results."
            )
        table_hits["reference_end"] = table_hits.apply(
            lambda x: x["reference_start"] + len(x["query_sequence"]), axis=1
        )
        bed = pd.DataFrame(
            {
                "chr": table_hits["reference"],
                "start": table_hits["reference_start"],
                "end": table_hits["reference_end"],
                "name": table_hits["query"],
                "score": 0,
                "strand": table_hits["strand"],
            }
        )
        file_bed = os.path.join(self.dir_output, f"references_{region_id}.bed")
        bed.to_csv(file_bed, sep="\t", index=False, header=False)

        references_fasta_file = os.path.join(self.dir_output, f"references_{region_id}.fasta")

        get_sequence_from_annotation(
            file_bed, file_reference, references_fasta_file, strand=True, nameOnly=True
        )
        references = [off_reference.seq for off_reference in SeqIO.parse(references_fasta_file, "fasta")]
        os.remove(references_fasta_file)
        os.remove(file_bed)
        return references

    def _add_alignment_gaps(
        self, table_hits: pd.DataFrame, queries: list, references: list
    ) -> tuple[list, list]:
        # bowtie does not support gaps
        return queries, references


############################################
# Oligo Bowtie2 Filter Classes
############################################


class Bowtie2Filter(AlignmentSpecificityFilter):
    def __init__(
        self,
        search_parameters: dict = {},
        hit_parameters: dict = {},
        names_search_output: list = [
            "query",
            "flags",
            "reference",
            "reference_start",
            "mapping_quality",
            "CIGAR_alignment",
            "mate_sequence_name",
            "mate_sequence_offset",
            "mate_sequence_fragment_length",
            "sequence",
            "read_qualities",
        ],
        filter_name: str = "bowtie2_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the Bowtie2Filter class."""
        super().__init__(filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters  # currently not used
        self.names_search_output = names_search_output

    def _create_index(self, file_reference: str, n_jobs: int) -> str:
        ## Create bowtie index
        file_reference = os.path.abspath(file_reference)
        filename_reference_index = os.path.basename(file_reference)

        # Check if bowtie database exists -> check for any of the bowtie index files, e.g. ".1.bt2" file
        cmd = (
            "bowtie2-build --quiet --offrate 4"
            + " --threads "
            + str(n_jobs)
            + " -f "
            + file_reference
            + " "
            + filename_reference_index
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output).wait()

        return filename_reference_index

    def _run_search(
        self,
        oligo_database: OligoDatabase,
        file_index: str,
        sequence_type: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ) -> pd.DataFrame:
        region_ids = check_if_list(obj=region_ids)
        if region_ids:
            region_name = "_".join(region_ids)
        else:
            region_name = "all_regions"

        file_oligo_database = oligo_database.write_database_to_fasta(
            filename=f"oligo_database_bowtie2_{region_name}",
            save_description=False,
            region_ids=region_ids,
            sequence_type=sequence_type,
        )
        file_bowtie_results = os.path.join(self.dir_output, f"bowtie2_results_{region_name}.txt")

        cmd_parameters = ""
        for parameter, value in self.search_parameters.items():
            cmd_parameters += f" {parameter} {value}"

        cmd = (
            "bowtie2 --quiet"
            + " --no-hd --no-unal"
            + " -x "
            + file_index
            + " -f"  # fast file is input
            + " -a"  # report all alignments -> TODO: does this make sense or set e.g. -k 100
            + cmd_parameters
            + " -U "
            + file_oligo_database
            + " -S "
            + file_bowtie_results
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output).wait()

        # read the reuslts of the bowtie seatch
        bowtie_results = self._read_search_output(
            file_search_results=file_bowtie_results,
            names_search_output=self.names_search_output,
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        )

        # remove temporary files
        os.remove(file_oligo_database)
        os.remove(file_bowtie_results)

        # return loaded results
        return bowtie_results

    def _find_hits(
        self,
        oligo_database: OligoDatabase,  # not used in this filter
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
        region_ids: Union[str, List[str]],  # not used in this filter
    ) -> pd.DataFrame:
        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        return search_results

    def _get_references(self, search_results: pd.DataFrame, file_reference: str, region_id: str) -> list:
        raise NotImplementedError("AI filters not supported for Bowtie2.")

    def _add_alignment_gaps(
        self, search_results: pd.DataFrame, queries: list, references: list
    ) -> tuple[list, list]:
        raise NotImplementedError("AI filters not supported for Bowtie2.")
