############################################
# imports
############################################

import os

import warnings
import subprocess
import numpy as np
import pandas as pd

from abc import abstractmethod
from typing import List, Union, Tuple
from Bio import SeqIO

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoAttributes, OligoDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import AlignmentSpecificityFilter

from ..utils._checkers_and_helpers import check_if_list
from ..utils._sequence_processor import get_sequence_from_annotation

############################################
# Oligo Blast Filter Classes
############################################


class BlastNFilter(AlignmentSpecificityFilter):
    def __init__(
        self,
        search_parameters: dict = {},
        hit_parameters: dict = {},
        names_search_output: list = [
            "query",
            "reference",
            "alignment_length",
            "query_start",
            "query_end",
            "query_length",
        ],
        filter_name: str = "blast_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BlastNFilter class."""
        super().__init__(filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters
        self.names_search_output = names_search_output

        # Define default output format for blast search filter. The fields are:
        # query, reference, alignment_length, query_start, query_end, query_length
        if "outfmt" not in self.search_parameters.keys():
            self.search_parameters["outfmt"] = "6 qseqid sseqid length qstart qend qlen"

    def _create_index(self, file_reference: str, n_jobs: int) -> str:
        ## Create blast index
        filename_reference_index = os.path.basename(file_reference)

        cmd = (
            "makeblastdb -dbtype nucl"
            + " -out "
            + os.path.join(self.dir_output, filename_reference_index)
            + " -in "
            + file_reference
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
            filename=f"oligo_database_blast_{region_name}",
            save_description=False,
            region_ids=region_ids,
            sequence_type=sequence_type,
        )
        file_blast_results = os.path.join(self.dir_output, f"blast_results_{region_name}.txt")

        cmd_parameters = ""
        for parameter, value in self.search_parameters.items():
            # add quotes if list of strings seperated by whitespace
            value = f'"{value}"' if " " in str(value) else value
            cmd_parameters += f" -{parameter} {value}"

        cmd = (
            "blastn"
            + " -query "
            + file_oligo_database
            + " -out "
            + file_blast_results
            + " -db "
            + os.path.join(self.dir_output, file_index)
            + " "
            + cmd_parameters
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output, stdout=subprocess.DEVNULL).wait()

        # read the reuslts of the blast seatch
        blast_results = self._read_search_output(
            file_search_results=file_blast_results,
            names_search_output=self.names_search_output,
        )

        # remove temporary files
        os.remove(file_oligo_database)
        os.remove(file_blast_results)

        # return loaded results
        return blast_results

    def _find_hits(
        self,
        oligo_database: OligoDatabase,  # not utilized in this filter
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
        region_ids: Union[str, List[str]],  # not utilized in this filter
    ) -> pd.DataFrame:
        if "min_alignment_length" in self.hit_parameters.keys():
            if "coverage" in self.hit_parameters.keys():
                warnings.warn(
                    "Both, 'min_alignment_length' and 'coverage' parameters were provided. Using 'min_alignment_length' parameter."
                )
            min_alignment_length = self.hit_parameters["min_alignment_length"]
        elif "coverage" in self.hit_parameters.keys():
            min_alignment_length = search_results["query_length"] * self.hit_parameters["coverage"] / 100
        else:
            raise KeyError("Please provide either 'coverage' or a 'min_alignment_length' in hit_parameters!")

        search_results["min_alignment_length"] = min_alignment_length

        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        blast_table_hits = search_results.loc[
            search_results.alignment_length > search_results.min_alignment_length
        ]

        return blast_table_hits

    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str) -> list:
        required_fields = [
            "query",
            "reference",
            "alignment_length",
            "query_start",
            "query_end",
            "query_length",
            "query_sequence",
            "reference_start",
            "reference_end",
            "reference_sequence",
            "reference_strand",
        ]
        if not all(field in self.names_search_output for field in required_fields):
            raise ValueError(
                f"Some of the required fields {required_fields} are missing in the search results."
            )
        # set the positions to a 0-based index
        table_hits = self._0_index_coordinates(table_hits)
        # Calculate adjusted "start" and "end" for BED format based on strand
        table_hits = self._extend_reference_start_end_coordinates(table_hits)

        # create bed file
        bed = pd.DataFrame(
            {
                "chr": table_hits["reference"],
                "start": table_hits["start"],
                "end": table_hits["end"],
                "name": table_hits["query"],
                "score": 0,
                "strand": table_hits["reference_strand"].map({"plus": "+", "minus": "-"}),
            }
        )

        # adjust for possible overflows (e.g. new coordinates are not included in the gene boundaries)
        # additionally we store how muchpadding we have to do to have two seqeunces of the same length
        bed = self._remove_overflows(bed, file_reference)
        file_bed = os.path.join(self.dir_output, f"references_{region_id}.bed")
        bed.to_csv(
            file_bed,
            sep="\t",
            index=False,
            header=False,
            columns=["chr", "start", "end", "name", "score", "strand"],
        )
        # generate the fasta file
        references_fasta_file = os.path.join(self.dir_output, f"references_{region_id}.fasta")
        get_sequence_from_annotation(
            file_bed, file_reference, references_fasta_file, strand=True, nameOnly=True
        )

        references = [off_reference.seq for off_reference in SeqIO.parse(references_fasta_file, "fasta")]
        references_padded = self._pad_overflows(bed, references)

        os.remove(references_fasta_file)
        os.remove(file_bed)
        return references_padded

    def _0_index_coordinates(self, table_hits: pd.DataFrame) -> pd.DataFrame:
        table_hits["query_start_corr"] = table_hits["query_start"] - 1
        table_hits["reference_start_corr"] = np.where(
            table_hits["reference_strand"] == "plus",
            table_hits["reference_start"] - 1,
            table_hits["reference_start"],
        )
        table_hits["reference_end_corr"] = np.where(
            table_hits["reference_strand"] == "minus",
            table_hits["reference_end"] - 1,
            table_hits["reference_end"],
        )
        return table_hits

    def _extend_reference_start_end_coordinates(self, table_hits: pd.DataFrame) -> pd.DataFrame:
        table_hits["start"] = np.where(
            table_hits["reference_strand"] == "plus",
            table_hits["reference_start_corr"] - table_hits["query_start_corr"],
            table_hits["reference_end_corr"] - (table_hits["query_length"] - table_hits["query_end"]),
        )
        table_hits["end"] = np.where(
            table_hits["reference_strand"] == "plus",
            table_hits["reference_end_corr"] + (table_hits["query_length"] - table_hits["query_end"]),
            table_hits["reference_start_corr"] + table_hits["query_start_corr"],
        )
        return table_hits

    def _remove_overflows(self, bed: pd.DataFrame, file_reference: str) -> pd.DataFrame:

        bed["overflow_start"] = bed["start"].apply(lambda x: -x if x < 0 else 0)
        bed["start"] = bed["start"].apply(lambda x: x if x >= 0 else 0)

        records = SeqIO.parse(file_reference, "fasta")
        regions_length = {record.id: len(record.seq) for record in records}
        bed["len_region"] = bed["chr"].map(regions_length)

        bed["overflow_end"] = bed[["end", "len_region"]].apply(
            lambda x: x["end"] - x["len_region"] if x["end"] > x["len_region"] else 0,
            axis=1,
        )
        bed["end"] = bed[["end", "len_region"]].apply(
            lambda x: x["end"] if x["end"] <= x["len_region"] else x["len_region"],
            axis=1,
        )
        return bed

    def _pad_overflows(self, bed: pd.DataFrame, references: list) -> list:
        references_padded = []
        for reference, overflow_start, overflow_end, strand in zip(
            references, bed["overflow_start"], bed["overflow_end"], bed["strand"]
        ):
            padding_start = "-" * overflow_start
            padding_end = "-" * overflow_end
            if strand == "+":
                reference = padding_start + reference + padding_end
            elif strand == "-":
                reference = padding_end + reference + padding_start

            references_padded.append(reference)
        return references_padded

    def _add_alignment_gaps(
        self, table_hits: pd.DataFrame, queries: list, references: list
    ) -> Tuple[list, list]:
        def add_gaps(seq, gaps):
            for gap in gaps:
                seq = seq[:gap] + "-" + seq[gap:]
            return seq

        table_hits["query_gaps"] = (
            table_hits["query_sequence"].apply(lambda x: np.where(np.array(list(x)) == "-")[0])
            + table_hits["query_start"]
            - 1  # blastn has 1-based indices
        )
        table_hits["reference_gaps"] = (
            table_hits["reference_sequence"].apply(lambda x: np.where(np.array(list(x)) == "-")[0])
            + table_hits["query_start"]
            - 1
        )
        gapped_queries = [add_gaps(query, gaps) for query, gaps in zip(queries, table_hits["query_gaps"])]
        gapped_references = [
            add_gaps(reference, gaps) for reference, gaps in zip(references, table_hits["reference_gaps"])
        ]
        return gapped_queries, gapped_references


############################################
# Oligo Blast Filter with Seedregion Classes
############################################


class BlastNSeedregionFilterBase(BlastNFilter):
    def __init__(
        self,
        search_parameters: dict = None,
        hit_parameters: dict = None,
        names_search_output: list = None,
        filter_name: str = "blast_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BlastNSeedregionFilterBase class."""
        if not search_parameters:
            search_parameters = {}
        if not hit_parameters:
            hit_parameters = {}
        if not names_search_output:
            names_search_output = [
                "query",
                "reference",
                "alignment_length",
                "query_start",
                "query_end",
                "query_length",
            ]
        super().__init__(search_parameters, hit_parameters, names_search_output, filter_name, dir_output)

    @abstractmethod
    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame, region_ids: Union[str, List[str]]
    ) -> pd.DataFrame:
        """ """

    def _find_hits(
        self,
        oligo_database: OligoDatabase,
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
        region_ids: Union[str, List[str]],
    ) -> pd.DataFrame:
        if "min_alignment_length" in self.hit_parameters.keys():
            if "coverage" in self.hit_parameters.keys():
                warnings.warn(
                    "Both, 'min_alignment_length' and 'coverage' parameters were provided. Using 'min_alignment_length' parameter."
                )
            min_alignment_length = self.hit_parameters["min_alignment_length"]
        elif "coverage" in self.hit_parameters.keys():
            min_alignment_length = search_results["query_length"] * self.hit_parameters["coverage"] / 100
        else:
            raise KeyError("Please provide either 'coverage' or a 'min_alignment_length' in hit_parameters!")

        search_results["min_alignment_length"] = min_alignment_length

        search_results = self._add_seed_region_information(
            oligo_database=oligo_database, search_results=search_results, region_ids=region_ids
        )

        # if seedregion not given
        search_results.seedregion_start.fillna(0, inplace=True)
        search_results.seedregion_end.fillna(search_results.query_length, inplace=True)

        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        blast_table_hits = search_results.loc[
            (search_results.alignment_length > search_results.min_alignment_length)
            | (
                (search_results.query_start < search_results.seedregion_start)
                & (search_results.query_end > search_results.seedregion_end)
            )
        ]

        return blast_table_hits


class BlastNSeedregionFilter(BlastNSeedregionFilterBase):
    def __init__(
        self,
        seedregion_start: Union[int, float],
        seedregion_end: Union[int, float],
        search_parameters: dict = None,
        hit_parameters: dict = None,
        names_search_output: list = None,
        filter_name: str = "blast_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BlastNSeedregionFilter class."""
        if not search_parameters:
            search_parameters = {}
        if not hit_parameters:
            hit_parameters = {}
        if not names_search_output:
            names_search_output = [
                "query",
                "reference",
                "alignment_length",
                "query_start",
                "query_end",
                "query_length",
            ]
        super().__init__(search_parameters, hit_parameters, names_search_output, filter_name, dir_output)

        self.seedregion_start = seedregion_start
        self.seedregion_end = seedregion_end

    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame, region_ids
    ) -> pd.DataFrame:
        oligo_attributes_calculator = OligoAttributes()
        oligo_database = oligo_attributes_calculator.calculate_seedregion(
            oligo_database=oligo_database,
            region_ids=region_ids,
            start=self.seedregion_start,
            end=self.seedregion_end,
        )

        seedregion = pd.merge(
            left=oligo_database.get_oligo_attribute_table("seedregion_start", region_ids),
            right=oligo_database.get_oligo_attribute_table("seedregion_end", region_ids),
            on="oligo_id",
        )
        search_results = pd.merge(
            left=search_results,
            right=seedregion,
            left_on="query",
            right_on="oligo_id",
            how="left",
        )
        return search_results


class BlastNSeedregionLigationsiteFilter(BlastNSeedregionFilterBase):
    def __init__(
        self,
        seedregion_size: int,
        search_parameters: dict = {},
        hit_parameters: dict = {},
        names_search_output: list = [
            "query",
            "reference",
            "alignment_length",
            "query_start",
            "query_end",
            "query_length",
        ],
        filter_name: str = "blast_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BlastNSeedregionLigationsiteFilter class."""
        super().__init__(search_parameters, hit_parameters, names_search_output, filter_name, dir_output)
        self.seedregion_size = seedregion_size

    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame, region_ids: Union[str, List[str]]
    ) -> pd.DataFrame:
        oligo_attributes_calculator = OligoAttributes()
        oligo_database = oligo_attributes_calculator.calculate_seedregion_ligationsite(
            oligo_database=oligo_database, region_ids=region_ids, seedregion_size=self.seedregion_size
        )

        seedregion = pd.merge(
            left=oligo_database.get_oligo_attribute_table("seedregion_start", region_ids),
            right=oligo_database.get_oligo_attribute_table("seedregion_end", region_ids),
            on="oligo_id",
        )
        search_results = pd.merge(
            left=search_results,
            right=seedregion,
            left_on="query",
            right_on="oligo_id",
            how="left",
        )
        return search_results
