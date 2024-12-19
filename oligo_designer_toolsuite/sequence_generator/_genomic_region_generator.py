############################################
# imports
############################################

import os
import copy
import warnings
import random

import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None

from pathlib import Path

from Bio import SeqIO

from oligo_designer_toolsuite._constants import (
    SEPARATOR_FASTA_HEADER_FIELDS,
    SEPARATOR_FASTA_HEADER_FIELDS_LIST,
    SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS,
)
from oligo_designer_toolsuite.utils import GffParser
from oligo_designer_toolsuite.sequence_generator import FtpLoaderEnsembl, FtpLoaderNCBI

from ..utils._sequence_processor import (
    get_complement_regions,
    get_sequence_from_annotation,
)


############################################
# Genomic Region Generator Classes
############################################


class CustomGenomicRegionGenerator:

    def __init__(
        self,
        annotation_file: str,
        sequence_file: str,
        files_source: str = None,
        species: str = None,
        annotation_release: str = None,
        genome_assembly: str = None,
        dir_output: str = "output",
    ) -> None:
        """Constructor for the CustomGenomicRegionGenerator class."""
        if not files_source:
            files_source = "custom"
            warnings.warn(f"No source defined. Using default source {files_source}!")

        if not species:
            species = "unknown"
            warnings.warn(f"No species defined. Using default species {species}!")

        if not annotation_release:
            annotation_release = "unknown"
            warnings.warn(f"No annotation release defined. Using default release {annotation_release}!")

        if not genome_assembly:
            genome_assembly = "unknown"
            warnings.warn(f"No genome assembly defined. Using default genome assembly {genome_assembly}!")

        self.dir_output = os.path.join(dir_output, "annotation")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.files_source = files_source
        self.species = species
        self.annotation_release = annotation_release
        self.genome_assembly = genome_assembly
        self.annotation_file = annotation_file
        self.parsed_annotation_file = os.path.join(
            self.dir_output,
            os.path.basename(f"{'.'.join(annotation_file.split('.')[:-1])}.pckl"),
        )
        self.sequence_file = sequence_file

        # load annotation file and store in pickel file
        self.gff_parser = GffParser()
        self.gff_parser.check_gff_format(self.annotation_file)
        self.gff_parser.parse_annotation_from_gff(
            annotation_file=self.annotation_file,
            file_pickle=self.parsed_annotation_file,
        )

        # columns required for bed12 split sequence format
        self.BED_HEADER = ["seqid", "start", "end", "fasta_header", "score", "strand"]
        self.BED12_HEADER = [
            "seqid",
            "start",
            "end",
            "fasta_header",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "block_count",
            "block_sizes",
            "blockStarts",
        ]
        self.FILE_INFO = f"source-{self.files_source}_species-{self.species}_annotation_release-{self.annotation_release}_genome_assemly-{self.genome_assembly}"

    def get_sequence_gene(self) -> str:

        # get gene annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "gene")

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            f"source={self.files_source}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"species={self.species}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"annotation_release={self.annotation_release}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"genome_assembly={self.genome_assembly}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"regiontype=gene{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + "gene_id="
            + annotation["gene_id"].astype("str")
        )
        annotation["region"] = self._get_annotation_region(annotation)

        # add BED12 fields
        annotation["start"] = annotation["start_0base"]
        annotation["score"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["add_inf"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        # get sequence from bed file
        file_fasta = os.path.join(self.dir_output, f"gene_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_intergenic(self) -> str:
        def _get_chromosome_length() -> str:

            dict_chromosome_length = {}
            for rec in SeqIO.parse(self.sequence_file, "fasta"):
                dict_chromosome_length[rec.id] = len(rec.seq)

            file_chromosome_length = os.path.join(self.dir_output, "annotation.genome")
            with open(file_chromosome_length, "w") as handle:
                for key, value in sorted(dict_chromosome_length.items()):
                    handle.write(f"{key}\t{value}\n")

            return file_chromosome_length

        def _compute_intergenic_annotation(annotation, file_chromosome_length) -> pd.DataFrame:

            intergenic_annotation = []

            annotation.rename(
                columns={"start_0base": "start", "gene_id": "fasta_header"},
                inplace=True,
            )
            annotation = annotation[self.BED_HEADER]

            for seqid, gene_annotation in annotation.groupby("seqid"):
                gene_annotation_plusstrand = gene_annotation[gene_annotation.strand == "+"]
                gene_annotation_minusstrand = gene_annotation[gene_annotation.strand == "-"]

                intergenic_annotation.append(
                    _compute_intergenic_annotation_strand(
                        seqid, gene_annotation_plusstrand, "+", file_chromosome_length
                    )
                )
                intergenic_annotation.append(
                    _compute_intergenic_annotation_strand(
                        seqid, gene_annotation_minusstrand, "-", file_chromosome_length
                    )
                )

            intergenic_annotation = pd.concat(intergenic_annotation, ignore_index=True)
            return intergenic_annotation

        def _compute_intergenic_annotation_strand(
            seqid, gene_annotatio, strand, file_chromosome_length
        ) -> pd.DataFrame:

            # case 1: no annotated genes on the respective chromosome and strand
            if gene_annotatio.empty:
                chromosome_length = pd.read_csv(
                    file_chromosome_length, sep="\t", comment="t", header=0, names=["seqid", "length"]
                )
                intergenic_annotation = pd.DataFrame(
                    {
                        "seqid": seqid,
                        "start_0base": 0,
                        "end": chromosome_length.length[chromosome_length.seqid == seqid],
                        "start_1base": 1,
                        "region_id": "InterRegPlus" + str(seqid) + "_1",
                        "score": ".",
                        "strand": strand,
                    }
                )
            # case 2: annotated genes on the respective chromosome and strand
            else:
                # define files
                file_bed_in = os.path.join(self.dir_output, "annotation_in.bed")
                file_bed_out = os.path.join(self.dir_output, "annotation_out.bed")

                # save the annotation as bed file
                gene_annotatio = gene_annotatio.sort_values(by="start")
                gene_annotatio.to_csv(file_bed_in, sep="\t", header=False, index=False)

                # get complementary regions
                get_complement_regions(file_bed_in, file_chromosome_length, file_bed_out)

                # load intergenic regions
                intergenic_annotation = pd.read_csv(
                    file_bed_out, sep="\t", comment="t", header=0, names=["seqid", "start_0base", "end"]
                )
                intergenic_annotation["start_1base"] = intergenic_annotation["start_0base"] + 1
                if strand == "+":
                    intergenic_annotation["region_id"] = (
                        "InterRegPlus" + str(seqid) + "_" + intergenic_annotation.index.astype("str")
                    )
                if strand == "-":
                    intergenic_annotation["region_id"] = (
                        "InterRegMinus" + str(seqid) + "_" + intergenic_annotation.index[::-1].astype("str")
                    )
                intergenic_annotation["score"] = "."
                intergenic_annotation["strand"] = strand

                os.remove(file_bed_in)
                os.remove(file_bed_out)

            return intergenic_annotation

        # save chromosome sizes as genome file
        file_chromosome_length = _get_chromosome_length()

        # get gene annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "gene")
        annotation = _compute_intergenic_annotation(annotation, file_chromosome_length)

        # generate region_id
        annotation["add_inf"] = (
            f"source={self.files_source}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"species={self.species}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"annotation_release={self.annotation_release}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"genome_assembly={self.genome_assembly}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"regiontype=intergenic"
        )
        annotation["region"] = self._get_annotation_region(annotation)

        # add BED12 fields
        annotation["start"] = annotation["start_0base"]
        annotation["score"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["add_inf"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        # get sequence from bed file
        file_fasta = os.path.join(self.dir_output, f"intergenic_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        os.remove(file_chromosome_length)
        del annotation

        return file_fasta

    def get_sequence_exon(self, collapse_duplicated_regions: bool = True) -> str:
        # get exon annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "exon")

        # add transcript counts for each gene
        number_total_transcripts = self._get_number_total_transcripts()
        annotation = pd.merge(annotation, number_total_transcripts, on="gene_id", how="left")

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "transcript_id="
            + annotation["transcript_id"].astype("str")
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS}exon_number="
            + annotation["exon_number"].astype("str")
        )
        annotation["region"] = self._get_annotation_region(annotation)

        # remove duplicated entries
        if collapse_duplicated_regions:
            annotation = self._collapse_duplicated_regions(annotation)

        # add BED12 fields
        annotation["start"] = annotation["start_0base"]
        annotation["score"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + f"source={self.files_source}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"species={self.species}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"annotation_release={self.annotation_release}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"genome_assembly={self.genome_assembly}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"regiontype=exon{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + "gene_id="
            + annotation["gene_id"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS_LIST
            + annotation["add_inf"]
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}number_total_transcripts="
            + annotation["transcript_count"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = os.path.join(self.dir_output, f"exon_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_intron(self, collapse_duplicated_regions: bool = True) -> str:

        def _compute_intron_annotation(annotation) -> pd.DataFrame:

            intron_list = []

            for transcript, transcript_annotation in annotation.groupby("transcript_id"):
                gene_id = transcript_annotation.iloc[0].gene_id
                seqid = transcript_annotation.iloc[0].seqid
                strand = transcript_annotation.iloc[0].strand

                num_exons = transcript_annotation.shape[0]
                transcript_annotation = transcript_annotation.sort_values(by="start_1base")

                for i, (start_0base, start_1base, end) in enumerate(
                    zip(
                        transcript_annotation["start_0base"],
                        transcript_annotation["start_1base"],
                        transcript_annotation["end"],
                    )
                ):
                    attributes = self._get_attributes(start_0base, start_1base, end)

                    if i == 0:
                        exon_upstream = attributes

                    else:
                        exon_downstream = attributes

                        intron_number = i if strand == "+" else num_exons - i
                        intron_list.append(
                            [
                                gene_id,
                                transcript,
                                f"intron_{intron_number}",
                                seqid,
                                strand,
                                int(exon_upstream.end),  # start 0-bse
                                int(exon_upstream.end + 1),  # start 1-base
                                int(exon_downstream.start_0base),  # end
                            ]
                        )
                        exon_upstream = attributes

            intron_annotation = pd.DataFrame(
                np.asarray(intron_list),
                columns=[
                    "gene_id",
                    "transcript_id",
                    "intron_number",
                    "seqid",
                    "strand",
                    "start_0base",
                    "start_1base",
                    "end",
                ],
            )

            return intron_annotation

        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "exon")
        annotation = _compute_intron_annotation(annotation)

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "transcript_id="
            + annotation["transcript_id"].astype("str")
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS}intron_number="
            + annotation["intron_number"].astype("str")
        )
        annotation["region"] = self._get_annotation_region(annotation)

        # remove duplicated entries
        if collapse_duplicated_regions:
            annotation = self._collapse_duplicated_regions(annotation)

        # add BED12 fields
        annotation["start"] = annotation["start_0base"]
        annotation["score"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + f"source={self.files_source}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"species={self.species}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"annotation_release={self.annotation_release}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"genome_assembly={self.genome_assembly}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"regiontype=intron{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + "gene_id="
            + annotation["gene_id"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS_LIST
            + annotation["add_inf"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        # get sequence from bed file
        file_fasta = os.path.join(self.dir_output, f"intron_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_CDS(self, collapse_duplicated_regions: bool = True) -> str:

        # get exon annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "CDS")

        # add transcript counts for each gene
        number_total_transcripts = self._get_number_total_transcripts()
        annotation = pd.merge(annotation, number_total_transcripts, on="gene_id", how="left")

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "transcript_id="
            + annotation["transcript_id"].astype("str")
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS}exon_number="
            + annotation["exon_number"].astype("str")
        )
        annotation["region"] = self._get_annotation_region(annotation)

        # remove duplicated entries
        if collapse_duplicated_regions:
            annotation = self._collapse_duplicated_regions(annotation)

        # add BED12 fields
        annotation["start"] = annotation["start_0base"]
        annotation["score"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + f"source={self.files_source}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"species={self.species}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"annotation_release={self.annotation_release}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"genome_assembly={self.genome_assembly}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"regiontype=CDS{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + "gene_id="
            + annotation["gene_id"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS_LIST
            + annotation["add_inf"]
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}number_total_transcripts="
            + annotation["transcript_count"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = os.path.join(self.dir_output, f"cds_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_UTR(
        self,
        five_prime: bool = True,
        three_prime: bool = True,
        collapse_duplicated_regions: bool = True,
    ) -> str:
        def _compute_UTR(annotation: pd.DataFrame) -> pd.DataFrame:
            utrs = []

            # get leftmost and rightmost CDS boundaries
            for transcript, transcript_annotation in annotation.groupby("transcript_id"):
                cds_start = transcript_annotation[transcript_annotation.type == "CDS"].start_1base.min()
                cds_end = transcript_annotation[transcript_annotation.type == "CDS"].end.max()

                # based on strand, set GFF record type
                if transcript_annotation.iloc[0].strand == "+":
                    UTR_left_type = "five_prime_UTR"
                    UTR_right_type = "three_prime_UTR"

                elif transcript_annotation.iloc[0].strand == "-":
                    UTR_left_type = "three_prime_UTR"
                    UTR_right_type = "five_prime_UTR"

                exons = transcript_annotation[transcript_annotation.type == "exon"]

                UTR_left = copy.deepcopy(exons)
                UTR_left = UTR_left[UTR_left.start_1base < cds_start]
                UTR_left.type = UTR_left_type
                UTR_left.end[UTR_left.end >= cds_start] = cds_start - 1
                utrs.append(UTR_left)

                UTR_right = copy.deepcopy(exons)
                UTR_right = UTR_right[UTR_right.end > cds_end]
                UTR_right.type = UTR_right_type
                UTR_right.start_1base[UTR_right.start_1base <= cds_end] = cds_end + 1
                UTR_right.start_0base[(UTR_right.start_0base + 1) <= cds_end] = cds_end
                utrs.append(UTR_right)

            utr_annotation = pd.concat(utrs, ignore_index=True)

            return utr_annotation

        annotation_exon = self._get_annotation_region_of_interest(self._load_annotation(), "exon")
        annotation_CDS = self._get_annotation_region_of_interest(self._load_annotation(), "CDS")

        transcripts_with_CDS = list(set(annotation_CDS.transcript_id))
        annotation_exon = annotation_exon[annotation_exon.transcript_id.isin(transcripts_with_CDS)]

        annotation = pd.concat([annotation_exon, annotation_CDS], ignore_index=True)
        annotation = _compute_UTR(annotation)

        if five_prime == False:
            annotation = annotation[annotation.type == "three_prime_UTR"]
        if three_prime == False:
            annotation = annotation[annotation.type == "five_prime_UTR"]

        # add transcript counts for each gene
        number_total_transcripts = self._get_number_total_transcripts()
        annotation = pd.merge(annotation, number_total_transcripts, on="gene_id", how="left")

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "gene_id="
            + annotation["gene_id"].astype("str")
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS}transcript_id="
            + annotation["transcript_id"].astype("str")
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS}exon_number="
            + annotation["exon_number"].astype("str")
        )
        annotation["region"] = self._get_annotation_region(annotation)

        # remove duplicated entries
        if collapse_duplicated_regions:
            annotation = self._collapse_duplicated_regions(annotation)

        # add BED12 fields
        annotation["start"] = annotation["start_0base"]
        annotation["score"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + f"source={self.files_source}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"species={self.species}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"annotation_release={self.annotation_release}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"genome_assembly={self.genome_assembly}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"regiontype="
            + annotation["type"]
            + SEPARATOR_FASTA_HEADER_FIELDS_LIST
            + annotation["add_inf"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = os.path.join(self.dir_output, f"utr_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_exon_exon_junction(
        self, block_size: int, collapse_duplicated_regions: bool = True
    ) -> str:
        def _compute_exon_exon_junction_annotation(annotation: pd.DataFrame, block_size: int) -> pd.DataFrame:

            junction_list = []

            for transcript, transcript_annotation in annotation.groupby("transcript_id"):
                gene_id = transcript_annotation.iloc[0].gene_id
                seqid = transcript_annotation.iloc[0].seqid
                strand = transcript_annotation.iloc[0].strand

                transcript_annotation = transcript_annotation.sort_values(by="start_1base")

                for i, (exon_number, start_0base, start_1base, end) in enumerate(
                    zip(
                        transcript_annotation["exon_number"],
                        transcript_annotation["start_0base"],
                        transcript_annotation["start_1base"],
                        transcript_annotation["end"],
                    )
                ):
                    attributes = self._get_attributes(
                        start_0base,
                        start_1base,
                        end,
                        exon_number=exon_number,
                        exon_size=(end - start_0base),
                    )

                    if i == 0:
                        exon_upstream = attributes
                        exons_small = []
                        regions_exons_small = ""

                    # if exon is not the last exon of transcript and shorter than oligo block_size but not the last exon -> create sequence with neighboring exons
                    elif ((i + 1) < transcript_annotation.shape[0]) & ((end - start_0base) < block_size):
                        exons_small.append(attributes)

                    else:
                        exon_downstream = attributes
                        # catch case that first or last exon < block_size
                        block_size_up = min(block_size, exon_upstream.exon_size)
                        block_size_down = min(block_size, exon_downstream.exon_size)
                        start_up = exon_upstream.end - block_size_up
                        end_down = exon_downstream.start_0base + block_size_down

                        if exons_small == []:
                            block_count = 2
                            block_size_length_entry = f"{block_size_up},{block_size_down}"
                            block_size_start_entry = f"{0},{exon_downstream.start_0base - start_up}"

                        # if we have exons that are smaller than the block size add block counts
                        else:
                            block_count = len(exons_small) + 2
                            block_size_length_entry = (
                                str(block_size_up)
                                + ","
                                + ",".join([str(attributes.exon_size) for attributes in exons_small])
                                + ","
                                + str(block_size_down)
                            )
                            block_size_start_entry = (
                                "0,"
                                + ",".join(
                                    [str(attributes.start_0base - start_up) for attributes in exons_small]
                                )
                                + ","
                                + str(exon_downstream.start_0base - start_up)
                            )
                            # return region in 1-base offset
                            regions_exons_small = SEPARATOR_FASTA_HEADER_FIELDS_LIST.join(
                                [
                                    f"{seqid}:{attributes.start_1base}-{attributes.end}({strand})"
                                    for attributes in exons_small
                                ]
                            )
                        # return region in 1-base offset
                        region_up = f"{seqid}:{start_up + 1}-{start_up+block_size_up}({strand})"
                        region_down = f"{seqid}:{(end_down-block_size_down) + 1}-{end_down}({strand})"
                        junction_list.append(
                            [
                                gene_id,
                                transcript,
                                f"{exon_upstream.exon_number}__JUNC__{exon_downstream.exon_number}",
                                seqid,
                                strand,
                                start_up,
                                end_down,
                                SEPARATOR_FASTA_HEADER_FIELDS_LIST.join(
                                    filter(
                                        None,
                                        [region_up, regions_exons_small, region_down],
                                    )
                                ),
                                block_count,
                                block_size_length_entry,
                                block_size_start_entry,
                            ]
                        )
                        exons_small = []
                        regions_exons_small = ""
                        exon_upstream = attributes

            junction_annotation = pd.DataFrame(
                junction_list,
                columns=[
                    "gene_id",
                    "transcript_id",
                    "exon_number",
                    "seqid",
                    "strand",
                    "start",
                    "end",
                    "region_junction",
                    "block_count",
                    "block_sizes",
                    "blockStarts",
                ],
            )

            return junction_annotation

        # get exon annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "exon")

        # compute exon junctions
        annotation = _compute_exon_exon_junction_annotation(annotation, block_size)

        # add transcript counts for each gene
        number_total_transcripts = self._get_number_total_transcripts()
        annotation = pd.merge(annotation, number_total_transcripts, on="gene_id", how="left")

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "transcript_id="
            + annotation["transcript_id"].astype("str")
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS}exon_number="
            + annotation["exon_number"].astype("str")
        )
        # generate regions -> taken from exon junction regions
        annotation["region"] = annotation["region_junction"]

        # remove duplicated entries
        if collapse_duplicated_regions:
            annotation = self._collapse_duplicated_regions(annotation)

        # add BED12 fields
        annotation["score"] = 0
        annotation["thickStart"] = annotation["start"]
        annotation["thickEnd"] = annotation["end"]
        annotation["itemRgb"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"]
            + SEPARATOR_FASTA_HEADER_FIELDS
            + f"source={self.files_source}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"species={self.species}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"annotation_release={self.annotation_release}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"genome_assembly={self.genome_assembly}{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + f"regiontype=exonexonjunction{SEPARATOR_FASTA_HEADER_FIELDS_LIST}"
            + "gene_id="
            + annotation["gene_id"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS_LIST
            + annotation["add_inf"]
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}number_total_transcripts="
            + annotation["transcript_count"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED12_HEADER]

        file_fasta = os.path.join(self.dir_output, f"exon_exon_junction_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta)

        del annotation

        return file_fasta

    def _load_annotation(self) -> pd.DataFrame:

        # read annotation file and store in dataframe
        annotation = self.gff_parser.load_annotation_from_pickle(self.parsed_annotation_file)

        # required to ensure that sorting is done correctly
        annotation.start = annotation.start.astype("int")
        annotation.end = annotation.end.astype("int")

        # add both annotations to dataframe: GFF 1-base offset and BED 0-base offset
        # since we read in a GFF file, the start coordinates are 1-base offset
        annotation.rename(columns={"start": "start_1base"}, inplace=True)
        annotation["start_0base"] = annotation.start_1base - 1

        return annotation

    def _get_annotation_region_of_interest(self, annotation: pd.DataFrame, region: str) -> pd.DataFrame:

        region_annotation = annotation.loc[annotation["type"] == region]
        region_annotation.reset_index(inplace=True, drop=True)
        return region_annotation

    def _get_attributes(
        self,
        start_0base: int = None,
        start_1base: int = None,
        end: int = None,
        gene_id: str = None,
        exon_number: int = None,
        exon_size: int = None,
    ) -> pd.Series:

        attributes = pd.Series(
            data={
                "start_0base": start_0base,
                "start_1base": start_1base,
                "end": end,
                "gene_id": gene_id,
                "exon_number": exon_number,
                "exon_size": exon_size,
            }
        )

        return attributes

    def _collapse_duplicated_regions(self, annotation: pd.DataFrame) -> pd.DataFrame:

        aggregate_function = {col: "first" for col in annotation.columns}
        aggregate_function["add_inf"] = SEPARATOR_FASTA_HEADER_FIELDS_LIST.join

        merged_annotation = annotation.groupby(annotation["region"]).agg(aggregate_function)
        merged_annotation.reset_index(inplace=True, drop=True)

        return merged_annotation

    def _get_annotation_region(self, annotation: pd.DataFrame) -> str:

        region = (
            annotation["seqid"].astype("str")
            + ":"
            + annotation["start_1base"].astype("str")
            + "-"
            + annotation["end"].astype("str")
            + "("
            + annotation["strand"].astype("str")
            + ")"
        )

        return region

    def _get_sequence_from_annotation(
        self,
        annotation: pd.DataFrame,
        file_fasta: str,
        split: bool = True,
        strand: bool = True,
    ) -> None:
        annotation = annotation.sort_values(by=["fasta_header"])
        annotation.reset_index(inplace=True, drop=True)

        # save the annotation as bed file
        id = random.randint(0, 10000000)
        file_bed = os.path.join(self.dir_output, f"annotation_{id}.bed")
        annotation.to_csv(file_bed, sep="\t", header=False, index=False)

        # create the fasta file
        get_sequence_from_annotation(
            file_bed,
            self.sequence_file,
            file_fasta,
            split=split,
            strand=strand,
            nameOnly=True,
        )
        os.remove(file_bed)

    def _get_number_total_transcripts(self) -> pd.DataFrame:

        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "transcript")
        number_total_transcripts = annotation["gene_id"].value_counts()

        number_total_transcripts = number_total_transcripts.reset_index()
        number_total_transcripts.columns = ["gene_id", "transcript_count"]

        return number_total_transcripts


class NcbiGenomicRegionGenerator(CustomGenomicRegionGenerator):
    def __init__(
        self,
        taxon: str = None,
        species: str = None,
        annotation_release: str = None,
        dir_output: str = "output",
    ) -> None:
        """Constructor for the NcbiGenomicRegionGenerator class."""
        files_source = "NCBI"
        if taxon is None:
            taxon = "vertebrate_mammalian"
            warnings.warn(f"No taxon defined. Using default taxon {taxon}!")

        if species is None:
            species = "Homo_sapiens"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(f"No annotation release defined. Using default release {annotation_release}!")

        self.dir_output = os.path.join(dir_output, "annotation")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ftp = FtpLoaderNCBI(self.dir_output, taxon, species, annotation_release)
        annotation_file, annotation_release, genome_assembly = ftp.download_files("gtf")
        sequence_file, _, _ = ftp.download_files("fasta")

        super().__init__(
            annotation_file,
            sequence_file,
            files_source,
            species,
            annotation_release,
            genome_assembly,
            dir_output,
        )


class EnsemblGenomicRegionGenerator(CustomGenomicRegionGenerator):
    def __init__(
        self,
        species: str = None,
        annotation_release: str = None,
        dir_output: str = "output",
    ) -> None:
        """Constructor for the EnsemblGenomicRegionGenerator class."""
        files_source = "Ensemble"
        if not species:
            species = "homo_sapiens"
            warnings.warn(f"No species defined. Using default species {species}!")

        if not annotation_release:
            annotation_release = "current"
            warnings.warn(f"No annotation release defined. Using default release {annotation_release}!")

        self.dir_output = os.path.join(dir_output, "annotation")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ftp = FtpLoaderEnsembl(self.dir_output, species, annotation_release)
        annotation_file, annotation_release, genome_assembly = ftp.download_files("gtf")
        sequence_file, _, _ = ftp.download_files("fasta")

        super().__init__(
            annotation_file,
            sequence_file,
            files_source,
            species,
            annotation_release,
            genome_assembly,
            dir_output,
        )
