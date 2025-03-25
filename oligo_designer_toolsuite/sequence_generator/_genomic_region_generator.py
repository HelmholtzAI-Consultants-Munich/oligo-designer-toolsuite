############################################
# imports
############################################

import copy
import os
import random
import warnings

import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None

from pathlib import Path

from Bio import SeqIO

from oligo_designer_toolsuite._constants import (
    SEPARATOR_FASTA_HEADER_FIELDS,
    SEPARATOR_FASTA_HEADER_FIELDS_LIST,
)
from oligo_designer_toolsuite.sequence_generator import FtpLoaderEnsembl, FtpLoaderNCBI
from oligo_designer_toolsuite.utils import GffParser

from ..utils._sequence_processor import (
    get_complement_regions,
    get_sequence_from_annotation,
)

############################################
# Genomic Region Generator Classes
############################################


class CustomGenomicRegionGenerator:
    """
    This class is designed to generate custom genomic regions based on provided annotation and sequence files.
    It supports defining species, annotation release and genome assembly.

    The generated sequences are saved as fasta file with region id, additional information and coordinates in header.
    The header of each sequence must start with '>' and contain the following information:
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
    where the region_id is compulsory and the other fileds are opional. Coordinated are saved in 1-base format.

    Output Format (per sequence):
    >{region_id}::{additional information}::{chromosome}:{start}-{end}({strand})
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    :param annotation_file: The path to the annotation file (e.g., GFF).
    :type annotation_file: str
    :param sequence_file: The path to the corresponding sequence file (e.g., FASTA).
    :type sequence_file: str
    :param files_source: The source of the files (e.g., Ensembl, NCBI), defaults to "custom".
    :type files_source: str, optional
    :param species: The species name related to the annotation and sequence files, defaults to "unknown".
    :type species: str, optional
    :param annotation_release: The annotation release version, defaults to "unknown".
    :type annotation_release: str, optional
    :param genome_assembly: The genome assembly version, defaults to "unknown".
    :type genome_assembly: str, optional
    :param dir_output: The directory path for storing output files, defaults to "output".
    :type dir_output: str, optional
    """

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
        """
        Generates gene sequences based on gene annotations. These sequences are then saved in a FASTA file.

        Output Format (per sequence):
        >{gene_id}::source={source};species={species};annotation_release={annotation_release};
        genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id}
        ::{chromosome}:{start}-{end}({strand})

        sequence

        :return: The path to the generated FASTA file containing the gene sequences.
        :rtype: str
        """
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
        """
        Generates intergenic sequences based on gene annotations.

        This function extracts intergenic sequences (regions between genes) for both the positive and negative strands of a chromosome from the gene annotation and chromosome length.
        These sequences are then saved in a FASTA file.

        Output Format (per sequence):
        >{intergenic_region_id}::source={source};species={species};annotation_release={annotation_release};
        genome_assembly={genome_assembly};regiontype={regiontype}::{chromosome}:{start}-{end}({strand})

        sequence

        :return: The path to the generated FASTA file containing the intergenic region sequences.
        :rtype: str
        """

        def _get_chromosome_length() -> str:
            """
            Calculates the length of each chromosome in the given FASTA file and writes the lengths to a file.

            :return: The path to the file containing chromosome lengths.
            :rtype: str
            """
            dict_chromosome_length = {}
            for rec in SeqIO.parse(self.sequence_file, "fasta"):
                dict_chromosome_length[rec.id] = len(rec.seq)

            file_chromosome_length = os.path.join(self.dir_output, "annotation.genome")
            with open(file_chromosome_length, "w") as handle:
                for key, value in sorted(dict_chromosome_length.items()):
                    handle.write(f"{key}\t{value}\n")

            return file_chromosome_length

        def _compute_intergenic_annotation(annotation, file_chromosome_length) -> pd.DataFrame:
            """
            Computes the intergenic regions based on gene annotations for each chromosome and for both positive and negative strands.
            It uses chromosome length data to determine the regions between annotated genes.

            :param annotation: DataFrame containing gene annotations.
            :type annotation: pd.DataFrame
            :param file_chromosome_length: Path to the file containing chromosome lengths for each chromosome.
            :type file_chromosome_length: str
            :return: DataFrame with intergenic annotations including chromosome ID, start, end, and other related fields.
            :rtype: pd.DataFrame
            """
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
                        seqid=seqid,
                        gene_annotatio=gene_annotation_plusstrand,
                        strand="+",
                        file_chromosome_length=file_chromosome_length,
                    )
                )
                intergenic_annotation.append(
                    _compute_intergenic_annotation_strand(
                        seqid=seqid,
                        gene_annotatio=gene_annotation_minusstrand,
                        strand="-",
                        file_chromosome_length=file_chromosome_length,
                    )
                )

            intergenic_annotation = pd.concat(intergenic_annotation, ignore_index=True)
            return intergenic_annotation

        def _compute_intergenic_annotation_strand(
            seqid, gene_annotatio, strand, file_chromosome_length
        ) -> pd.DataFrame:
            """
            Computes the intergenic regions for a given chromosome and strand based on gene annotations.

            This method handles two cases:
            1. When there are no annotated genes on the respective chromosome and strand, it generates an intergenic region spanning from the start to the end of the chromosome.
            2. When there are annotated genes, it calculates intergenic regions between these genes using BED format files to determine the gaps.

            :param seqid: Identifier of the chromosome or sequence.
            :type seqid: str
            :param gene_annotatio: DataFrame containing gene annotations for the chromosome and strand.
            :type gene_annotatio: pd.DataFrame
            :param strand: The strand of interest ('+' or '-') to compute intergenic regions.
            :type strand: str
            :param file_chromosome_length: Path to the file containing chromosome lengths for each chromosome.
            :type file_chromosome_length: str
            :return: DataFrame with intergenic annotations including chromosome ID, start, end, and additional fields such as region_id and strand.
            :rtype: pd.DataFrame
            """

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
        """
        Generates exon sequences based on exon annotations.

        It optionally collapses duplicated regions, originating from different transcripts with the same start and end coordinates.
        It includes additional information in the FASTA header, such as transcript IDs and exon numbers.
        These sequences are then saved in a FASTA file.

        Output Format (per sequence):
        >{gene_id}::source={source};species={species};annotation_release={annotation_release};
        genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
        exon_number={exon_number_x};transcript_id={transcript_id_b},exon_number={exon_number_y};
        number_total_transcripts={number_total_transcripts}::{chromosome}:{start}-{end}({strand})
        sequence

        :param collapse_duplicated_regions: Whether to collapse duplicated regions into a single entry, defauls to True.
        :type collapse_duplicated_regions: bool
        :return: The path to the generated FASTA file containing the exon sequences.
        :rtype: str
        """

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
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}exon_number="
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
        """
        Generates intron sequences based on exon annotations.

        It optionally collapses duplicated regions, originating from different transcripts with the same start and end coordinates.
        It includes additional information in the FASTA header, such as transcript IDs and intron numbers.
        These sequences are then saved in a FASTA file.

        Output Format (per sequence):
        >{gene_id}::source={source};species={species};annotation_release={annotation_release};
        genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
        intron_number={intron_number_x};transcript_id={transcript_id_b},intron_number={intron_number_y};
        number_total_transcripts={number_total_transcripts}::{chromosome}:{start}-{end}({strand})
        sequence

        :param collapse_duplicated_regions: Whether to collapse duplicated regions into a single entry, defauls to True.
        :type collapse_duplicated_regions: bool
        :return: The path to the generated FASTA file containing the intron sequences.
        :rtype: str
        """

        def _compute_intron_annotation(annotation) -> pd.DataFrame:
            """
            Computes intron annotations based on exon annotations for each transcript.

            This method calculates intron regions by identifying the gaps between consecutive exons within each transcript.
            The introns are classified by their position relative to the exons and strand orientation.
            It generates a DataFrame containing intron annotations with details such as gene ID, transcript ID, and intron number.

            :param annotation: DataFrame containing exon annotations.
            :type annotation: pd.DataFrame
            :return: DataFrame with computed intron annotations.
            :rtype: pd.DataFrame
            """

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
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}intron_number="
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
        """
        Generates Coding Sequence (CDS) sequences based on CDS annotations.

        It optionally collapses duplicated regions, originating from different transcripts with the same start and end coordinates.
        It includes additional information in the FASTA header, such as transcript IDs and exon numbers.
        These sequences are then saved in a FASTA file.

        Output Format (per sequence):
        >{gene_id}::source={source};species={species};annotation_release={annotation_release};
        genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
        exon_number={exon_number_x};transcript_id={transcript_id_b},exon_number={exon_number_y};
        number_total_transcripts={number_total_transcripts}::{chromosome}:{start}-{end}({strand})
        sequence

        :param collapse_duplicated_regions: Whether to collapse duplicated regions into a single entry, defauls to True.
        :type collapse_duplicated_regions: bool
        :return: The path to the generated FASTA file containing the CDS sequences.
        :rtype: str
        """
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
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}exon_number="
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
        """
        Generates UTR (Untranslated Region) sequences based on exon and CDS (Coding Sequence) annotations.

        This method calculates both 5' and 3' UTR regions for transcripts based on the provided annotations.
        It allows for selective generation of 5' and/or 3' UTR sequences by specifying flags.
        It optionally collapses duplicated regions, originating from different transcripts with the same start and end coordinates.
        It includes additional information in the FASTA header, such as transcript IDs and exon numbers.
        These sequences are then saved in a FASTA file.

        Output Format (per sequence):
        >{gene_id}::source={source};species={species};annotation_release={annotation_release};
        genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
        exon_number={exon_number_x};transcript_id={transcript_id_b},exon_number={exon_number_y};
        number_total_transcripts={number_total_transcripts}::{chromosome}:{start}-{end}({strand})
        sequence

        :param five_prime: Boolean flag indicating whether to include 5' UTR sequences. Defaults to True.
        :type five_prime: bool
        :param three_prime: Boolean flag indicating whether to include 3' UTR sequences. Defaults to True.
        :type three_prime: bool
        :param collapse_duplicated_regions: Whether to collapse duplicated regions into a single entry, defauls to True.
        :type collapse_duplicated_regions: bool
        :return: The path to the generated FASTA file containing the UTR sequences.
        :rtype: str
        """

        def _compute_UTR(annotation: pd.DataFrame) -> pd.DataFrame:
            """
            Computes the 5' and 3' UTR (Untranslated Region) annotations from a DataFrame of exon and CDS annotations.

            This function processes exon and CDS annotations to determine UTR regions based on the boundaries of CDS (Coding Sequence) annotations.
            It distinguishes between 5' and 3' UTRs depending on the strand orientation and generates the corresponding annotations.

            :param annotation: DataFrame containing exon and CDS annotations.
            :type annotation: pd.DataFrame
            :return: DataFrame with computed UTR annotations, including 5' and/or 3' UTRs.
            :rtype: pd.DataFrame
            """
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
                UTR_left.loc[UTR_left["end"] >= cds_start, "end"] = cds_start - 1
                utrs.append(UTR_left)

                UTR_right = copy.deepcopy(exons)
                UTR_right = UTR_right[UTR_right.end > cds_end]
                UTR_right.type = UTR_right_type
                UTR_right.loc[UTR_right["start_1base"] <= cds_end, "start_1base"] = cds_end + 1
                UTR_right.loc[(UTR_right["start_0base"] + 1) <= cds_end, "start_0base"] = cds_end
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
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}transcript_id="
            + annotation["transcript_id"].astype("str")
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}exon_number="
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
        """
        Generates exon-exon junction sequences based on exon annotations.

        This function identifies junctions between consecutive exons within transcripts, considering specified block sizes to ensure that the generated sequences are of sufficient length.
        It computes the junction annotations and sequences based on exon boundaries.
        It optionally collapses duplicated regions, originating from different transcripts with the same start and end coordinates.
        It includes additional information in the FASTA header, such as transcript IDs and exon numbers.
        These sequences are then saved in a FASTA file.

        Output Format (per sequence):
        >{gene_id}::source={source};species={species};annotation_release={annotation_release};
        genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
        exon_number={exon_number_x__JUNC__exon_number_y};transcript_id={transcript_id_b},
        exon_number={exon_number_y__JUNC__exon_number_z};number_total_transcripts={number_total_transcripts}
        ::{chromosome}:{start}-{end}({strand})
        sequence

        :param block_size: The size of the sequence block used to generate junctions between exons. This parameter determines the length of the sequence blocks flanking the junctions.
        :type block_size: int
        :param collapse_duplicated_regions: Whether to collapse duplicated regions into a single entry, defauls to True.
        :type collapse_duplicated_regions: bool
        :return: The path to the generated FASTA file containing the UTR sequences.
        :rtype: str
        """

        def _compute_exon_exon_junction_annotation(annotation: pd.DataFrame, block_size: int) -> pd.DataFrame:
            """
            Computes exon-exon junction annotations from exon data.

            This function identifies junctions between consecutive exons within transcripts. It handles cases where exons are shorter than a specified block size by
            incorporating neighboring exons into the junction annotation, ensuring that the generated sequences are of sufficient length.

            :param annotation: A DataFrame containing exon annotations.
            :type annotation: pd.DataFrame
            :param block_size: The length of the sequence blocks flanking the junctions, i.e. +/- "block_size" bp around the junction.
            :type block_size: int
            :return: A DataFrame with annotations for exon-exon junctions.
            :rtype: pd.DataFrame
            """
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
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}exon_number="
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
        """
        Loads annotation data from a file and prepares it for processing.

        This function reads a GFF annotation file, loads it into a DataFrame, and processes the coordinates to include both 1-based and 0-based offsets.
        The resulting DataFrame includes columns for both coordinate systems, which are necessary for various genomic analyses.

        :return: A DataFrame containing the loaded and processed annotation data with columns for 1-based and 0-based start coordinates.
        :rtype: pd.DataFrame
        """
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
        """
        Filters the provided annotation DataFrame to retain only the rows corresponding to a specific region type.

        :param annotation: A DataFrame containing genomic annotations. It must include a 'type' column to specify the type of region.
        :type annotation: pd.DataFrame
        :param region: The type of region to filter for in the annotation DataFrame. For example, "exon" or "CDS".
        :type region: str
        :return: A DataFrame containing only the rows where the 'type' column matches the specified region.
        :rtype: pd.DataFrame
        """
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
        """
        Generates a pandas Series containing attributes related to genomic features.

        :param start_0base: The start coordinate of the feature in 0-based indexing.
        :type start_0base: int, optional
        :param start_1base: The start coordinate of the feature in 1-based indexing.
        :type start_1base: int, optional
        :param end: The end coordinate of the feature.
        :type end: int, optional
        :param gene_id: The identifier for the gene associated with the feature.
        :type gene_id: str, optional
        :param exon_number: The number of the exon within the transcript.
        :type exon_number: int, optional
        :param exon_size: The size of the exon.
        :type exon_size: int, optional
        :return: A pandas Series containing the provided attributes.
        :rtype: pd.Series
        """
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
        """
        Merges duplicate genomic regions in the annotation DataFrame by aggregating the information for each region.

        :param annotation: A DataFrame containing genomic annotations with possible duplicate regions. Each row represents a genomic feature, and there must be a 'region' column to identify the regions to be collapsed.
        :type annotation: pd.DataFrame
        :return: A DataFrame with duplicate regions merged. The values in columns other than 'region' are aggregated, with a special aggregation for the 'add_inf' column.
        :rtype: pd.DataFrame
        """
        aggregate_function = {col: "first" for col in annotation.columns}
        aggregate_function["add_inf"] = SEPARATOR_FASTA_HEADER_FIELDS_LIST.join

        merged_annotation = annotation.groupby(annotation["region"]).agg(aggregate_function)
        merged_annotation.reset_index(inplace=True, drop=True)

        return merged_annotation

    def _get_annotation_region(self, annotation: pd.DataFrame) -> str:
        """
        Generates a formatted string representing the genomic region for each annotation entry.

        :param annotation: A DataFrame containing genomic annotations. It must include columns 'seqid', 'start_1base', 'end', and 'strand' to construct the region string.
        :type annotation: pd.DataFrame
        :return: A string formatted as 'seqid:start_1base-end(strand)' for each annotation entry.
        :rtype: str
        """
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
        """
        Generates a FASTA file from genomic annotations by first converting the annotations to a BED file and then extracting sequences.

        :param annotation: A DataFrame containing genomic annotations. Must include columns used to generate BED files and FASTA sequences.
        :type annotation: pd.DataFrame
        :param file_fasta: The path to the output FASTA file where the sequences will be saved.
        :type file_fasta: str
        :param split: A flag indicating hether to split the sequences. Given BED12 input, extract and concatenate the sequences from the BED “blocks” (e.g., exons), defaults to True.
        :type split: bool
        :param strand: A flag indicating whether to force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented, defaults to True.
        :type strand: bool
        :return: None
        """
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
        """
        Calculates the total number of transcripts per gene from the annotation data.

        :return: A DataFrame where each row represents a gene with the total count of its transcripts.
        :rtype: pd.DataFrame
        """
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "transcript")
        number_total_transcripts = annotation["gene_id"].value_counts()

        number_total_transcripts = number_total_transcripts.reset_index()
        number_total_transcripts.columns = ["gene_id", "transcript_count"]

        return number_total_transcripts


class NcbiGenomicRegionGenerator(CustomGenomicRegionGenerator):
    """
    This class generates custom genomic regions using data from NCBI.
    It automatically downloads and processes annotation and sequence files based on the provided taxon, species, and annotation release.

    :param taxon: The taxonomic classification of the species, used to locate the appropriate NCBI files, defaults to "vertebrate_mammalian".
    :type taxon: str, optional
    :param species: The species name for which genomic regions will be generated, defaults to "Homo_sapiens".
    :type species: str, optional
    :param annotation_release: The version of the annotation release to use, defaults to "current".
    :type annotation_release: str, optional
    :param dir_output: The directory where output files will be stored, defaults to "output".
    :type dir_output: str, optional
    """

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
    """
    This class generates custom genomic regions using data from Ensembl.
    It automates the process of downloading and processing annotation and sequence files for the specified species and annotation release.

    :param species: The species name for which genomic regions will be generated, defaults to "homo_sapiens".
    :type species: str, optional
    :param annotation_release: The version of the annotation release to use, defaults to "current".
    :type annotation_release: str, optional
    :param dir_output: The directory where output files will be stored, defaults to "output".
    :type dir_output: str, optional
    """

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
