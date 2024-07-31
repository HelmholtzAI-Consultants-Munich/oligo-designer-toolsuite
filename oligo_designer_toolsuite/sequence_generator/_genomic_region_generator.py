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
    """CustomGenomicRegionGenerator class.

    This class is designed for generating custom genomic regions based on annotation and sequence files. It provides
    functionality to parse GFF annotation files, store them in a pickle file, and generate custom genomic regions.

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

    :param annotation_file: Path to the GFF annotation file.
    :type annotation_file: str
    :param sequence_file: Path to the sequence file.
    :type sequence_file: str
    :param files_source: Source identifier for the files, defaults to "custom".
    :type files_source: str, optional
    :param species: Species identifier, defaults to "unknown".
    :type species: str, optional
    :param annotation_release: Annotation release version, defaults to "unknown".
    :type annotation_release: str, optional
    :param genome_assembly: Genome assembly version, defaults to "unknown".
    :type genome_assembly: str, optional
    :param dir_output: Output directory, defaults to "output".
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
    ):
        """Constructor for the CustomGenomicRegionGenerator class."""
        if files_source is None:
            files_source = "custom"
            warnings.warn(f"No source defined. Using default source {files_source}!")

        if species is None:
            species = "unknown"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "unknown"
            warnings.warn(f"No annotation release defined. Using default release {annotation_release}!")

        if genome_assembly is None:
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

    def get_sequence_gene(self):
        """Generate a sequence file for gene annotations.

        This method retrieves gene annotation entries, generates a unique region ID based on gene ID,
        and creates a BED file with the necessary fields. It then retrieves the sequence from the BED file and
        saves it as a FASTA file.

        Output Format (per sequence):
        >{gene_id}::{chromosome}:{start}-{end}({strand})
        sequence

        :return: Path to the generated FASTA file containing gene sequences.
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

    def get_sequence_intergenic(self):
        """Generate intergenic sequence annotations.

        This method computes intergenic annotations based on the provided gene annotations.
        By utilizing the "gene" attribute in the GTF annotation type field,
        it identifies gene regions and extracts intergenic regions from the entire chromosome.
        The sequence for these intergenic regions is then retrieved separately for both the plus and minus strands,
        respecting the strandness of the gene regions.

        Output Format (per sequence):
        >{intergenic_region_id}::{chromosome}:{start}-{end}({strand})
        sequence

        :return: Fasta file with gene sequences.
        :rtype: str
        """

        def _compute_intergenic_annotation(annotation):
            """Compute intergenic annotations from gene annotations.

            This function takes a DataFrame of gene annotations and computes intergenic annotations.
            It renames relevant columns for processing, separates annotations based on strand,
            and computes intergenic annotations for both the plus and minus strands.

            :param annotation: DataFrame containing gene annotations.
            :type annotation: pd.DataFrame
            :return: DataFrame containing intergenic annotations.
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
                    _compute_intergenic_annotation_strand(gene_annotation_plusstrand, "+")
                )
                intergenic_annotation.append(
                    _compute_intergenic_annotation_strand(gene_annotation_minusstrand, "-")
                )

            intergenic_annotation = pd.concat(intergenic_annotation, ignore_index=True)
            return intergenic_annotation

        def _compute_intergenic_annotation_strand(gene_annotatio, strand):
            """Compute intergenic annotations for a specific strand.

            This function takes a DataFrame of gene annotations for a specific strand and computes intergenic annotations.
            It saves the gene annotations and chromosome sizes as bed and genome files, respectively.
            It then calculates complementary regions and loads the resulting intergenic regions.
            The intergenic regions are formatted into a DataFrame with relevant columns.

            :param gene_annotation: DataFrame containing gene annotations for a specific strand.
            :type gene_annotation: pd.DataFrame
            :param strand: Strand of the gene annotations, either "+" or "-".
            :type strand: str
            :return: DataFrame containing intergenic annotations for the specified strand.
            :rtype: pd.DataFrame
            """
            # define files
            file_bed_in = os.path.join(self.dir_output, "annotation_in.bed")
            file_chromosome_length = os.path.join(self.dir_output, "annotation.genome")
            file_bed_out = os.path.join(self.dir_output, "annotation_out.bed")

            # save the annotation as bed file
            gene_annotatio = gene_annotatio.sort_values(by="start")
            gene_annotatio.to_csv(file_bed_in, sep="\t", header=False, index=False)

            # save chromosome sizes as genome file
            _get_chromosome_length(file_chromosome_length)

            # get complementary regions
            get_complement_regions(file_bed_in, file_chromosome_length, file_bed_out)

            # load intergenic regions
            intergenic_annotation = pd.read_csv(file_bed_out, sep="\t", comment="t", header=None)
            intergenic_annotation.columns = ["seqid", "start_0base", "end"]
            intergenic_annotation["start_1base"] = intergenic_annotation["start_0base"] + 1
            if strand == "+":
                intergenic_annotation["region_id"] = "InterRegPlus" + intergenic_annotation.index.astype(
                    "str"
                )
            if strand == "-":
                intergenic_annotation["region_id"] = "InterRegMinus" + intergenic_annotation.index[
                    ::-1
                ].astype("str")
            intergenic_annotation["score"] = "."
            intergenic_annotation["strand"] = strand

            os.remove(file_bed_in)
            os.remove(file_chromosome_length)
            os.remove(file_bed_out)

            return intergenic_annotation

        def _get_chromosome_length(file_chromosome_length):
            """Get the length of each chromosome and save it to a file.

            This function reads the sequence file in FASTA format to determine the length of each chromosome.
            The chromosome lengths are stored in a dictionary and written to the specified file.

            :param file_chromosome_length: Path to the file to store chromosome lengths.
            :type file_chromosome_length: str
            :return: Dictionary containing chromosome IDs and their respective lengths.
            :rtype: Dict[str, int]
            """
            dict_chromosome_length = {}
            for rec in SeqIO.parse(self.sequence_file, "fasta"):
                dict_chromosome_length[rec.id] = len(rec.seq)

            with open(file_chromosome_length, "w") as handle:
                for key, value in sorted(dict_chromosome_length.items()):
                    handle.write(f"{key}\t{value}\n")

            return dict_chromosome_length

        # get gene annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "gene")
        annotation = _compute_intergenic_annotation(annotation)

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

        del annotation

        return file_fasta

    def get_sequence_exon(self, collapse_duplicated_regions: bool = True):
        """Retrieve sequences for annotated exon regions.

        This method extracts exon annotation entries, generates unique region IDs, and retrieves the
        corresponding sequences for the exonic regions. The generated sequences are stored in a FASTA file.
        If `collapse_duplicated_regions` is set to True (default), then exons from different transcripts who share
        the exact same start and end coordinates, are merged into one sequence entry.

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},exon_number={exon_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param collapse_duplicated_regions: Flag to collapse duplicated regions, default is True.
        :type collapse_duplicated_regions: bool
        :return: Path to the generated FASTA file containing exon sequences.
        :rtype: str
        """
        # get exon annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "exon")

        # add transcript counts for each gene
        number_transcripts = self._get_number_transcripts()
        annotation = pd.merge(annotation, number_transcripts, on="gene_id", how="left")

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
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}number_transcripts="
            + annotation["transcript_count"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = os.path.join(self.dir_output, f"exon_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_intron(self, collapse_duplicated_regions: bool = True):
        """Generate and retrieve intron sequences based on the provided annotation.

        This method computes intron annotations from exon annotations and retrieves the intron sequences.
        The generated intron sequences are saved in a FASTA file. Intron are specified per transcript and
        calculated from the regions between exons. If `collapse_duplicated_regions` is set to True (default),
        then introns from different transcripts who share the exact same start and end coordinates,
        are merged into one sequence entry.

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},intron_number={intron_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param collapse_duplicated_regions: Flag to collapse duplicated regions in the annotation, defaults to True.
        :type collapse_duplicated_regions: bool, optional
        :return: Path to the generated FASTA file containing intron sequences.
        :rtype: str
        """

        def _compute_intron_annotation(annotation):
            """Compute intron annotations from transcript annotations.

            Given a DataFrame of transcript annotations, this method computes and returns intron annotations.
            The introns are determined based on the coordinates of exons in the transcripts.

            :param annotation: DataFrame containing transcript annotations.
            :type annotation: pd.DataFrame
            :return: DataFrame containing intron annotations.
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

    def get_sequence_CDS(self, collapse_duplicated_regions: bool = True):
        """Generate and retrieve CDS (coding sequence) sequences based on the provided annotation.

        This method retrieves CDS annotations, generates region IDs, and retrieves the CDS sequences.
        The generated CDS sequences are saved in a FASTA file. If `collapse_duplicated_regions` is set to True (default),
        then CDS from different transcripts who share the exact same start and end coordinates,
        are merged into one sequence entry.

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},exon_number={exon_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param collapse_duplicated_regions: Flag to collapse duplicated regions in the annotation, defaults to True.
        :type collapse_duplicated_regions: bool, optional
        :return: Path to the generated FASTA file containing CDS sequences.
        :rtype: str
        """
        # get exon annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "CDS")

        # add transcript counts for each gene
        number_transcripts = self._get_number_transcripts()
        annotation = pd.merge(annotation, number_transcripts, on="gene_id", how="left")

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
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}number_transcripts="
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
    ):
        """Generate sequence annotations for UTRs (untranslated regions).

        This method combines exon and CDS annotations, computes UTRs, and filters them based on the specified
        parameters (five_prime and three_prime). The UTR represents regions from the transcription start/end
        site or beginning of the known UTR to the base before the start codon of the transcript. If this region
        is interrupted  by introns then each exon or partial exon is annotated as a separate UTR feature.
        3 prime and 5 prime UTRs are defined dependent on the strand of the gene. If `collapse_duplicated_regions`
        is set to True (default), then UTRs from different transcripts who share the exact same start and end coordinates,
        are merged into one sequence entry. The resulting UTR annotations are stored in a BED file, and
        the corresponding sequence is retrieved and saved in a FASTA file.

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},exon_number={exon_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param five_prime: Include 5' UTRs in the output, defaults to True.
        :type five_prime: bool
        :param three_prime: Include 3' UTRs in the output, defaults to True.
        :type three_prime: bool
        :param collapse_duplicated_regions: Collapse duplicated regions, defaults to True.
        :type collapse_duplicated_regions: bool
        :return: Path to the generated FASTA file containing UTR sequences.
        :rtype: str
        """

        def _compute_UTR(annotation):
            """Compute UTR (untranslated region) annotations based on the provided annotation.

            This method calculates the UTRs for each transcript based on the CDS (coding sequence) boundaries.
            The resulting UTR annotations are returned as a DataFrame.

            :param annotation: DataFrame containing the transcript annotations.
            :type annotation: pd.DataFrame
            :return: DataFrame containing UTR annotations.
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

    def get_sequence_exon_exon_junction(self, block_size: int, collapse_duplicated_regions: bool = True):
        """Generate sequence annotations for exon-exon junctions based on the input exon annotation.

        This method computes exon-exon junction annotations for transcripts, considering the specified block size.
        To compute the exon junctions, the "block_size" parameter defines the size of the exon junction region,
        i.e. +/- "block_size" bp around the junction. If `collapse_duplicated_regions` is set to True (default),
        then junctions from  different transcripts who share the exact same start and end coordinates, are merged into one sequence entry.


        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},exon_number={exon_number}__JUNC__{exon_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param block_size: Size of the block for generating exon-exon junctions, i.e. +/- "block_size" bp around the junction.
        :type block_size: int
        :param collapse_duplicated_regions: Flag to indicate whether to collapse duplicated regions or not (default is True).
        :type collapse_duplicated_regions: bool, optional
        :return: File path of the generated exon-exon junction annotation in FASTA format.
        :rtype: str
        """

        def _compute_exon_exon_junction_annotation(annotation, block_size):
            """Compute exon-exon junction annotations based on the input exon annotation.

            This method generates junction annotations for exon-exon junctions in transcripts.
            The annotations include details such as gene_id, transcript_id, exon_number, sequence identifier (seqid),
            strand, start position, end position, region junction, block count, block sizes, and block starts.

            :param annotation: Input DataFrame containing exon annotations.
            :type annotation: pd.DataFrame
            :param block_size: Block size for generating junctions.
            :type block_size: int
            :return: DataFrame containing exon-exon junction annotations.
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
        number_transcripts = self._get_number_transcripts()
        annotation = pd.merge(annotation, number_transcripts, on="gene_id", how="left")

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
            + f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}number_transcripts="
            + annotation["transcript_count"].astype("str")
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED12_HEADER]

        file_fasta = os.path.join(self.dir_output, f"exon_exon_junction_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta)

        del annotation

        return file_fasta

    def _load_annotation(self):
        """Load annotation from the parsed GFF file into a DataFrame.

        This method reads the annotation file stored in a pickled format and returns a DataFrame.
        The DataFrame includes columns for gene ID, transcript ID, sequence identifier (seqid), source,
        type, start (0-base and 1-base offset), end, score, strand, and additional attributes.

        :return: DataFrame containing the loaded annotation.
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

    def _get_annotation_region_of_interest(self, annotation: pd.DataFrame, region: str):
        """Extract annotations for a specific genomic region type.

        Given an annotation DataFrame and a specified genomic region type (e.g., "exon" or "CDS"),
        this method filters the DataFrame to retain only entries corresponding to the specified region.

        :param annotation: DataFrame containing the genomic annotations.
        :type annotation: pd.DataFrame
        :param region: Genomic region type of interest (e.g., "exon" or "CDS").
        :type region: str
        :return: DataFrame containing annotations for the specified genomic region type.
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
    ):
        """Generate a Pandas Series containing attributes for a genomic feature.

        This method creates a Pandas Series with attributes such as start and end coordinates,
        gene ID, exon number, and exon size. These attributes are commonly used in genomic annotations.

        :param start_0base: Start coordinate of the feature in 0-based offset.
        :type start_0base: int, optional
        :param start_1base: Start coordinate of the feature in 1-based offset.
        :type start_1base: int, optional
        :param end: End coordinate of the feature.
        :type end: int, optional
        :param gene_id: Identifier of the gene to which the feature belongs.
        :type gene_id: str, optional
        :param exon_number: Number identifying the exon in the gene.
        :type exon_number: int, optional
        :param exon_size: Size of the exon.
        :type exon_size: int, optional
        :return: Pandas Series containing the specified attributes.
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

    def _collapse_duplicated_regions(self, annotation: pd.DataFrame):
        """Collapse duplicated regions in the annotation DataFrame.

        This method collapses duplicated regions in the annotation DataFrame based on the "region" column,
        merging entries with the same region into a single entry.

        :param annotation: Input DataFrame containing genomic annotation.
        :type annotation: pd.DataFrame
        :return: DataFrame with collapsed duplicated regions.
        :rtype: pd.DataFrame
        """
        aggregate_function = {col: "first" for col in annotation.columns}
        aggregate_function["add_inf"] = SEPARATOR_FASTA_HEADER_FIELDS_LIST.join

        merged_annotation = annotation.groupby(annotation["region"]).agg(aggregate_function)
        merged_annotation.reset_index(inplace=True, drop=True)

        return merged_annotation

    def _get_annotation_region(self, annotation: pd.DataFrame):
        """Generate a region string from the given annotation DataFrame.

        This method creates a region string based on the "seqid," "start_1base," "end," and "strand" columns
        in the input annotation DataFrame.

        :param annotation: Input DataFrame containing genomic annotation.
        :type annotation: pd.DataFrame
        :return: Series containing region strings.
        :rtype: pd.Series
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
    ):
        """Retrieve sequence from the annotation and save it in a FASTA file.

        This method sorts the annotation DataFrame by "fasta_header," saves it as a BED file,
        and then uses the BED file along with the reference sequence file to generate a FASTA file
        containing the annotated sequences.

        :param annotation: DataFrame containing annotation information.
        :type annotation: pd.DataFrame
        :param file_fasta: Path to the output FASTA file.
        :type file_fasta: str
        :param split: Flag to split the sequences, default is True.
        :type split: bool, optional
        :param strand: Flag to consider strand information, default is True.
        :type strand: bool, optional
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

    def _get_number_transcripts(self):
        """Get the number of transcripts associated with each gene.

        This function loads the annotation, extracts transcripts, and counts the number of transcripts
        for each gene. It returns a DataFrame with 'gene_id' as the first column and 'transcript_count'
        as the second column, where 'gene_id' contains unique gene IDs and 'transcript_count' contains
        the corresponding number of transcripts for each gene.

        :return: DataFrame with each gene ID and its associated transcript count.
        :rtype: pandas.DataFrame
        """
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "transcript")
        number_transcripts = annotation["gene_id"].value_counts()

        number_transcripts = number_transcripts.reset_index()
        number_transcripts.columns = ["gene_id", "transcript_count"]

        return number_transcripts


class NcbiGenomicRegionGenerator(CustomGenomicRegionGenerator):
    """A class for generating genomic regions using NCBI as the data source.

    This class provides functionality to download genomic annotation and sequence files
    from NCBI for a specified taxon, species, and annotation release. It inherits from
    CustomGenomicRegionGenerator and utilizes an FTP loader to fetch the necessary files.

    :param taxon: Taxonomic group for the genomic region, default is "vertebrate_mammalian".
    :type taxon: str, optional
    :param species: Species for the genomic region, default is "Homo_sapiens".
    :type species: str, optional
    :param annotation_release: Annotation release version, default is "current".
    :type annotation_release: str, optional
    :param dir_output: Directory path for output files, default is "output/annotation".
    :type dir_output: str, optional
    """

    def __init__(
        self,
        taxon: str = None,
        species: str = None,
        annotation_release: str = None,
        dir_output: str = "output",
    ):
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
    """A class for generating genomic regions using Ensembl as the data source.

    This class provides functionality to download genomic annotation and sequence files
    from Ensembl for a specified species and annotation release. It inherits from
    CustomGenomicRegionGenerator and utilizes an FTP loader to fetch the necessary files.

    :param species: Species for the genomic region, default is "homo_sapiens".
    :type species: str, optional
    :param annotation_release: Annotation release version, default is "current".
    :type annotation_release: str, optional
    :param dir_output: Directory path for output files, default is "output/annotation".
    :type dir_output: str, optional
    """

    def __init__(
        self,
        species: str = None,
        annotation_release: str = None,
        dir_output: str = "output",
    ):
        """Constructor for the EnsemblGenomicRegionGenerator class."""
        files_source = "Ensemble"
        if species is None:
            species = "homo_sapiens"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
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
