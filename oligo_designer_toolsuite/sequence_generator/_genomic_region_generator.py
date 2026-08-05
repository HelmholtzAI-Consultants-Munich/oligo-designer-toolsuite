"""
Build FASTA sequences for annotated genomic regions from GFF and genome files.

Given a genome annotation and the matching FASTA, this module extracts intervals
such as genes, exons, introns, CDS, UTRs, intergenic stretches, and exon-exon junctions.
Use :class:`CustomGenomicRegionGenerator` with local files, or
:class:`NcbiGenomicRegionGenerator` / :class:`EnsemblGenomicRegionGenerator` to download
annotation and sequence from NCBI or Ensembl first.
"""

############################################
# imports
############################################

import copy
import os
import random
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from oligo_designer_toolsuite.utils._checkers_and_helpers import safe_append_filename

pd.options.mode.chained_assignment = None

from pathlib import Path

from Bio import SeqIO

from oligo_designer_toolsuite._constants import (
    SEPARATOR_FASTA_HEADER_FIELDS,
    SEPARATOR_FASTA_HEADER_FIELDS_LIST,
)
from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.sequence_generator import FtpLoaderEnsembl, FtpLoaderNCBI
from oligo_designer_toolsuite.utils import (
    BedParser,
    GffParser,
    get_complement_regions,
    get_sequence_from_annotation,
    logger,
)

############################################
# Genomic Region Generator Classes
############################################


class CustomGenomicRegionGenerator:
    """
    Build FASTA sequences for custom genomic regions from local annotation and genome files.

    Provide a GFF annotation and matching genome FASTA, plus optional species and assembly
    labels. Region extractors (gene, exon, intron, CDS, UTR, intergenic, exon-exon junction)
    write one FASTA file per region type.

    Coordinate indexing:

    Annotation files (GFF/GTF) use 1-based starts. Sequence extraction writes temporary
    BED files, which use 0-based starts. FASTA headers written by this class also use
    1-based starts so they stay aligned with the annotation.

    Each sequence header starts with ``>`` and holds region ID, optional metadata, and
    coordinates (chromosome, start, end, strand). The region ID is required; other
    fields are optional.

    Output Format (per sequence):

        >{region_id}::{additional information}::{chromosome}:{start}-{end}({strand})
        sequence

    Example::

        >ASR1::transcript_id=XM456,exon_number=5::16:54552-54786(+)
        AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    :param annotation_file: Path to the annotation file (for example GFF/GTF).
    :type annotation_file: str
    :param sequence_file: Path to the matching genome FASTA.
    :type sequence_file: str
    :param files_source: Label for the file source (for example Ensembl or NCBI). Defaults to ``custom``.
    :type files_source: str, optional
    :param species: Species name for the annotation and genome. Defaults to ``unknown``.
    :type species: str, optional
    :param annotation_release: Annotation release version. Defaults to ``unknown``.
    :type annotation_release: str, optional
    :param genome_assembly: Genome assembly version. Defaults to ``unknown``.
    :type genome_assembly: str, optional
    :param dir_output: Directory for output files. Defaults to ``output``.
    :type dir_output: str, optional
    """

    def __init__(
        self,
        annotation_file: str,
        sequence_file: str,
        files_source: str | None = None,
        species: str | None = None,
        annotation_release: str | None = None,
        genome_assembly: str | None = None,
        dir_output: str = "output",
    ) -> None:
        """Constructor for the CustomGenomicRegionGenerator class."""
        if files_source is None:
            files_source = "custom"
            logger.warning(f"No source defined. Using default source {files_source}!")

        if species is None:
            species = "unknown"
            logger.warning(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "unknown"
            logger.warning(f"No annotation release defined. Using default release {annotation_release}!")

        if genome_assembly is None:
            genome_assembly = "unknown"
            logger.warning(f"No genome assembly defined. Using default genome assembly {genome_assembly}!")

        self.dir_output = os.path.join(dir_output, "annotation")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.files_source = files_source
        self.species = species
        self.annotation_release = annotation_release
        self.genome_assembly = genome_assembly
        self.annotation_file = annotation_file
        self.parsed_annotation_file = safe_append_filename(
            self.dir_output,
            os.path.basename(f"{'.'.join(annotation_file.split('.')[:-1])}.pckl"),
        )
        self.sequence_file = sequence_file

        # load annotation file and store in pickel file
        self.gff_parser = GffParser()
        self.bed_parser = BedParser()
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
        self.FILE_INFO = f"source__{self.files_source}__species__{self.species}__annotation_release__{self.annotation_release}__genome_assemly__{self.genome_assembly}"

    def get_sequence_gene(self) -> str:
        """
        Build FASTA sequences for gene intervals from gene annotations.

        Loads gene features from the annotation, labels each with source and gene ID
        metadata, and writes the sequences to a FASTA file. FASTA headers use
        1-based coordinates.

        Output Format (per sequence):

            >{gene_id}::source={source};species={species};annotation_release={annotation_release};
            genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id}
            ::{chromosome}:{start}-{end}({strand})
            sequence

        :return: Path to the FASTA file with gene sequences.
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

        # BED12 fields: start must be 0-based for BED sequence extraction.
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
        file_fasta = safe_append_filename(self.dir_output, f"{self.FILE_INFO}__annotation_type__gene.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_intergenic(self) -> str:
        """
        Build FASTA sequences for intergenic regions from gene annotations.

        For each chromosome, gaps between genes are taken separately on the plus
        and minus strands. Chromosome lengths come from the genome FASTA. The
        resulting intervals are written to a FASTA file. FASTA headers use
        1-based coordinates.

        Output Format (per sequence):

            >{intergenic_region_id}::source={source};species={species};annotation_release={annotation_release};
            genome_assembly={genome_assembly};regiontype={regiontype}::{chromosome}:{start}-{end}({strand})
            sequence

        :return: Path to the FASTA file with intergenic sequences.
        :rtype: str
        """

        def _compute_intergenic_annotation(
            annotation: pd.DataFrame, file_chromosome_length: str
        ) -> pd.DataFrame:
            """
            Find intergenic intervals for every chromosome on both strands.

            Uses gene annotations and chromosome lengths to define the regions
            between genes (or a whole chromosome when that strand has no genes).

            :param annotation: Gene annotations for the genome.
            :type annotation: pd.DataFrame
            :param file_chromosome_length: Path to the chrom.sizes file with sequence lengths.
            :type file_chromosome_length: str
            :return: Intergenic intervals with coordinates, region IDs, and strand.
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
                        gene_annotation=gene_annotation_plusstrand,
                        strand="+",
                        file_chromosome_length=file_chromosome_length,
                    )
                )
                intergenic_annotation.append(
                    _compute_intergenic_annotation_strand(
                        seqid=seqid,
                        gene_annotation=gene_annotation_minusstrand,
                        strand="-",
                        file_chromosome_length=file_chromosome_length,
                    )
                )

            intergenic_annotation_df: pd.DataFrame = pd.concat(intergenic_annotation, ignore_index=True)
            return intergenic_annotation_df

        def _compute_intergenic_annotation_strand(
            seqid: str, gene_annotation: pd.DataFrame, strand: str, file_chromosome_length: str
        ) -> pd.DataFrame:
            """
            Find intergenic intervals on one chromosome and strand.

            If that strand has no genes, the interval spans the full chromosome.
            Otherwise, gene intervals are written as BED, and the gaps between
            them are taken with bedtools complement via :class:`~oligo_designer_toolsuite.utils.BedParser`.

            :param seqid: Chromosome or sequence name.
            :type seqid: str
            :param gene_annotation: Gene intervals on this chromosome and strand.
            :type gene_annotation: pd.DataFrame
            :param strand: Strand to process (``+`` or ``-``).
            :type strand: str
            :param file_chromosome_length: Path to the chrom.sizes file with sequence lengths.
            :type file_chromosome_length: str
            :return: Intergenic intervals for this chromosome and strand.
            :rtype: pd.DataFrame
            :raises ConfigurationError: If ``strand`` is not ``+`` or ``-``.
            """

            # case 1: no annotated genes on the respective chromosome and strand
            if gene_annotation.empty:
                chromosome_length = pd.read_csv(
                    file_chromosome_length, sep="\t", header=None, names=["seqid", "length"]
                )
                if strand == "+":
                    region_id_name = "InterRegPlus"
                elif strand == "-":
                    region_id_name = "InterRegMinus"
                else:
                    raise ConfigurationError(f"Invalid strand value: '{strand}'. Expected '+' or '-'.")
                # Whole chromosome is intergenic: BED start 0, header start 1.
                intergenic_annotation = pd.DataFrame(
                    {
                        "seqid": seqid,
                        "start_0base": 0,
                        "end": chromosome_length.length[chromosome_length.seqid == seqid],
                        "start_1base": 1,
                        "region_id": region_id_name + str(seqid) + "_0",
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
                gene_annotation = gene_annotation.sort_values(by="start")
                self.bed_parser.write_bed(gene_annotation, file_bed_in)

                # get complementary regions
                get_complement_regions(file_bed_in, file_chromosome_length, file_bed_out)

                # load intergenic regions
                intergenic_annotation = self.bed_parser.read_bed(
                    file_bed_out, names=["seqid", "start_0base", "end"]
                )
                # bedtools complement returns BED (0-based); convert for FASTA headers.
                intergenic_annotation["start_1base"] = self.bed_parser.convert_start(
                    intergenic_annotation["start_0base"], "0-based", "1-based"
                )
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
        file_chromosome_length = self._get_chromosome_length()

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

        # BED12 fields: start must be 0-based for BED sequence extraction.
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
        file_fasta = safe_append_filename(
            self.dir_output, f"{self.FILE_INFO}__annotation_type__intergenic.fna"
        )
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        os.remove(file_chromosome_length)
        del annotation

        return file_fasta

    def get_sequence_exon(self, collapse_duplicated_regions: bool = True) -> str:
        """
        Build FASTA sequences for exon intervals from exon annotations.

        Optionally merges exons that share the same start and end but come from
        different transcripts. FASTA headers include transcript IDs, exon numbers,
        and the total transcript count per gene. Headers use 1-based coordinates.

        Output Format (per sequence):

            >{gene_id}::source={source};species={species};annotation_release={annotation_release};
            genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
            exon_number={exon_number_x};transcript_id={transcript_id_b},exon_number={exon_number_y};
            number_total_transcripts={number_total_transcripts}::{chromosome}:{start}-{end}({strand})
            sequence

        :param collapse_duplicated_regions: If ``True``, merge exons with identical coordinates into one entry. Defaults to ``True``.
        :type collapse_duplicated_regions: bool
        :return: Path to the FASTA file with exon sequences.
        :rtype: str
        """

        # get exon annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "exon")

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

        # add transcript counts for each gene
        annotation, annotation_transcript_inf = self._add_transcript_counts(annotation)

        # BED12 fields: start must be 0-based for BED sequence extraction.
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
            + annotation_transcript_inf
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = safe_append_filename(self.dir_output, f"{self.FILE_INFO}__annotation_type__exon.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_intron(self, collapse_duplicated_regions: bool = True) -> str:
        """
        Build FASTA sequences for intron intervals derived from exon annotations.

        Introns are the gaps between consecutive exons within each transcript.
        Optionally merges introns that share the same start and end but come from
        different transcripts. FASTA headers include transcript IDs and intron
        numbers. Headers use 1-based coordinates.

        Output Format (per sequence):

            >{gene_id}::source={source};species={species};annotation_release={annotation_release};
            genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
            intron_number={intron_number_x};transcript_id={transcript_id_b},intron_number={intron_number_y};
            number_total_transcripts={number_total_transcripts}::{chromosome}:{start}-{end}({strand})
            sequence

        :param collapse_duplicated_regions: If ``True``, merge introns with identical coordinates into one entry. Defaults to ``True``.
        :type collapse_duplicated_regions: bool
        :return: Path to the FASTA file with intron sequences.
        :rtype: str
        """

        def _compute_intron_annotation(annotation: pd.DataFrame) -> pd.DataFrame:
            """
            Derive intron intervals from exon annotations for each transcript.

            Finds the gaps between consecutive exons within a transcript and
            numbers them by position, accounting for strand orientation. Returns
            intron rows with gene ID, transcript ID, and intron number.

            :param annotation: Exon annotations for the genome.
            :type annotation: pd.DataFrame
            :return: Intron intervals with gene, transcript, and intron number.
            :rtype: pd.DataFrame
            :raises ConfigurationError: If no introns can be derived from the exons.
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

            if len(intron_list) == 0:
                raise ConfigurationError("Could not calculate introns.")

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

        # BED12 fields: start must be 0-based for BED sequence extraction.
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
        file_fasta = safe_append_filename(self.dir_output, f"{self.FILE_INFO}__annotation_type__intron.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_CDS(self, collapse_duplicated_regions: bool = True) -> str:
        """
        Build FASTA sequences for coding sequence (CDS) intervals from CDS annotations.

        Optionally merges CDS intervals that share the same start and end but come
        from different transcripts. FASTA headers include transcript IDs, exon
        numbers, and the total transcript count per gene. Headers use 1-based
        coordinates.

        Output Format (per sequence):

            >{gene_id}::source={source};species={species};annotation_release={annotation_release};
            genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
            exon_number={exon_number_x};transcript_id={transcript_id_b},exon_number={exon_number_y};
            number_total_transcripts={number_total_transcripts}::{chromosome}:{start}-{end}({strand})
            sequence

        :param collapse_duplicated_regions: If ``True``, merge CDS intervals with identical coordinates into one entry. Defaults to ``True``.
        :type collapse_duplicated_regions: bool
        :return: Path to the FASTA file with CDS sequences.
        :rtype: str
        """
        # get exon annotation entries
        annotation = self._load_annotation()
        annotation = self._get_annotation_region_of_interest(annotation, "CDS")

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

        # add transcript counts for each gene
        annotation, annotation_transcript_inf = self._add_transcript_counts(annotation)

        # BED12 fields: start must be 0-based for BED sequence extraction.
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
            + annotation_transcript_inf
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = safe_append_filename(self.dir_output, f"{self.FILE_INFO}__annotation_type__cds.fna")
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
        Build FASTA sequences for UTR intervals from exon and CDS annotations.

        Derives 5' and/or 3' untranslated regions from exon spans outside the CDS
        boundaries of each transcript. Use ``five_prime`` and ``three_prime`` to
        choose which ends to keep. Optionally merges UTRs with identical coordinates.
        FASTA headers use 1-based coordinates.

        Output Format (per sequence):

            >{gene_id}::source={source};species={species};annotation_release={annotation_release};
            genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
            exon_number={exon_number_x};transcript_id={transcript_id_b},exon_number={exon_number_y};
            number_total_transcripts={number_total_transcripts}::{chromosome}:{start}-{end}({strand})
            sequence

        :param five_prime: If ``True``, include 5' UTR sequences. Defaults to ``True``.
        :type five_prime: bool
        :param three_prime: If ``True``, include 3' UTR sequences. Defaults to ``True``.
        :type three_prime: bool
        :param collapse_duplicated_regions: If ``True``, merge UTRs with identical coordinates into one entry. Defaults to ``True``.
        :type collapse_duplicated_regions: bool
        :return: Path to the FASTA file with UTR sequences.
        :rtype: str
        :raises ConfigurationError: If no UTR intervals can be derived.
        """

        def _compute_UTR(annotation: pd.DataFrame) -> pd.DataFrame:
            """
            Derive 5' and 3' UTR intervals from exon and CDS annotations.

            Uses CDS boundaries within each transcript to trim exons outside the
            coding span. Strand orientation decides which side is 5' versus 3'.

            :param annotation: Combined exon and CDS annotations for transcripts with CDS.
            :type annotation: pd.DataFrame
            :return: UTR intervals labeled as five_prime_UTR or three_prime_UTR.
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
                # Keep 1-based and 0-based starts aligned when trimming to the CDS edge.
                UTR_right.loc[UTR_right["start_1base"] <= cds_end, "start_1base"] = cds_end + 1
                UTR_right.loc[(UTR_right["start_0base"] + 1) <= cds_end, "start_0base"] = cds_end
                utrs.append(UTR_right)

            utr_annotation: pd.DataFrame = pd.concat(utrs, ignore_index=True)

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

        if annotation.shape[0] == 0:
            raise ConfigurationError("Could not calculate the UTR.")

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

        # add transcript counts for each gene
        annotation, annotation_transcript_inf = self._add_transcript_counts(annotation)

        # BED12 fields: start must be 0-based for BED sequence extraction.
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
            + annotation_transcript_inf
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        utr_suffix = "__".join(
            ["utr"] + (["five_prime"] if five_prime else []) + (["three_prime"] if three_prime else [])
        )

        file_fasta = safe_append_filename(
            self.dir_output, f"{self.FILE_INFO}__annotation_type__{utr_suffix}.fna"
        )
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_exon_exon_junction(
        self, block_size: int, collapse_duplicated_regions: bool = True
    ) -> str:
        """
        Build FASTA sequences spanning exon-exon junctions from exon annotations.

        For consecutive exons in a transcript, takes ``block_size`` bases on each
        side of the junction (or the full exon if shorter). Optionally merges
        junctions with identical coordinates from different transcripts. FASTA
        headers include transcript IDs and joined exon numbers. Headers use
        1-based coordinates.

        Output Format (per sequence):

            >{gene_id}::source={source};species={species};annotation_release={annotation_release};
            genome_assembly={genome_assembly};regiontype={regiontype};gene_id={gene_id};transcript_id={transcript_id_a},
            exon_number={exon_number_x__JUNC__exon_number_y};transcript_id={transcript_id_b},
            exon_number={exon_number_y__JUNC__exon_number_z};number_total_transcripts={number_total_transcripts}
            ::{chromosome}:{start}-{end}({strand})
            sequence

        :param block_size: Number of bases to take on each side of the junction.
        :type block_size: int
        :param collapse_duplicated_regions: If ``True``, merge junctions with identical coordinates into one entry. Defaults to ``True``.
        :type collapse_duplicated_regions: bool
        :return: Path to the FASTA file with exon-exon junction sequences.
        :rtype: str
        :raises ConfigurationError: If no exon-exon junctions can be derived.
        """

        def _compute_exon_exon_junction_annotation(annotation: pd.DataFrame, block_size: int) -> pd.DataFrame:
            """
            Derive exon-exon junction intervals from exon annotations.

            For consecutive exons in a transcript, builds a junction spanning
            ``block_size`` bases on each side. If an exon is shorter than
            ``block_size``, neighboring exons are included so the junction sequence
            stays long enough.

            :param annotation: Exon annotations for the genome.
            :type annotation: pd.DataFrame
            :param block_size: Number of bases to take on each side of the junction.
            :type block_size: int
            :return: Junction intervals with block sizes for BED12-style extraction.
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
                        # Header coords are 1-based; start_up/end_down below stay 0-based for BED.
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

        if annotation.shape[0] == 0:
            raise ConfigurationError("Could not calculate exon-exon junctions.")

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

        # add transcript counts for each gene
        annotation, annotation_transcript_inf = self._add_transcript_counts(annotation)

        # BED12 fields: junction "start" is already 0-based from the exon walk above.
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
            + annotation_transcript_inf
            + SEPARATOR_FASTA_HEADER_FIELDS
            + annotation["region"]
        )
        annotation = annotation[self.BED12_HEADER]

        file_fasta = safe_append_filename(
            self.dir_output,
            f"{self.FILE_INFO}__annotation_type__exon_exon_junction__block_size__{block_size}.fna",
        )
        self._get_sequence_from_annotation(annotation, file_fasta)

        del annotation

        return file_fasta

    def _get_chromosome_length(self) -> str:
        """
        Write chromosome lengths from the genome FASTA to a chrom.sizes file.

        The file is used by bedtools when finding regions outside annotated
        genes. Each line is ``seqid`` and length, tab-separated, with no header.

        :return: Path to the chromosome-length file.
        :rtype: str
        """
        chromosome_lengths = {}
        for rec in SeqIO.parse(self.sequence_file, "fasta"):
            chromosome_lengths[rec.id] = len(rec.seq)

        file_chromosome_length = os.path.join(self.dir_output, "annotation.genome")
        with open(file_chromosome_length, "w") as handle:
            for key, value in sorted(chromosome_lengths.items()):
                handle.write(f"{key}\t{value}\n")

        return file_chromosome_length

    def _load_annotation(self) -> pd.DataFrame:
        """
        Load the parsed GFF annotation and add 0-based and 1-based start columns.

        Reads the pickled annotation written at construction time. GFF starts are
        1-based and kept as ``start_1base`` for FASTA headers. A matching
        ``start_0base`` column is added for BED-style sequence extraction.

        :return: Annotation table with ``start_1base`` and ``start_0base`` columns.
        :rtype: pd.DataFrame
        """
        annotation: pd.DataFrame = self.gff_parser.load_annotation_from_pickle(self.parsed_annotation_file)

        # Required so numeric sorting of genomic coordinates is correct.
        annotation.start = annotation.start.astype("int")
        annotation.end = annotation.end.astype("int")

        # GFF starts are 1-based; keep both so headers and BED extraction stay aligned.
        annotation.rename(columns={"start": "start_1base"}, inplace=True)
        annotation["start_0base"] = annotation.start_1base - 1

        return annotation

    def _get_annotation_region_of_interest(self, annotation: pd.DataFrame, region: str) -> pd.DataFrame:
        """
        Keep only annotation rows of a given feature type.

        Filters on the ``type`` column (for example ``gene``, ``exon``, or ``CDS``).

        :param annotation: Full annotation table.
        :type annotation: pd.DataFrame
        :param region: Feature type to keep.
        :type region: str
        :return: Rows whose ``type`` matches ``region``.
        :rtype: pd.DataFrame
        :raises ConfigurationError: If ``region`` is not present in the annotation.
        """
        region_annotation = annotation.loc[annotation["type"] == region]
        if region_annotation.shape[0] == 0:
            raise ConfigurationError(f"type {region} is not contained in the annotation file.")
        region_annotation.reset_index(inplace=True, drop=True)
        return region_annotation

    def _get_attributes(
        self,
        start_0base: int | None = None,
        start_1base: int | None = None,
        end: int | None = None,
        gene_id: str | None = None,
        exon_number: int | None = None,
        exon_size: int | None = None,
    ) -> pd.Series:
        """
        Bundle optional genomic feature fields into a single series.

        Used when walking exons within a transcript to carry coordinates and
        exon metadata between steps.

        :param start_0base: Feature start in 0-based coordinates.
        :type start_0base: int, optional
        :param start_1base: Feature start in 1-based coordinates.
        :type start_1base: int, optional
        :param end: Feature end coordinate.
        :type end: int, optional
        :param gene_id: Gene identifier for the feature.
        :type gene_id: str, optional
        :param exon_number: Exon number within the transcript.
        :type exon_number: int, optional
        :param exon_size: Length of the exon in bases.
        :type exon_size: int, optional
        :return: Series with the provided attribute fields.
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
        Merge annotation rows that share the same genomic region.

        Groups by the ``region`` string and keeps one row per unique interval.
        Extra info strings from duplicate rows are joined so transcript and exon
        labels from all isoforms stay in the FASTA header.

        :param annotation: Annotation rows that may contain duplicate intervals.
        :type annotation: pd.DataFrame
        :return: Annotation with one row per unique ``region``.
        :rtype: pd.DataFrame
        """
        aggregate_function: dict[str, str | Callable[[Iterable[str]], str]] = {
            col: "first" for col in annotation.columns
        }
        aggregate_function["add_inf"] = SEPARATOR_FASTA_HEADER_FIELDS_LIST.join

        merged_annotation = annotation.groupby(annotation["region"]).agg(aggregate_function)
        merged_annotation.reset_index(inplace=True, drop=True)

        return merged_annotation

    def _get_annotation_region(self, annotation: pd.DataFrame) -> pd.Series:
        """
        Format each annotation row as a 1-based coordinate string for FASTA headers.

        Builds ``seqid:start-end(strand)`` from ``seqid``, ``start_1base``, ``end``,
        and ``strand``. Uses 1-based starts so the header matches GFF/GTF coordinates.

        :param annotation: Annotation rows with coordinate and strand columns.
        :type annotation: pd.DataFrame
        :return: Series of region strings, one per annotation row.
        :rtype: pd.Series
        """
        region: pd.Series = (
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
        Write sequences for annotated intervals to a FASTA file.

        Intervals are written as a temporary BED file with
        :class:`~oligo_designer_toolsuite.utils.BedParser`, then sequences are
        pulled from the genome FASTA. The BED ``start`` column must already be
        0-based. The temporary BED file is removed afterward.

        :param annotation: Intervals to extract. Must include the columns needed for BED and FASTA headers.
        :type annotation: pd.DataFrame
        :param file_fasta: Path of the FASTA file to create.
        :type file_fasta: str
        :param split: If ``True`` and the BED is BED12-style, concatenate block sequences (for example exons). Defaults to ``True``.
        :type split: bool
        :param strand: If ``True``, reverse-complement intervals on the minus strand. Defaults to ``True``.
        :type strand: bool
        :return: None
        :rtype: None
        """
        annotation = annotation.sort_values(by=["fasta_header"])
        annotation.reset_index(inplace=True, drop=True)

        # BED extraction expects 0-based starts (usually annotation["start"] = start_0base upstream).
        id = random.randint(0, 10000000)
        file_bed = safe_append_filename(self.dir_output, f"annotation_{id}.bed")
        self.bed_parser.write_bed(annotation, file_bed)

        get_sequence_from_annotation(
            file_bed,
            self.sequence_file,
            file_fasta,
            split=split,
            strand=strand,
            nameOnly=True,
        )
        os.remove(file_bed)

    def _get_number_total_transcripts(self, gene_ids: set[str]) -> pd.DataFrame | None:
        """
        Count transcripts per gene for the requested gene IDs.

        Uses ``transcript`` features in the annotation when available. Returns
        ``None`` if transcript counts cannot be computed for all requested genes.

        :param gene_ids: Gene IDs that need a transcript count. Some assemblies lack transcript rows for every gene.
        :type gene_ids: set[str]
        :return: Table of ``gene_id`` and ``transcript_count``, or ``None`` if counts are unavailable.
        :rtype: pd.DataFrame | None
        """
        annotation = self._load_annotation()
        try:
            annotation_interest = self._get_annotation_region_of_interest(annotation, "transcript")
            number_total_transcripts = annotation_interest["gene_id"].value_counts()
            number_total_transcripts_df = number_total_transcripts.reset_index()
            number_total_transcripts_df.columns = ["gene_id", "transcript_count"]
            if not set(gene_ids).issubset(set(number_total_transcripts_df["gene_id"])):
                raise ConfigurationError
        except ConfigurationError:
            logger.warning("Could not calculate the number of total transcripts.")
            number_total_transcripts_df = None
        finally:
            return number_total_transcripts_df

    def _add_transcript_counts(self, annotation: pd.DataFrame) -> tuple[pd.DataFrame, str | pd.Series]:
        """
        Attach per-gene transcript counts for FASTA header metadata.

        Looks up total transcripts per gene and builds a header fragment of the
        form ``number_total_transcripts=...``. If counts are unavailable, the
        fragment is an empty string.

        :param annotation: Annotation rows that include a ``gene_id`` column.
        :type annotation: pd.DataFrame
        :return: Updated annotation and a string or series for the FASTA header.
        :rtype: tuple[pd.DataFrame, str | pd.Series]
        """
        # add transcript counts for each gene
        number_total_transcripts = self._get_number_total_transcripts(set(annotation["gene_id"]))
        if number_total_transcripts is not None:
            annotation = pd.merge(annotation, number_total_transcripts, on="gene_id", how="left")
            annotation_transcript_inf = (
                f"{SEPARATOR_FASTA_HEADER_FIELDS_LIST}number_total_transcripts="
                + annotation["transcript_count"].astype("str")
            )
        else:
            annotation_transcript_inf = ""

        return annotation, annotation_transcript_inf


class NcbiGenomicRegionGenerator(CustomGenomicRegionGenerator):
    """
    Build genomic region FASTA files from NCBI annotation and genome downloads.

    Downloads GTF annotation and genome FASTA via :class:`~oligo_designer_toolsuite.sequence_generator.FtpLoaderNCBI`,
    then reuses :class:`CustomGenomicRegionGenerator` extractors. Choose the genome
    either by taxon, species, and annotation release (``mode='species'``), or by
    RefSeq assembly accession and assembly name (``mode='assembly'``). Only one
    mode may be used; unused parameters should be ``None``.

    :param mode: How to select the genome: ``species`` or ``assembly``.
    :type mode: str
    :param taxon: NCBI taxon folder (for example ``vertebrate_mammalian``). Used when ``mode='species'``.
    :type taxon: str, optional
    :param species: Species name (for example ``Homo_sapiens``). Used when ``mode='species'``.
    :type species: str, optional
    :param annotation_release: Annotation release version. Defaults to ``current`` when ``mode='species'``.
    :type annotation_release: str, optional
    :param assembly_source: NCBI assembly source. Supported values are ``auto``, ``annotation_releases``,
        ``latest_assembly_versions``, and ``reference``. Defaults to ``auto``.
    :type assembly_source: str, optional
    :param refseq_assembly_accession: RefSeq assembly accession (for example ``GCF_000001405.38``). Used when ``mode='assembly'``.
    :type refseq_assembly_accession: str, optional
    :param assembly_name: Assembly name (for example ``GRCh38.p12``). Used when ``mode='assembly'``.
    :type assembly_name: str, optional
    :param dir_output: Directory for output files. Defaults to ``output``.
    :type dir_output: str, optional
    """

    def __init__(
        self,
        mode: str | None = None,
        taxon: str | None = None,
        species: str | None = None,
        annotation_release: str | None = None,
        assembly_source: str | None = None,
        refseq_assembly_accession: str | None = None,
        assembly_name: str | None = None,
        dir_output: str = "output",
    ) -> None:
        """Constructor for the NcbiGenomicRegionGenerator class."""
        files_source = "NCBI"

        if mode is None:
            raise ConfigurationError("For source='ncbi', parameter 'mode' must be provided.")

        if mode == "species":
            if taxon is None:
                raise ConfigurationError(f"No taxon defined.")

            if species is None:
                raise ConfigurationError(f"No species defined.")

            if annotation_release is None:
                annotation_release = "current"
                logger.warning(f"No annotation release defined. Using default release {annotation_release}!")

            if assembly_source is None:
                assembly_source = "auto"
                logger.warning(
                    f"No assembly source defined. Using default assembly source {assembly_source}!"
                )
        elif mode == "assembly":
            if assembly_source is None:
                assembly_source = "auto"
        else:
            raise ConfigurationError("For source='ncbi', mode must be either 'species' or 'assembly'.")

        self.dir_output = os.path.join(dir_output, "annotation")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ftp = FtpLoaderNCBI(
            self.dir_output,
            mode=mode,
            taxon=taxon,
            species=species,
            annotation_release=annotation_release,
            assembly_source=assembly_source if assembly_source is not None else "auto",
            refseq_assembly_accession=refseq_assembly_accession,
            assembly_name=assembly_name,
        )
        annotation_file, annotation_release, genome_assembly = ftp.download_files("gtf")
        sequence_file, _, _ = ftp.download_files("fasta")

        super().__init__(
            annotation_file,
            sequence_file,
            files_source,
            species if species is not None else "unknown",
            annotation_release,
            genome_assembly,
            dir_output,
        )


class EnsemblGenomicRegionGenerator(CustomGenomicRegionGenerator):
    """
    Build genomic region FASTA files from Ensembl annotation and genome downloads.

    Downloads GTF annotation and genome FASTA via
    :class:`~oligo_designer_toolsuite.sequence_generator.FtpLoaderEnsembl`, then
    reuses :class:`CustomGenomicRegionGenerator` extractors for the chosen species
    and annotation release.

    :param species: Species name (for example ``homo_sapiens``). Defaults to ``homo_sapiens``.
    :type species: str, optional
    :param annotation_release: Annotation release version. Defaults to ``current``.
    :type annotation_release: str, optional
    :param dir_output: Directory for output files. Defaults to ``output``.
    :type dir_output: str, optional
    """

    def __init__(
        self,
        species: str | None = None,
        annotation_release: str | None = None,
        dir_output: str = "output",
    ) -> None:
        """Constructor for the EnsemblGenomicRegionGenerator class."""
        files_source = "Ensemble"
        if species is None:
            species = "homo_sapiens"
            logger.warning(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "current"
            logger.warning(f"No annotation release defined. Using default release {annotation_release}!")

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
