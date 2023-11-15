############################################
# imports
############################################

import os
import copy
import shutil
import warnings
import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None

from pathlib import Path
from Bio import SeqIO

from ..utils._data_parser import get_sequence_from_annotation, get_complement_regions
from ..utils._ftp_loader import FtpLoaderEnsembl, FtpLoaderNCBI
from ..utils._gff_parser import GffParser


############################################
# Genomic Region Generator Classes
############################################


class CustomGenomicRegionGenerator:
    """Class to generate sequences for different types of regions from annotation in GTF format
    and genomic fasta file.

    Sequences are saved as fasta file with region id, additional information and coordinates in header.
    The header of each sequence must start with '>' and contain the following information:
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
    where the region_id is compulsory and the other fileds are opional. Coordinated are saved in 1-base format.

    Output Format (per sequence):
    >{region_id}::{additional information}::{chromosome}:{start}-{end}({strand})
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    :param annotation_file: Annotation of genomic regions. Has to be provided in GTF format.
    :type annotation_file: str
    :param sequence_file: Fasta file with genome sequence.
    :type sequence_file: str
    :param files_source: Source of annotations, e.g. NCBI, defaults to None.
    :type files_source: str, optional
    :param species: Species of annotation, e.g. Homo_sapiens, defaults to None.
    :type species: str, optional
    :param annotation_release: Release number of annotation, e.g. 110, defaults to None.
    :type annotation_release: str, optional
    :param genome_assembly: Genome assembly of annotation, e.g. GRCh38, defaults to None.
    :type genome_assembly: str, optional
    :param dir_output: Output directory, defaults to 'output'.
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
        """Constructor"""
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

        self.files_source = files_source
        self.species = species
        self.annotation_release = annotation_release
        self.genome_assembly = genome_assembly
        self.annotation_file = annotation_file
        self.sequence_file = sequence_file

        self.dir_output = os.path.join(dir_output, "annotation")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # read annotation file and store in dataframe
        parser = GffParser()
        self.annotation = parser.read_gff(annotation_file)

        # required to ensure that sorting is done correctly
        self.annotation.start = self.annotation.start.astype("int")
        self.annotation.end = self.annotation.end.astype("int")

        # add both annotations to dataframe: GFF 1-base offset and BED 0-base offset
        # since we read in a GFF file, the start coordinates are 1-base offset
        self.annotation.rename(columns={"start": "start_1base"}, inplace=True)
        self.annotation["start_0base"] = self.annotation.start_1base - 1

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
        """Generate a fasta file with all gene sequences, using the attribute "gene" in the GTF annotation type field.
        The header of each sequence in the fasta file contains the following information:

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id}::{chromosome}:{start}-{end}({strand})
        sequence

        :return: Fasta file with gene sequences.
        :rtype: str
        """
        # get gene annotation entries
        annotation = copy.deepcopy(self.annotation)
        annotation = self._get_annotation_region_of_interest(annotation, "gene")

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = "regiontype=gene;gene_id=" + annotation["gene_id"].astype("str")
        annotation["region"] = self._get_annotation_region(annotation)

        # add BED12 fields
        annotation["start"] = annotation["start_0base"]
        annotation["score"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"] + "::" + annotation["add_inf"] + "::" + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        # get sequence from bed file
        file_fasta = os.path.join(self.dir_output, f"gene_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_intergenic(self):
        """Generate a fasta file with all intergenic sequences. Using the attribute "gene" in the GTF annotation type field,
        we substract those regions from the region of the whole chromosome and retrieve the sequence for those regions
        on the plus and minus strand seperately, according to the gene regions strandness.
        The header of each sequence in the fasta file contains the following information:

        Output Format (per sequence):
        >{intergenic_region_id}::{chromosome}:{start}-{end}({strand})
        sequence

        :return: Fasta file with gene sequences.
        :rtype: str
        """

        def _compute_intergenic_annotation(annotation):
            """Retrieve the annotation of intergenic regions by retrieving the complement to the
            gene regions, taking the strandness of genes into account.

            :param annotation: Gene annotation in GTF format.
            :type annotation: pandas.DataFrame
            :return: Intergenic regions annotation.
            :rtype: pandas.DataFrame
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
            """Retrieve strandspecific annotation of intergenic regions.

            :param gene_annotatio: Gene annotation in GTF format.
            :type gene_annotatio: pandas.DataFrame
            :param strand: Genomic strand of genes in annotation.
            :type strand: str
            :return: Intergenic Region annotation for one strand.
            :rtype: list
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
            """Get the length, i.e. last coordinate, of each chromosome in the input fasta file.

            :return: Chromosome length.
            :rtype: dict
            """
            dict_chromosome_length = {}
            for rec in SeqIO.parse(self.sequence_file, "fasta"):
                dict_chromosome_length[rec.id] = len(rec.seq)

            with open(file_chromosome_length, "w") as handle:
                for key, value in sorted(dict_chromosome_length.items()):
                    handle.write(f"{key}\t{value}\n")

            return dict_chromosome_length

        # get gene annotation entries
        annotation = copy.deepcopy(self.annotation)
        annotation = self._get_annotation_region_of_interest(annotation, "gene")
        annotation = _compute_intergenic_annotation(annotation)

        # generate region_id
        annotation["add_inf"] = "regiontype=intergenic"
        annotation["region"] = self._get_annotation_region(annotation)

        # add BED12 fields
        annotation["start"] = annotation["start_0base"]
        annotation["score"] = 0
        annotation["fasta_header"] = (
            annotation["region_id"] + "::" + annotation["add_inf"] + "::" + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        # get sequence from bed file
        file_fasta = os.path.join(self.dir_output, f"intergenic_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_exon(self, collapse_duplicated_regions=True):
        """Generate a fasta file with all exon sequences, using the attribute "exon" in the GTF annotation type field.
        If `collapse_duplicated_regions` is set to True (default), then exons from different transcripts who share
        the exact same start and end coordinates, are merged into one sequence entry.
        The header of each sequence in the fasta file contains the following information:

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},exon_number={exon_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param collapse_duplicated_regions: merge region entries with the exact same start and end coordinates (strandspecific), defaults to True
        :type collapse_duplicated_regions: bool, optional
        :return: Fasta file with exon sequences.
        :rtype: str
        """
        # get exon annotation entries
        annotation = copy.deepcopy(self.annotation)
        annotation = self._get_annotation_region_of_interest(annotation, "exon")

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "gene_id="
            + annotation["gene_id"].astype("str")
            + ",transcript_id="
            + annotation["transcript_id"].astype("str")
            + ",exon_number="
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
            + "::regiontype=exon;"
            + annotation["add_inf"]
            + "::"
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = os.path.join(self.dir_output, f"exon_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_intron(self, collapse_duplicated_regions=True):
        """Generate a fasta file with all intron sequences, using the attribute "exon" in the GTF annotation type field.
        Intron are specified per transcript and calculated from the regions between exons. If `collapse_duplicated_regions`
        is set to True (default), then introns from different transcripts who share the exact same start and end coordinates,
        are merged into one sequence entry. The header of each sequence in the fasta file contains the following information:

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},intron_number={intron_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :return: Fasta file with exon sequences.
        :rtype: str

        :param collapse_duplicated_regions: merge region entries with the exact same start and end coordinates (strandspecific), defaults to True
        :type collapse_duplicated_regions: bool, optional
        :return: Fasta file with intron sequences.
        :rtype: str
        """

        def _compute_intron_annotation(annotation):
            """Compute the start and end coordinates of introns from the exon annotation.

            :param annotation: Exon annotation in GTF format.
            :type annotation: pandas.DataFrame
            :return: Intron annotation in GTF format.
            :rtype: pandas.DataFrame
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

        annotation = copy.deepcopy(self.annotation)
        annotation = self._get_annotation_region_of_interest(annotation, "exon")
        annotation = _compute_intron_annotation(annotation)

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "gene_id="
            + annotation["gene_id"].astype("str")
            + ",transcript_id="
            + annotation["transcript_id"].astype("str")
            + ",intron_number="
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
            + "::regiontype=intron;"
            + annotation["add_inf"]
            + "::"
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        # get sequence from bed file
        file_fasta = os.path.join(self.dir_output, f"intron_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_CDS(self, collapse_duplicated_regions=True):
        """Generate a fasta file with all CDS sequences, using the attribute "CDS" in the GTF annotation type field.
        If `collapse_duplicated_regions` is set to True (default), then CDS from different transcripts who share
        the exact same start and end coordinates, are merged into one sequence entry.
        The header of each sequence in the fasta file contains the following information:

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},exon_number={exon_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param collapse_duplicated_regions: merge region entries with the exact same start and end coordinates (strandspecific), defaults to True
        :type collapse_duplicated_regions: bool, optional
        :return: Fasta file with CDS sequences.
        :rtype: str
        """
        # get exon annotation entries
        annotation = copy.deepcopy(self.annotation)
        annotation = self._get_annotation_region_of_interest(annotation, "CDS")

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "gene_id="
            + annotation["gene_id"].astype("str")
            + ",transcript_id="
            + annotation["transcript_id"].astype("str")
            + ",exon_number="
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
            + "::regiontype=CDS;"
            + annotation["add_inf"]
            + "::"
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = os.path.join(self.dir_output, f"CDS_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_UTR(self, five_prime=True, three_prime=True, collapse_duplicated_regions=True):
        """Generate a fasta file with all UTR sequences, using the attribute "CDS" and "exon" in the GTF
        annotation type field. The UTR represents regions from the transcription start/end site or beginning
        of the known UTR to the base before the start codon of the transcript. If this region is interrupted
        by introns then each exon or partial exon is annotated as a separate UTR feature. 3 prime and 5 prime UTRS
        are defined dependent on the strand of the gene. If `collapse_duplicated_regions` is set to True (default),
        then UTRs from different transcripts who share the exact same start and end coordinates,are merged into one
        sequence entry. The header of each sequence in the fasta file contains the following information:

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},exon_number={exon_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param five_prime: include 5 prime UTR regions, defaults to True
        :type five_prime: bool, optional
        :param three_prime: include 3 prime UTR regions, defaults to True
        :type three_prime: bool, optional
        :param collapse_duplicated_regions: merge region entries with the exact same start and end coordinates (strandspecific), defaults to True
        :type collapse_duplicated_regions: bool, optional
        :return: Fasta file with UTR sequences.
        :rtype: str
        """

        def _compute_UTR(annotation):
            """Compute the UTR boundaries, and define 3 prime and 5 prime UTR dependent on the strandness of the gene.

            :param annotation: Region annotation.
            :type annotation: pandas.DataFrame
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

                UTR_left = copy.deepcopy(exons[exons.start_1base < cds_start])
                UTR_left.type = UTR_left_type
                UTR_left.end[UTR_left.end >= cds_start] = cds_start - 1
                utrs.append(UTR_left)

                UTR_right = copy.deepcopy(exons[exons.end > cds_end])
                UTR_right.type = UTR_right_type
                UTR_right.start_1base[UTR_right.start_1base <= cds_end] = cds_end + 1
                UTR_right.start_0base[UTR_right.start_1base <= cds_end] = cds_end
                utrs.append(UTR_right)

            utr_annotation = pd.concat(utrs, ignore_index=True)

            return utr_annotation

        annotation_exon = self._get_annotation_region_of_interest(copy.deepcopy(self.annotation), "exon")
        annotation_CDS = self._get_annotation_region_of_interest(copy.deepcopy(self.annotation), "CDS")

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
            + ",transcript_id="
            + annotation["transcript_id"].astype("str")
            + ",exon_number="
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
            + "::regiontype="
            + annotation["type"]
            + ";"
            + annotation["add_inf"]
            + "::"
            + annotation["region"]
        )
        annotation = annotation[self.BED_HEADER]

        file_fasta = os.path.join(self.dir_output, f"UTR_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta, split=False)

        del annotation

        return file_fasta

    def get_sequence_exon_exon_junction(self, block_size, collapse_duplicated_regions=True):
        """Generate a fasta file with all exon/exon junction sequences, using the attribute "exon" in the GTF
        annotation type field. To compute the exon junctions, the "block_size" parameter defines the size
        of the exon junction region, i.e. +/- "block_size" bp around the junction. If `collapse_duplicated_regions`
        is set to True (default), then junctions from different transcripts who share the exact same start and end coordinates,
        are merged into one sequence entry. The header of each sequence in the fasta file contains the following information:

        Output Format (per sequence):
        >{gene_id}::gene_id={gene_id},transcript_id={transcript_id},exon_number={exon_number}__JUNC__{exon_number}::{chromosome}:{start}-{end}({strand})
        sequence

        :param block_size: Size of region, i.e. +/- "block_size" bp around the junction
        :type block_size: int
        :param collapse_duplicated_regions: merge region entries with the exact same start and end coordinates (strandspecific), defaults to True
        :type collapse_duplicated_regions: bool, optional
        :return: Fasta file with exon/exon junction sequences.
        :rtype: str
        """

        def _compute_exon_exon_junction_annotation(annotation, block_size):
            """Compute the exon/exon junction region from the exon start and end coordinates, depending
            on the strandness of the gene.

            :param annotation: Exon annotation.
            :type annotation: pandas.DataFrame
            :param block_size: Size of region, i.e. +/- "block_size" bp around the junction
            :type block_size: int
            :return: Exon/exon junction annotation.
            :rtype: pandas.DataFrame
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
                            regions_exons_small = ";".join(
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
                                ";".join(
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
        annotation = copy.deepcopy(self.annotation)
        annotation = self._get_annotation_region_of_interest(annotation, "exon")

        # compute exon junctions
        annotation = _compute_exon_exon_junction_annotation(annotation, block_size)

        # generate region_id
        annotation["region_id"] = annotation["gene_id"].astype("str")
        annotation["add_inf"] = (
            "gene_id="
            + annotation["gene_id"].astype("str")
            + ",transcript_id="
            + annotation["transcript_id"].astype("str")
            + ",exon_number="
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
            + "::regiontype=exonexonjunction;"
            + annotation["add_inf"]
            + "::"
            + annotation["region"]
        )
        annotation = annotation[self.BED12_HEADER]

        file_fasta = os.path.join(self.dir_output, f"exon_exon_junction_annotation_{self.FILE_INFO}.fna")
        self._get_sequence_from_annotation(annotation, file_fasta)

        del annotation

        return file_fasta

    def _get_annotation_region_of_interest(self, annotation, region):
        """Retrieve annotation ofr region of interest from loaded annotation dataframe.
        Annotation is filtered based on the 3rd "type" column in the GFT file.

        :param region: Region to extract.
        :type region: str {'gene', 'transcript', 'exon', 'CDS'}
        :return: Region annotation.
        :rtype: pandas.DataFrame
        """
        region_annotation = annotation.loc[annotation["type"] == region]
        region_annotation.reset_index(inplace=True, drop=True)
        return region_annotation

    def _get_attributes(
        self,
        start_0base=None,
        start_1base=None,
        end=None,
        gene_id=None,
        exon_number=None,
        exon_size=None,
    ):
        """Helper function to define attributes for a specific region.

        :param start_0base: Value for start_0base, defaults to None
        :type start_0base: int, optional
        :param start_1base: Value for start_1base, defaults to None
        :type start_1base: int, optional
        :param end: Value for end, defaults to None
        :type end: int, optional
        :param gene_id: Value for gene_id, defaults to None
        :type gene_id: str, optional
        :param exon_number: Value for exon_number, defaults to None
        :type exon_number: int, optional
        :param exon_size: Value for exon_size, defaults to None
        :type exon_size: int, optional
        :return: _description_
        :rtype: _type_
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

    def _collapse_duplicated_regions(self, annotation):
        """Merge overlapping regions, which have the exact same start and end coordinates.
        Save all additional information for merged regions.

        :param annotation: Region annotation.
        :type annotation: pandas.DataFrame
        :return: Merged region annotation.
        :rtype: pandas.DataFrame
        """
        aggregate_function = {col: "first" for col in annotation.columns}
        aggregate_function["add_inf"] = ";".join

        merged_annotation = annotation.groupby(annotation["region"]).aggregate(aggregate_function)
        merged_annotation.reset_index(inplace=True, drop=True)

        return merged_annotation

    def _get_annotation_region(self, annotation):
        """Get coordinates of region in string format as: 'chromosome':'start'-'end'('strand')

        :param annotation: Region annotation.
        :type annotation: pandas.DataFrame
        :return: Region coordinates in string format.
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

    def _get_sequence_from_annotation(self, annotation, file_fasta, split=True, strand=True):
        """Retrieve fasta sequence from bed file annotation.

        :param annotation: Region annotation in bed or bed12 format.
        :type annotation: pandas.DataFrame
        :param file_fasta: Name of output fasta file.
        :type file_fasta: str
        :param split: Enable split read as defined in bed12 format, defaults to True
        :type split: bool, optional
        :param strand: Get sequence strand specific, defaults to True
        :type strand: bool, optional
        """
        annotation = annotation.sort_values(by=["fasta_header"])
        annotation.reset_index(inplace=True, drop=True)

        # save the annotation as bed file
        file_bed = os.path.join(self.dir_output, "annotation.bed")
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


class NcbiGenomicRegionGenerator(CustomGenomicRegionGenerator):
    """Class to generate sequences from annotated regions in GTF format and genomic fasta file.
    Sequences are safed as fasta file with region, coordinate and additional information in header.
    GTF and fasta files downloaded from NCBI server via FTP. Taxon, species and annotation release
    are set to default values if not provided.

    :param taxon: Taxon of the files to download, defaults to None
    :type taxon: str, optional
    :param species: Species of the files to dowload, defaults to None
    :type species: str, optional
    :param annotation_release: Annotation release of the files to dowload, defaults to None
    :type annotation_release: str, optional
    :param dir_output: directory where the files are saved, defaults to './output/annotation'
    :type dir_output: str, optional
    """

    def __init__(
        self,
        taxon: str = None,
        species: str = None,
        annotation_release: str = None,
        dir_output: str = "output",
    ):
        """Constructor"""
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
    """Class to generate sequences from annotated regions in GTF format and genomic fasta file.
    Sequences are safed as fasta file with region, coordinate and additional information in header.
    GTF and fasta files downloaded from Ensemble server via FTP. Species and annotation release
    are set to default values if not provided.

    :param species: Species of the files to dowload, defaults to None
    :type species: str, optional
    :param annotation_release: Annotation release of the files to dowload, defaults to None
    :type annotation_release: str, optional
    :param dir_output: directory where the files are saved, defaults to './output/annotation'
    :type dir_output: str, optional
    """

    def __init__(
        self,
        species: str = None,
        annotation_release: str = None,
        dir_output: str = "output",
    ):
        """Constructor"""
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
