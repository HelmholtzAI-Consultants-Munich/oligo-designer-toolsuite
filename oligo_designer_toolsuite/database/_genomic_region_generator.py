############################################
# imports
############################################

import os
import copy
import shutil
import warnings
import pandas as pd

pd.options.mode.chained_assignment = None

from pathlib import Path

from ..utils._data_parser import get_sequence_from_annotation
from ..utils._ftp_loader import FtpLoaderEnsembl, FtpLoaderNCBI
from ..utils._gff_parser import GffParser


############################################
# Genomic Region Generator Classes
############################################


class CustomGenomicRegionGenerator:
    """Class to generate sequences from annotated regions in GTF format and genomic fasta file.
    Sequences are saved as fasta file with region id, additional information and coordinates in header.

    Output Format (per sequence):
    >'region_id'::'additional information'::'chromosome':'start'-'end'('strand')
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    :param annotation_file: Annotation of genomic regions. Has to be provided in GTF format.
    :type annotation_file: str
    :param sequence_file: Fasta file with genome sequence.
    :type sequence_file: str
    :param source: Source of annotations, e.g. NCBI, defaults to None.
    :type source: str, optional
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
        source: str = None,
        species: str = None,
        annotation_release: str = None,
        genome_assembly: str = None,
        dir_output: str = "output",
    ):
        """Constructor"""
        if source is None:
            source = "custom"
            warnings.warn(f"No source defined. Using default source {source}!")

        if species is None:
            species = "unknown"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "unknown"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        if genome_assembly is None:
            genome_assembly = "unknown"
            warnings.warn(
                f"No genome assembly defined. Using default genome assembly {genome_assembly}!"
            )

        self.source = source
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
        self.annotation.start = self.annotation.start.astype("int")
        self.annotation.end = self.annotation.end.astype("int")

        # columns required for bed12 split sequence format
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
        self.FASTA_HEDAER = f"source={self.source};species={self.species};annotation_release={self.annotation_release};genome_assemly={self.genome_assembly}"

    def generate_genome(self):
        """Generate a fasta file with the genome sequences.

        :return: Name of the output fasta file, where the generated transcript sequences are stored.
        :rtype: str
        """
        file_genome_fasta = self._get_fasta_file_name("genome")

        shutil.copyfile(self.sequence_file, file_genome_fasta)
        self._add_header_to_fasta_file(file_genome_fasta)
        return file_genome_fasta

    def generate_transcript_reduced_representation(
        self,
        include_exon_junctions: bool = True,
        exon_junction_size: int = 100,
    ):
        """Generate a reduced representation of the transcriptome. This representation merges overlapping exons from the same
        gene but different transcripts into one sequence entry, to reduce duplicated sequence records in the fasta file.
        In addition to exon, it is also possible to include exon junction information. The "exon_junction_size" parameter
        defines the size of the exon junction region, i.e. +/- "exon_junction_size" bp around the junction.
        The annotation is taken from the classes GTF file, which must include the value "exon" in the 3rd column and the following
        fileds in the attributes column: "gene_id", "transcript_id" and "exon_number".

        :param include_exon_junctions: Include exon junctions in the transcriptome sequences, defaults to True
        :type include_exon_junctions: bool, optional
        :param exon_junction_size: size of the exon junction region, i.e. +/- "exon_junction_size" bp around the junction , defaults to 100
        :type exon_junction_size: int, optional
        :return: Name of the output fasta file, where the generated transcript sequences are stored.
        :rtype: str
        """
        file_transcriptome_fasta = self._get_fasta_file_name(
            "transcriptome", include_exon_junctions, exon_junction_size
        )

        # get annotation of exons
        exons = self._get_annotation_region_of_interest("exon")

        # get exon and (optionally) exon junction sequences
        self._get_exon_sequence(
            exons,
            file_transcriptome_fasta,
            include_exon_junctions,
            exon_junction_size,
        )
        self._add_header_to_fasta_file(file_transcriptome_fasta)
        return file_transcriptome_fasta

    def generate_CDS_reduced_representation(
        self,
        include_exon_junctions: bool = True,
        exon_junction_size: int = 100,
    ):
        """Generate a reduced representation of the coding region (CDS). This representation merges overlapping exons from the same
        gene but different transcripts into one sequence entry, to reduce duplicated sequence records in the fasta file.
        In addition to exon, it is also possible to include exon junction information. The "exon_junction_size" parameter
        defines the size of the exon junction region, i.e. +/- "exon_junction_size" bp around the junction.
        The annotation is taken from the classes GTF file, which must include the value "CDS" in the 3rd column and the following
        fileds in the attributes column: "gene_id", "transcript_id" and "exon_number".

        :param include_exon_junctions: Include exon junctions in the coding region (CDS) sequences, defaults to True
        :type include_exon_junctions: bool, optional
        :param exon_junction_size: size of the exon junction region, i.e. +/- "exon_junction_size" bp around the junction , defaults to 100
        :type exon_junction_size: int, optional
        :return: Name of the output fasta file, where the generated coding region sequences are stored.
        :rtype: str
        """
        file_CDS_fasta = self._get_fasta_file_name(
            "CDS", include_exon_junctions, exon_junction_size
        )

        # get annotation of exons
        exons = self._get_annotation_region_of_interest("CDS")

        self._get_exon_sequence(
            exons,
            file_CDS_fasta,
            include_exon_junctions,
            exon_junction_size,
        )
        self._add_header_to_fasta_file(file_CDS_fasta)
        return file_CDS_fasta

    def _get_fasta_file_name(
        self, region, include_exon_junctions=False, exon_junction_size=None
    ):
        """Get name of output fasta file.

        :param region: Region to extract, which is the 3rd "type" column in the GFT file.
        :type region: str {'exon', 'CDS', 'genome'}
        :param include_exon_junctions: Include exon junctions in the coding region (CDS) sequences, defaults to True
        :type include_exon_junctions: bool, optional
        :param exon_junction_size: size of the exon junction region, i.e. +/- "exon_junction_size" bp around the junction , defaults to None
        :type exon_junction_size: int, optional
        :return: Fasta file name.
        :rtype: str
        """
        FASTA_HEDAER = self.FASTA_HEDAER
        FASTA_HEDAER = FASTA_HEDAER.replace("=", "_").replace(";", "_")
        file_fasta = f"{region}_{FASTA_HEDAER}"
        if include_exon_junctions:
            file_fasta += f"_incl_exonjunctions_of_size_{exon_junction_size}"
        file_fasta += ".fna"
        file_fasta = os.path.join(self.dir_output, file_fasta)
        return file_fasta

    def _add_header_to_fasta_file(self, file_fasta):
        """Add additional information to fasta header comment line.

        :param file_fasta: Fasta file name.
        :type file_fasta: str
        """
        with open(file_fasta, "r") as original:
            data = original.read()
        with open(file_fasta, "w") as modified:
            modified.write(f"# {self.FASTA_HEDAER}\n" + data)

    def _get_annotation_region_of_interest(self, region):
        """Retrieve exon annotation from loaded annotation dataframe.

        :param region: Region to extract, which is the 3rd "type" column in the GFT file.
        :type region: str {'exon', 'CDS'}
        :return: Exon annotation.
        :rtype: pandas.DataFrame
        """
        annotation = copy.deepcopy(self.annotation)
        exons = annotation.loc[annotation["type"] == region]
        exons.reset_index(inplace=True, drop=True)
        return exons

    def _get_exon_sequence(
        self,
        exons,
        file_transcriptome_fasta,
        include_exon_junctions,
        exon_junction_size,
    ):
        """Get sequence in fasta format for exons of transcriptome or CDS. Retrieve fasta file from
        generated bed12 and genome fasta file.

        :param exons: Exon annotation.
        :type exons: pandas.DataFrame
        :param file_transcriptome_fasta: Name of the output fasta file, where the generated coding region sequences are stored.
        :type file_transcriptome_fasta: str
        :param include_exon_junctions: Include exon junctions, defaults to True
        :type include_exon_junctions: bool, optional
        :param exon_junction_size: size of the exon junction region, i.e. +/- "exon_junction_size" bp around the junction , defaults to 100
        :type exon_junction_size: int, optional
        """
        exons.start = exons.start - 1  # convert GFF 1-base offset to BED 0-base offset
        list_annotations = []
        # merge exon annotations for the same region
        unique_exons = self._get_unique_exons(exons)
        list_annotations.append(unique_exons)
        # merges exons with with same start but different end or different start but same end coordinates
        # since oligos are also created from this transcriptome, disable this function
        # unique_exons = self._merge_containing_exons(unique_exons)

        # get exon junction annotation and merge annotation for the same region
        if include_exon_junctions:
            exon_junctions = self._get_exon_junction_annotation(
                exons, exon_junction_size
            )
            unique_exon_junctions = self._get_unique_exon_junctions(exon_junctions)
            list_annotations.append(unique_exon_junctions)

        # concatenate the two annotations
        annotation = pd.concat(list_annotations)
        annotation = annotation.sort_values(by=["fasta_header"])
        annotation.reset_index(inplace=True, drop=True)

        # save the gene transcript annotation
        file_transcriptome_bed = os.path.join(self.dir_output, "transcriptome.bed")
        annotation[self.BED12_HEADER].to_csv(
            file_transcriptome_bed, sep="\t", header=False, index=False
        )

        # create the fasta file
        get_sequence_from_annotation(
            file_transcriptome_bed,
            self.sequence_file,
            file_transcriptome_fasta,
            split=True,
            strand=True,
            nameOnly=True,
        )
        # os.remove(file_transcriptome_bed)

    def _get_unique_exons(self, exons):
        """Merge overlapping exons, which have the same start and end coordinates.
        Save transcript information for those exons.

        :param exons: Exon annotation
        :type exons: pandas.DataFrame
        :return: Merged exon annotation
        :rtype: pandas.DataFrame
        """
        exons["region"] = (
            exons["seqid"]
            + ":"
            + exons["start"].astype("str")
            + "-"
            + exons["end"].astype("str")
            + "("
            + exons["strand"]
            + ")"
        )

        exons["transcript_exon_id"] = (
            "transcript_id="
            + exons["transcript_id"]
            + ",exon_number="
            + exons["exon_number"]
        )
        aggregate_function = {
            "region": "first",
            "gene_id": "first",
            "transcript_exon_id": ";".join,
            "seqid": "first",
            "start": "first",
            "end": "first",
            "score": "first",
            "strand": "first",
        }
        merged_exons = exons.groupby(exons["region"]).aggregate(aggregate_function)

        merged_exons["score"] = 0
        merged_exons["thickStart"] = merged_exons["start"]
        merged_exons["thickEnd"] = merged_exons["end"]
        merged_exons["itemRgb"] = 0
        merged_exons["block_count"] = 1
        merged_exons["block_sizes"] = merged_exons["end"] - merged_exons["start"]
        merged_exons["blockStarts"] = 0
        merged_exons["fasta_header"] = (
            merged_exons["gene_id"]
            + "::"
            + merged_exons["transcript_exon_id"]
            + "::"
            + merged_exons["region"]
        )
        merged_exons = merged_exons[self.BED12_HEADER]

        return merged_exons

    def _get_exon_junction_annotation(
        self,
        exons,
        block_size,
    ):
        """Get list of all possible exon junctions and save all required information for bed12 file.

        :param exons: Exon annotation
        :type exons: pandas.DataFrame
        :param block_size: Size of the exon junction regions, i.e. <block_size> bp upstream of first exon
            and <block_size> bp downstream of second exon.
        :type block_size: int
        :return: Dataframe with exon junctions.
        :rtype: pandas.DataFrame
        """
        transcript_exons, transcript_info = self._get_transcript_exons_and_info(exons)
        exon_junction_list = []

        for transcript, exons in transcript_exons.items():
            gene_id = transcript_info[transcript][0]
            seqid = transcript_info[transcript][1]
            strand = transcript_info[transcript][2]

            if strand == "+":
                exons = [
                    entry[1] for entry in sorted(exons.items())
                ]  # return only exon attributes sorted by exon number
            elif strand == "-":
                exons = [
                    entry[1] for entry in sorted(exons.items(), reverse=True)
                ]  # return only exon attributes sorted by exon number -> reverse sort on minus strand

            for idx, attributes in enumerate(exons):
                # first exon of transcript
                if idx == 0:
                    exon_upstream = attributes
                    exons_small = []
                    regions_exons_small = ""

                # if exon is not the last exon of transcript and shorter than oligo block_size but not the last exon -> create sequence with neighboring exons
                elif ((idx + 1) < len(exons)) & (
                    (attributes[2] - attributes[1]) < block_size
                ):
                    exons_small.append(attributes)

                else:
                    exon_downstream = attributes
                    # catch case that first or last exon < block_size
                    block_size_up = min(
                        block_size, (exon_upstream[2] - exon_upstream[1])
                    )
                    block_size_down = min(
                        block_size, (exon_downstream[2] - exon_downstream[1])
                    )
                    start_up = exon_upstream[2] - block_size_up
                    end_down = exon_downstream[1] + block_size_down

                    if exons_small == []:
                        block_count = 2
                        block_size_length_entry = "{},{}".format(
                            block_size_up, block_size_down
                        )
                        block_size_start_entry = "{},{}".format(
                            0, exon_downstream[1] - start_up
                        )

                    else:
                        block_count = len(exons_small) + 2
                        block_size_length_entry = (
                            str(block_size_up)
                            + ","
                            + ",".join(
                                [
                                    str(attributes[2] - attributes[1])
                                    for attributes in exons_small
                                ]
                            )
                            + ","
                            + str(block_size_down)
                        )
                        block_size_start_entry = (
                            "0,"
                            + ",".join(
                                [
                                    str(attributes[1] - start_up)
                                    for attributes in exons_small
                                ]
                            )
                            + ","
                            + str(exon_downstream[1] - start_up)
                        )
                        regions_exons_small = ";".join(
                            [
                                f"{seqid}:{attributes[1]}-{attributes[2]}({strand})"
                                for attributes in exons_small
                            ]
                        )

                    region_up = f"{seqid}:{start_up}-{start_up+block_size_up}({strand})"
                    region_down = (
                        f"{seqid}:{end_down-block_size_down}-{end_down}({strand})"
                    )

                    exon_junction_list.append(
                        [
                            ";".join(
                                filter(
                                    None, [region_up, regions_exons_small, region_down]
                                )
                            ),
                            gene_id,
                            f"transcript_id={transcript},exon_number={exon_upstream[0]}__JUNC__{exon_downstream[0]}",
                            seqid,
                            start_up,
                            end_down,
                            0,
                            strand,
                            start_up,
                            end_down,
                            0,
                            block_count,
                            block_size_length_entry,
                            block_size_start_entry,
                        ]
                    )
                    exons_small = []
                    regions_exons_small = ""
                    exon_upstream = attributes

            exon_junctions = pd.DataFrame(
                exon_junction_list,
                columns=[
                    "region",
                    "gene_id",
                    "transcript_exon_id",
                    "seqid",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "thickStart",
                    "thickEnd",
                    "itemRgb",
                    "block_count",
                    "block_sizes",
                    "blockStarts",
                ],
            )

        return exon_junctions

    def _get_transcript_exons_and_info(self, exons):
        """Get list of exons that belong to a transcript.
        Save information about each transcript, i.e. gene ID, chromosome and strand.

        :param exons: Exon annotation.
        :type exons: pandas.DataFrame
        :return: dict with transcript - exon mapping, dict with transcript information.
        :rtype: dict, dict
        """
        transcript_ids = exons.loc[exons["transcript_id"].notna()]
        transcript_ids = sorted(list(transcript_ids["transcript_id"].unique()))
        transcript_exons = {key: {} for key in transcript_ids}
        transcript_info = dict()

        for (
            transcript_id,
            exon_number,
            start,
            end,
            gene_id,
            seqid,
            strand,
        ) in zip(
            exons["transcript_id"],
            exons["exon_number"],
            exons["start"],
            exons["end"],
            exons["gene_id"],
            exons["seqid"],
            exons["strand"],
        ):
            transcript_exons[transcript_id][int(exon_number)] = [
                exon_number,
                start,
                end,
            ]
            transcript_info[transcript_id] = [gene_id, seqid, strand]

        if not ("transcript" in self.annotation["type"].values):
            delete_ids = []
            for transcript_id in transcript_exons:
                if transcript_id not in transcript_info:
                    delete_ids.append(transcript_id)
            for transcript_id in delete_ids:
                del transcript_exons[transcript_id]

        return transcript_exons, transcript_info

    def _get_unique_exon_junctions(self, exon_junctions):
        """Get all possible exon junctions from the transcript annotation.
        Merge overlapping exons jucntions, which have the same satrt and end coordinates.
        Save transcript information for those exons.

        :param exon_junctions: Dataframe with exon junctions.
        :type exon_junctions: pandas.DataFrame
        :return: Dataframe with annotation of exons junctions, where overlapping exons junctions are merged.
        :rtype: pandas.DataFrame
        """
        aggregate_function = {
            "region": "first",
            "gene_id": "first",
            "transcript_exon_id": ";".join,
            "seqid": "first",
            "start": "first",
            "end": "first",
            "score": "first",
            "strand": "first",
            "thickStart": "first",
            "thickEnd": "first",
            "itemRgb": "first",
            "block_count": "first",
            "block_sizes": "first",
            "blockStarts": "first",
        }
        merged_exon_junctions = exon_junctions.groupby(
            exon_junctions["region"]
        ).aggregate(aggregate_function)

        merged_exon_junctions["fasta_header"] = (
            merged_exon_junctions["gene_id"]
            + "::"
            + merged_exon_junctions["transcript_exon_id"]
            + "::"
            + merged_exon_junctions["region"]
        )
        merged_exon_junctions = merged_exon_junctions[self.BED12_HEADER]

        return merged_exon_junctions


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
        dir_output: str = "output/annotation",
    ):
        """Constructor"""
        source = "NCBI"
        if taxon is None:
            taxon = "vertebrate_mammalian"
            warnings.warn(f"No taxon defined. Using default taxon {taxon}!")

        if species is None:
            species = "Homo_sapiens"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        Path(dir_output).mkdir(parents=True, exist_ok=True)

        ftp = FtpLoaderNCBI(dir_output, taxon, species, annotation_release)
        annotation_file, annotation_release, genome_assembly = ftp.download_files("gtf")
        sequence_file, _, _ = ftp.download_files("fasta")

        super().__init__(
            annotation_file,
            sequence_file,
            source,
            species,
            annotation_release,
            genome_assembly,
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
        dir_output: str = "output/annotation",
    ):
        """Constructor"""
        source = "Ensemble"
        if species is None:
            species = "homo_sapiens"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        Path(dir_output).mkdir(parents=True, exist_ok=True)

        ftp = FtpLoaderEnsembl(dir_output, species, annotation_release)
        annotation_file, annotation_release, genome_assembly = ftp.download_files("gtf")
        sequence_file, _, _ = ftp.download_files("fasta")

        super().__init__(
            annotation_file,
            sequence_file,
            source,
            species,
            annotation_release,
            genome_assembly,
        )
