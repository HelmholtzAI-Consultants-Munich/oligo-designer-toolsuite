import os

import IO._data_parser as data_parser
import pandas as pd


class GeneTranscript:
    """Class that creates the transcriptome for the whole genome or for a set of genes."""

    def __init__(self, file_sequence, file_annotation):
        self.file_sequence = file_sequence
        self.file_annotation = file_annotation
        self.annotation = data_parser.read_gtf(
            self.file_annotation
        )  # dataframe with annotation file

    def generate_for_reference(self, block_size, dir_output="./output/annotation"):
        # does also teh merging ad generates the fasta file
        # returns the path of the annotation and fasta files

        # get annotation of exons and
        exons_annotation = self.__load_exon_annotation()
        # merge exon annotations for the same region
        unique_exons = self.__load_unique_exons(exons_annotation)
        unique_exons = self.__merge_containing_exons(unique_exons)
        # get exon junction annotation
        exon_junctions_probes = self.__load_exon_junctions(block_size, exons_annotation)
        # concatenate the two annotations
        gene_transcript_annotation = pd.concat([unique_exons, exon_junctions_probes])
        gene_transcript_annotation = gene_transcript_annotation.sort_values(
            by=["gene_id"]
        )
        gene_transcript_annotation.reset_index(inplace=True, drop=True)
        # save the gene transcript annotation
        file_gene_transcript_annotation = os.path.join(
            dir_output, "gene_transcript.bed"
        )
        gene_transcript_annotation[
            [
                "seqname",
                "start",
                "end",
                "gene_id",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "block_sizes",
                "blockStarts",
            ]
        ].to_csv(file_gene_transcript_annotation, sep="\t", header=False, index=False)
        # create the fasta file
        file_gene_transcript_fasta = os.path.join(dir_output, "gene_transcript.fasta")
        data_parser.get_sequence_from_annotation(
            file_gene_transcript_annotation,
            self.file_sequence,
            file_gene_transcript_fasta,
            split=True,
        )

        return file_gene_transcript_annotation, file_gene_transcript_fasta

    def generate_for_oligos(
        self, block_size, dir_output="./output/annotation", genes=None
    ):
        # can take in input a list of genes for witch the transcriptome will be generated
        # returns the path of the annotation file
        if genes is None:
            annotation_genes = None
        else:
            annotation_genes = self.annotation.loc[
                self.annotation["gene_id"].isin(genes)
            ].copy()
        # get annotation of exons and
        exons_annotation = self.__load_exon_annotation(annotation_genes)
        # merge exon annotations for the same region
        unique_exons = self.__load_unique_exons(exons_annotation)
        # get exon junction annotation
        exon_junctions_probes = self.__load_exon_junctions(
            block_size, exons_annotation, annotation_genes
        )
        # concatenate the two annotations
        gene_transcript_annotation = pd.concat([unique_exons, exon_junctions_probes])
        gene_transcript_annotation = gene_transcript_annotation.sort_values(
            by=["gene_id"]
        )
        gene_transcript_annotation.reset_index(inplace=True, drop=True)
        # save the gene transcript annotation
        file_gene_transcript_annotation = os.path.join(
            dir_output, "gene_transcript.bed"
        )
        # for the oligos we don't need the fasta file and we need the gene, transcript, exon informations
        gene_transcript_annotation[
            [
                "seqname",
                "start",
                "end",
                "gene_id",
                "gene_transcript_exon_id",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "block_sizes",
                "blockStarts",
            ]
        ].to_csv(file_gene_transcript_annotation, sep="\t", header=False, index=False)
        return file_gene_transcript_annotation

    def __load_exon_annotation(self, annotation_genes=None):
        """Retrive exon annotation from loaded gene annotation.

        :return: Exon annotation
        :rtype: pandas.DataFrame
        """
        if annotation_genes is None:
            exon_annotation = self.annotation.loc[self.annotation["feature"] == "exon"]
        else:
            exon_annotation = annotation_genes.loc[
                annotation_genes["feature"] == "exon"
            ]
        exon_annotation = exon_annotation.assign(source="unknown")

        if not "exon_id" in exon_annotation.columns:
            exon_annotation["exon_id"] = (
                exon_annotation["transcript_id"]
                + "_exon"
                + exon_annotation["exon_number"]
            )

        # there are some exon annotations which have the same start and end coordinates and can't be saved as fasta from bedtools
        exon_annotation = exon_annotation[
            (exon_annotation.end - exon_annotation.start) > 0
        ]
        exon_annotation.reset_index(inplace=True, drop=True)
        exon_annotation["start"] -= 1

        return exon_annotation

    def __load_unique_exons(self, exon_annotation):
        """Merge overlapping exons, which have the same start and end coordinates.
        Save transcript information for those exons.

        :param exon_annotation: Exon annotation
        :type exon_annotation: pandas.DataFrame
        :return: Merged exon annotation
        :rtype: pandas.DataFrame
        """

        exon_annotation["region"] = (
            exon_annotation["seqname"]
            + "_"
            + exon_annotation["start"].astype("str")
            + "_"
            + exon_annotation["end"].astype("str")
            + "_"
            + exon_annotation["strand"]
        )

        aggregate_function = {
            "gene_id": "first",
            "transcript_id": ":".join,
            "exon_id": ":".join,
            "seqname": "first",
            "start": "first",
            "end": "first",
            "score": "first",
            "strand": "first",
        }
        merged_exons = exon_annotation.groupby(exon_annotation["region"]).aggregate(
            aggregate_function
        )

        merged_exons["score"] = 0
        merged_exons["thickStart"] = merged_exons["start"]
        merged_exons["thickEnd"] = merged_exons["end"]
        merged_exons["itemRgb"] = 0
        merged_exons["blockCount"] = 1
        merged_exons["block_sizes"] = merged_exons["end"] - merged_exons["start"]
        merged_exons["blockStarts"] = 0
        merged_exons["gene_transcript_exon_id"] = (
            merged_exons["gene_id"]
            + "_tid"
            + merged_exons["transcript_id"]
            + "_eid"
            + merged_exons["exon_id"]
        )
        merged_exons = merged_exons[
            [
                "gene_id",
                "gene_transcript_exon_id",
                "seqname",
                "start",
                "end",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "block_sizes",
                "blockStarts",
            ]
        ]

        return merged_exons

    def __merge_containing_exons(self, unique_exons):
        """Merge exons that are contained in a larger exon, e.g. have the same start coordinates but different end coordinates, into one entry.

        :param unique_exons: Dataframe with annotation of unique exons, where overlapping exons are merged.
        :type unique_exons: pandas.DataFrame
        :return: Dataframe with annotation of merged exons, where containing exons are merged.
        :rtype: pandas.DataFrame
        """

        aggregate_function = {
            "gene_id": "first",
            "gene_transcript_exon_id": ":".join,
            "seqname": "first",
            "start": "min",
            "end": "max",
            "score": "first",
            "strand": "first",
            "thickStart": "min",
            "thickEnd": "max",
            "itemRgb": "first",
            "blockCount": "first",
            "block_sizes": "max",
            "blockStarts": "first",
        }

        unique_exons["region_start"] = (
            unique_exons["seqname"]
            + "_"
            + unique_exons["start"].astype("str")
            + "_"
            + unique_exons["strand"]
        )
        merged_unique_exons = unique_exons.groupby(
            unique_exons["region_start"]
        ).aggregate(aggregate_function)

        merged_unique_exons["region_end"] = (
            merged_unique_exons["seqname"]
            + "_"
            + merged_unique_exons["end"].astype("str")
            + "_"
            + merged_unique_exons["strand"]
        )
        merged_unique_exons = merged_unique_exons.groupby(
            merged_unique_exons["region_end"]
        ).aggregate(aggregate_function)

        return merged_unique_exons

    def __load_exon_junctions(self, block_size, exon_annotation, annotation_genes=None):
        """Get all possible exon junctions from the transcript annotation.
        Merge overlapping exons jucntions, which have the same satrt and end coordinates.
        Save transcript information for those exons.

        :param block_size: Size of the exon junction regions, i.e. <block_size> bp upstream of first exon
            and <block_size> bp downstream of second exon.
        :type block_size: int
        :return: Dataframe with annotation of exons junctions, where overlapping exons junctions are merged.
        :rtype: pandas.DataFrame
        """
        transcript_exons, transcript_info = self.__get_transcript_exons_and_info(
            exon_annotation, annotation_genes
        )
        exon_junction_list = self.__get_exon_junction_list(
            block_size, transcript_exons, transcript_info
        )

        exon_junctions = pd.DataFrame(
            exon_junction_list,
            columns=[
                "region",
                "gene_id",
                "transcript_id",
                "exon_id",
                "seqname",
                "start",
                "end",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "block_sizes",
                "blockStarts",
            ],
        )
        aggregate_function = {
            "gene_id": "first",
            "transcript_id": ":".join,
            "exon_id": ":".join,
            "seqname": "first",
            "start": "first",
            "end": "first",
            "score": "first",
            "strand": "first",
            "thickStart": "first",
            "thickEnd": "first",
            "itemRgb": "first",
            "blockCount": "first",
            "block_sizes": "first",
            "blockStarts": "first",
        }
        merged_exon_junctions = exon_junctions.groupby(
            exon_junctions["region"]
        ).aggregate(aggregate_function)

        merged_exon_junctions["gene_transcript_exon_id"] = (
            merged_exon_junctions["gene_id"]
            + "_tid"
            + merged_exon_junctions["transcript_id"]
            + "_eid"
            + merged_exon_junctions["exon_id"]
        )
        merged_exon_junctions = merged_exon_junctions[
            [
                "gene_id",
                "gene_transcript_exon_id",
                "seqname",
                "start",
                "end",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "block_sizes",
                "blockStarts",
            ]
        ]

        return merged_exon_junctions

    def __get_transcript_exons_and_info(self, exon_annotation, annotation_genes=None):
        """Get list of exons that belong to a transcript.
        Save information about each transcript, i.e. gene ID, chromosome and strand.

        :param exon_annotation: Exon annotation.
        :type exon_annotation: pandas.DataFrame
        :return: Two dictionaries: transcript - exon mapping and transcript information.
        :rtype: dict
        """
        if annotation_genes is None:
            transcripts = self.annotation.loc[self.annotation["transcript_id"] != ""]
        else:
            transcripts = annotation_genes.loc[annotation_genes["transcript_id"] != ""]
        transcripts = sorted(list(transcripts["transcript_id"].unique()))
        transcript_exons = {key: {} for key in transcripts}
        transcript_info = dict()

        for (
            transcript_id,
            exon_number,
            exon_id,
            start,
            end,
            gene_id,
            seqname,
            strand,
        ) in zip(
            exon_annotation["transcript_id"],
            exon_annotation["exon_number"],
            exon_annotation["exon_id"],
            exon_annotation["start"],
            exon_annotation["end"],
            exon_annotation["gene_id"],
            exon_annotation["seqname"],
            exon_annotation["strand"],
        ):
            transcript_exons[transcript_id][int(exon_number)] = [
                exon_id,
                start,
                end,
            ]
            transcript_info[transcript_id] = [gene_id, seqname, strand]

        if not ("transcript" in self.annotation["feature"].values):
            delete_ids = []
            for transcript_id in transcript_exons:
                if transcript_id not in transcript_info:
                    delete_ids.append(transcript_id)
            for transcript_id in delete_ids:
                del transcript_exons[transcript_id]

        return transcript_exons, transcript_info

    def __get_exon_junction_list(self, block_size, transcript_exons, transcript_info):
        """Get list of all possible exon junctions and save all required information for bed12 file.

        :param block_size: Size of the exon junction regions, i.e. <block_size> bp upstream of first exon
            and <block_size> bp downstream of second exon.
        :type block_size: int
        :param transcript_exons: Transcript - exon mapping.
        :type transcript_exons: dict
        :param transcript_info: Transcript information.
        :type transcript_info: dict
        :return: List of exon junctions.
        :rtype: list
        """
        exon_junction_list = []

        for transcript, exons in transcript_exons.items():
            gene_id = transcript_info[transcript][0]
            seqname = transcript_info[transcript][1]
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
                if idx == 0:  # first exon of transcript
                    exon_upstream = attributes
                    exons_small = []
                elif ((idx + 1) < len(exons)) & (
                    (attributes[2] - attributes[1]) < block_size
                ):  # if exon is not the last exon of transcript and shorter than probe block_size but not the last exon -> create sequence with neighboring exons
                    exons_small.append(attributes)
                else:
                    exon_downstream = attributes
                    block_size_up = min(
                        block_size, (exon_upstream[2] - exon_upstream[1])
                    )  # catch case that first or last exon < block_size
                    block_size_down = min(
                        block_size, (exon_downstream[2] - exon_downstream[1])
                    )  # catch case that first or last exon < block_size
                    start_up = exon_upstream[2] - block_size_up
                    end_down = exon_downstream[1] + block_size_down

                    if exons_small == []:
                        blockCount = 2
                        block_size_length_entry = "{},{}".format(
                            block_size_up, block_size_down
                        )
                        block_size_start_entry = "{},{}".format(
                            0, exon_downstream[1] - start_up
                        )
                    else:
                        blockCount = len(exons_small) + 2
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
                                    str(
                                        attributes[1] - start_up,
                                    )
                                    for attributes in exons_small
                                ]
                            )
                            + ","
                            + str(exon_downstream[1] - start_up)
                        )
                        exons_small = []

                    exon_junction_list.append(
                        [
                            "{}_{}_{}_{}".format(seqname, start_up, end_down, strand),
                            gene_id,
                            transcript,
                            "{}_{}".format(exon_upstream[0], exon_downstream[0]),
                            seqname,
                            start_up,
                            end_down,
                            0,
                            strand,
                            start_up,
                            end_down,
                            0,
                            blockCount,
                            block_size_length_entry,
                            block_size_start_entry,
                        ]
                    )
                    exon_upstream = attributes

        return exon_junction_list
