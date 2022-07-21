import os

import IO._data_parser as data_parser
import pandas as pd
from Bio import SeqIO


class Oligos:
    """Creates the oligo sequences belonging to the given genes and filters them according to a list of filters given in input"""

    def __init__(self, probe_length_min, probe_length_max, file_sequence, filters):
        """Initialize the class

        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int_
        :param file_sequence: path to the fasta file of the whole genome
        :type file_sequence: str
        :param filters: list of filters classes already initialized
        :type filters: list of classes
        """
        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max
        self.file_sequence = file_sequence
        self.filters = filters

    def generate(self, file_region_annotation, genes, number_batchs, dir_output):
        """Get the fasta sequence of all possible probes with user-defined length for all input genes and store the in a dictionary.
        Generated probes are filtered by the filters give in input and the features computed for the  filters are added to the dictionary.

        :param file_region_annotation: path to the gtf annotaiton file of the region
        :type file_region_annotation: str
        :param genes: genes for which the probes are computed
        :type genes: list of str
        :param number_batchs: probes are computed for a batch of genes
        :type number_batchs: int
        :param dir_output: directory where temporary diles are written, defaults to './output/annotation'
        :type dir_output: str, optional
        :return: the oligos DB containing all the probes created
        :rtype: dict
        """

        def _get_probes(batch_id, genes_batch):
            """Get the fasta sequence of all possible probes for all genes in the batch.

            :param batch_id: Batch ID.
            :type batch_id: int
            :param genes_batch: List of genes for which probes should be designed.
            :type genes_batch: list
            """

            file_region_bed_batch = os.path.join(
                dir_output, "region_batch{}.bed".format(batch_id)
            )
            file_region_fasta_batch = os.path.join(
                dir_output, "region_batch{}.fna".format(batch_id)
            )

            _get_region_fasta(
                genes_batch, file_region_bed_batch, file_region_fasta_batch
            )
            gene_probes_batch = _get_probes_info(genes_batch, file_region_fasta_batch)

            os.remove(file_region_bed_batch)
            os.remove(file_region_fasta_batch)

            return gene_probes_batch

        def _get_region_fasta(
            genes_batch, file_region_bed_batch, file_region_fasta_batch
        ):
            """Extract transcripts for the current batch and write transcript regions to bed file.
            Get sequence for annotated transcript regions (bed file) from genome sequence (fasta file) and write transcriptome sequences to fasta file.

            :param genes_batch: List of genes for which probes should be designed.
            :type genes_batch: list
            :param file_region_bed_batch: Path to bed transcriptome annotation output file.
            :type file_region_bed_batch: string
            :param file_region_fasta_batch: Path to fasta transcriptome sequence output file.
            :type file_region_fasta_batch: string
            """

            region_annotation_batch = region_annotation.loc[
                region_annotation["gene_id"].isin(genes_batch)
            ].copy()
            region_annotation_batch = region_annotation_batch.sort_values(
                by=["gene_id"]
            )
            region_annotation_batch[
                [
                    "seqname",
                    "start",
                    "end",
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
            ].to_csv(file_region_bed_batch, sep="\t", header=False, index=False)

            # get sequence for exons
            data_parser.get_sequence_from_annotation(
                file_region_bed_batch,
                self.file_sequence,
                file_region_fasta_batch,
                split=True,
            )

        def _parse_header(header):
            """Helper function to parse the header of each exon fasta entry.

            :param header: Header of fasta entry.
            :type header: string
            :return: Parsed header information, i.e.
                gene_id, transcript_id, exon_id and position (chromosome, start, strand)
            :rtype: string
            """

            identifier = header.split("::")[0]
            gene_id = identifier.split("_tid")[0]
            transcript_id = identifier.split("_tid")[1].split("_eid")[0]
            exon_id = identifier.split("_eid")[1]

            coordinates = header.split("::")[1]
            chrom = coordinates.split(":")[0]
            start = int(coordinates.split(":")[1].split("-")[0])
            strand = coordinates.split("(")[1].split(")")[0]

            return gene_id, transcript_id, exon_id, chrom, start, strand

        def _get_probes_info(genes_batch, file_region_fasta_batch):
            """Merge all probes with identical sequence that come from the same gene into one fasta entry.
            Filter all probes based on user-defined filters and collect the additional information about each probe.

            :param genes_batch: List of genes for which probes should be designed.
            :type genes_batch: list
            :param file_region_fasta_batch: Path to fasta transcriptome sequence output file.
            :type file_region_fasta_batch: string
            :return: Mapping of probes to corresponding genes with additional information about each probe, i.e.
                position (chromosome, start, end, strand), gene_id, transcript_id, exon_id
            :rtype: dict
            """

            gene_probes = {key: {} for key in genes_batch}
            total_probes = 0

            # parse the exon fasta sequence file
            for exon in SeqIO.parse(file_region_fasta_batch, "fasta"):
                sequence = exon.seq

                for probe_length in range(
                    self.probe_length_min, self.probe_length_max + 1
                ):
                    if len(sequence) > probe_length:
                        number_probes = len(sequence) - (probe_length - 1)
                        probes_sequence = [
                            sequence[i : i + probe_length] for i in range(number_probes)
                        ]

                        for i in range(number_probes):
                            total_probes += 1
                            probe_sequence = probes_sequence[i]

                            (
                                gene_id,
                                transcript_id,
                                exon_id,
                                chrom,
                                start,
                                strand,
                            ) = _parse_header(exon.id)
                            probe_start = start + i
                            probe_end = start + i + probe_length

                            if probe_sequence in gene_probes[gene_id]:
                                gene_probes[gene_id][probe_sequence][
                                    "transcript_id"
                                ].append(transcript_id)
                                gene_probes[gene_id][probe_sequence]["exon_id"].append(
                                    exon_id
                                )
                                gene_probes[gene_id][probe_sequence]["start"].append(
                                    probe_start
                                )
                                gene_probes[gene_id][probe_sequence]["end"].append(
                                    probe_end
                                )

                            else:
                                gene_probes[gene_id][probe_sequence] = {
                                    "transcript_id": [transcript_id],
                                    "exon_id": [exon_id],
                                    "chromosome": chrom,
                                    "start": [probe_start],
                                    "end": [probe_end],
                                    "strand": strand,
                                    "length": probe_length,
                                }
            return gene_probes

        region_annotation = pd.read_csv(
            file_region_annotation,
            sep="\t",
            names=[
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
            ],
        )

        # keep the batch structure for parellalization
        batch_size = int(len(genes) / number_batchs) + (len(genes) % number_batchs > 0)
        probes = {}
        for batch_id in range(number_batchs):
            genes_batch = genes[
                (batch_size * batch_id) : (
                    min(batch_size * (batch_id + 1), len(genes) + 1)
                )
            ]
            probes.update(_get_probes(batch_id, genes_batch))
        probes = self.filter_oligos_DB(probes)

        return probes

    def filter(self, sequence):
        """Applies the used-defined filters and returns the result and the additional computed features

        :param sequence: sequence to check
        :type sequence: Bio.Seq
        :return: if the filters are fulfilled and the additional features computed
        :rtype: bool, dict
        """

        fulfills = True
        additional_features = {}
        for filt in self.filters:
            out, feat = filt.apply(sequence)
            if not out:  # stop at the first false we obtain
                return False, {}
            additional_features.update(feat)
        return fulfills, additional_features

    def filter_oligos_DB(self, oligos_DB):
        """Filter the database of probes based on teh given filters

        :param oligos_DB: da5tabase of probes
        :type oligos_DB: dict
        :return: filtered oligos DB
        :rtype: dict
        """

        loaded_probes = 0
        gene_ids = list(oligos_DB.keys())
        for gene_id in gene_ids:
            probes_sequences = list(oligos_DB[gene_id].keys())
            for probe_sequence in probes_sequences:
                fulfills, additional_features = self.filter(probe_sequence)
                if fulfills:
                    oligos_DB[gene_id][probe_sequence].update(additional_features)
                    loaded_probes += 1
                else:
                    del oligos_DB[gene_id][probe_sequence]
        return oligos_DB
