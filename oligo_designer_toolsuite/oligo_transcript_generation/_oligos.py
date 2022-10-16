import os

import pandas as pd
from Bio import SeqIO
from joblib import Parallel, delayed

import oligo_designer_toolsuite.IO._data_parser as _data_parser


class Oligos:
    """Creates the oligo sequences belonging to the given genes and filters them according to a list of filters given in input"""

    def __init__(self, probe_length_min, probe_length_max, file_sequence):
        """Initialize the class

        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int_
        :param file_sequence: path to the fasta file of the whole genome
        :type file_sequence: str

        """
        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max
        self.file_sequence = file_sequence

    def generate(self, file_region_annotation, genes, n_jobs, dir_output):
        """Get the fasta sequence of all possible probes with user-defined length for all input genes and store the in a dictionary.
        Generated probes are filtered by the filters give in input and the features computed for the  filters are added to the dictionary.

        :param file_region_annotation: path to the gtf annotaiton file of the region
        :type file_region_annotation: str
        :param genes: Genes for which the probes are computed
        :type genes: list of str
        :param n_jobs: nr of cores used to compute the probes
        :type n_jobs: int
        :param dir_output: directory where temporary diles are written, defaults to './output/annotation'
        :type dir_output: str, optional
        :return: the oligos DB containing all the probes created
        :rtype: dict
        """

        def _get_probes(gene, region_annotation):
            """Get the fasta sequence of all possible probes for all genes in the batch.

            :param gene: gene for which probes should be designed.
            :type gene: str
            :param region_annotation: gtf annotation of the region
            :type region_annotation: pandas.DataFrame
            :return: the oligos DB containing all the probes created for the gene
            :rtype: dict
            """

            file_region_bed = os.path.join(dir_output, "region_gene{}.bed".format(gene))
            file_region_fasta = os.path.join(
                dir_output, "region_gene{}.fna".format(gene)
            )

            _get_region_fasta(
                gene, region_annotation, file_region_bed, file_region_fasta
            )
            gene_probes_batch = _get_probes_info(gene, file_region_fasta)

            os.remove(file_region_bed)
            os.remove(file_region_fasta)

            return gene_probes_batch

        def _get_region_fasta(
            gene, region_annotation, file_region_bed, file_region_fasta
        ):
            """Extract transcripts for the current batch and write transcript regions to bed file.
            Get sequence for annotated transcript regions (bed file) from genome sequence (fasta file) and write transcriptome sequences to fasta file.

            :param genes_batch: Genes for which probes should be designed.
            :type genes_batch: list
            :param region_annotation: gtf annotation of the region
            :type region_annotation: pandas.DataFrame
            :param file_region_bed: Path to bed transcriptome annotation output file.
            :type file_region_bed: string
            :param file_region_fasta: Path to fasta transcriptome sequence output file.
            :type file_region_fasta: string
            """

            region_annotation_gene = region_annotation.loc[
                region_annotation["gene_id"] == gene
            ].copy()
            region_annotation_gene[
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
            ].to_csv(file_region_bed, sep="\t", header=False, index=False)

            # get sequence for exons
            _data_parser.get_sequence_from_annotation(
                file_region_bed,
                self.file_sequence,
                file_region_fasta,
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

        def _get_probes_info(gene, file_region_fasta):
            """Merge all probes with identical sequence that come from the same gene into one fasta entry.
            Filter all probes based on user-defined filters and collect the additional information about each probe.

            :param gene: Gene for which probes should be designed.
            :type genes_batch: str
            :param file_region_fasta: Path to fasta transcriptome sequence output file.
            :type file_region_fasta: string
            :return: Mapping of probes to corresponding genes with additional information about each probe, i.e.
                position (chromosome, start, end, strand), gene_id, transcript_id, exon_id
            :rtype: dict
            """

            total_probes = 0
            gene_probes = {}  # assiciate each id to the additional features
            gene_sequences = {}  # associate each sequence to its id
            id = 1

            # parse the exon fasta sequence file
            for exon in SeqIO.parse(file_region_fasta, "fasta"):
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

                            if probe_sequence in gene_sequences:
                                probe_id = gene_sequences[probe_sequence]
                                gene_probes[probe_id]["transcript_id"].append(
                                    transcript_id
                                )
                                gene_probes[probe_id]["exon_id"].append(exon_id)
                                gene_probes[probe_id]["start"].append(probe_start)
                                gene_probes[probe_id]["end"].append(probe_end)

                            else:
                                probe_id = f"{gene_id}_{id}"
                                gene_sequences[probe_sequence] = probe_id
                                id += 1
                                gene_probes[probe_id] = {
                                    "probe_sequence": probe_sequence,
                                    "transcript_id": [transcript_id],
                                    "exon_id": [exon_id],
                                    "chromosome": chrom,
                                    "start": [probe_start],
                                    "end": [probe_end],
                                    "strand": strand,
                                    "length": probe_length,
                                }
            # filter probes based on user-defined filters
            gene_oligo_DB = {}
            gene_oligo_DB[gene] = gene_probes

            return gene_oligo_DB

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

        probes = {}
        results = Parallel(n_jobs=n_jobs)(
            delayed(_get_probes)(gene, region_annotation) for gene in genes
        )

        for result in results:
            probes.update(result)

        return probes
