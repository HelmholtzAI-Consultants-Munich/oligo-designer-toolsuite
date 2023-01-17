import os

import pandas as pd
from Bio import SeqIO
from joblib import Parallel, delayed

from ..utils._data_parser import get_sequence_from_annotation


class OligosGenerator:
    """Creates the oligo sequences belonging to the given genes.

    :param oligo_length_min: minimum length of the oligos created
    :type oligo_length_min: int
    :param oligo_length_max: maximum length of the oligos created
    :type oligo_length_max: int
    :param file_sequence: path to the fasta file of the whole genome
    :type file_sequence: str
    :param n_jobs: nr of cores used to compute the oligos
    :type n_jobs: int
    """

    def __init__(
        self,
        oligo_length_min: int,
        oligo_length_max: int,
        file_sequence: str,
        n_jobs: int,
    ):
        """Initialize the class."""
        self.oligo_length_min = oligo_length_min
        self.oligo_length_max = oligo_length_max
        self.file_sequence = file_sequence
        self.n_jobs = n_jobs

    def generate(self, file_region_annotation: str, genes: list[str], dir_output: str):
        """Get the fasta sequence of all possible oligos with user-defined length for all input genes and store the in a dictionary.
        Generated oligos are filtered by the filters give in input and the features computed for the  filters are added to the dictionary.

        :param file_region_annotation: path to the gtf annotaiton file of the region
        :type file_region_annotation: str
        :param genes: Genes for which the oligos are computed
        :type genes: list of str
        :param dir_output: directory where temporary diles are written, defaults to './output/annotation'
        :type dir_output: str, optional
        :return: the oligos DB containing all the oligos created
        :rtype: dict
        """

        self.dir_output = dir_output
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
        # compute the oligos
        oligos = {}
        results = Parallel(n_jobs=self.n_jobs)(
            delayed(self._get_oligos)(gene, region_annotation) for gene in genes
        )
        for result in results:
            oligos.update(result)

        return oligos

    def _get_oligos(self, gene, region_annotation):
        """Get the fasta sequence of all possible oligos in the current gene.

        :param gene: gene for which oligos should be designed.
        :type gene: str
        :param region_annotation: gtf annotation of the region
        :type region_annotation: pandas.DataFrame
        :return: the oligos DB containing all the oligos created for the gene
        :rtype: dict
        """

        file_region_bed = os.path.join(
            self.dir_output, "region_gene{}.bed".format(gene)
        )
        file_region_fasta = os.path.join(
            self.dir_output, "region_gene{}.fna".format(gene)
        )

        self._get_region_fasta(
            gene, region_annotation, file_region_bed, file_region_fasta
        )
        gene_oligos_batch = self._get_oligos_info(gene, file_region_fasta)

        os.remove(file_region_bed)
        os.remove(file_region_fasta)

        return gene_oligos_batch

    def _get_region_fasta(
        self, gene, region_annotation, file_region_bed, file_region_fasta
    ):
        """Extract transcripts for the current batch and write transcript regions to bed file.
        Get sequence for annotated transcript regions (bed file) from genome sequence (fasta file) and write transcriptome sequences to fasta file.

        :param genes_batch: Genes for which oligos should be designed.
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
        get_sequence_from_annotation(
            file_region_bed,
            self.file_sequence,
            file_region_fasta,
            split=True,
        )

    def _parse_header(self, header):
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

    def _get_oligos_info(self, gene, file_region_fasta):
        """Merge all oligos with identical sequence that come from the same gene into one fasta entry.
        Filter all oligos based on user-defined filters and collect the additional information about each oligo.

        :param gene: Gene for which oligos should be designed.
        :type genes_batch: str
        :param file_region_fasta: Path to fasta transcriptome sequence output file.
        :type file_region_fasta: string
        :return: Mapping of oligos to corresponding genes with additional information about each oligo, i.e.
            position (chromosome, start, end, strand), gene_id, transcript_id, exon_id
        :rtype: dict
        """

        total_oligos = 0
        gene_oligos = {}  # assiciate each id to the additional features
        gene_sequences = {}  # associate each sequence to its id
        id = 1

        # parse the exon fasta sequence file
        for exon in SeqIO.parse(file_region_fasta, "fasta"):
            sequence = exon.seq

            for oligo_length in range(self.oligo_length_min, self.oligo_length_max + 1):
                if len(sequence) > oligo_length:
                    number_oligos = len(sequence) - (oligo_length - 1)
                    oligos_sequence = [
                        sequence[i : i + oligo_length] for i in range(number_oligos)
                    ]

                    for i in range(number_oligos):
                        total_oligos += 1
                        oligo_sequence = oligos_sequence[i]

                        (
                            gene_id,
                            transcript_id,
                            exon_id,
                            chrom,
                            start,
                            strand,
                        ) = self._parse_header(exon.id)
                        oligo_start = start + i
                        oligo_end = start + i + oligo_length

                        if oligo_sequence in gene_sequences:
                            oligo_id = gene_sequences[oligo_sequence]
                            gene_oligos[oligo_id]["transcript_id"].append(transcript_id)
                            gene_oligos[oligo_id]["exon_id"].append(exon_id)
                            gene_oligos[oligo_id]["start"].append(oligo_start)
                            gene_oligos[oligo_id]["end"].append(oligo_end)

                        else:
                            oligo_id = f"{gene_id}_{id}"
                            gene_sequences[oligo_sequence] = oligo_id
                            id += 1
                            gene_oligos[oligo_id] = {
                                "sequence": oligo_sequence,
                                "transcript_id": [transcript_id],
                                "exon_id": [exon_id],
                                "chromosome": chrom,
                                "start": [oligo_start],
                                "end": [oligo_end],
                                "strand": strand,
                                "length": oligo_length,
                            }
        # filter oligos based on user-defined filters
        gene_oligo_DB = {}
        gene_oligo_DB[gene] = gene_oligos

        return gene_oligo_DB
