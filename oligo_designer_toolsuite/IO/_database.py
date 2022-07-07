import os
import random
import shutil
import warnings
from pathlib import Path

import IO._data_parser as data_parser
import IO._ftp_loader as ftp_loader
import pandas as pd
import pyfaidx
from Bio import SeqIO
from oligo_transcript_generation._gene_transcript import GeneTranscript


class BaseDB:
    """This class generates all possible guides that can be designed for a given list of genes,
    based on the transcriptome annotation or of those genes and the reference fasta file.
    """

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        filters,
        species=None,
        genome_assembly=None,
        annotation_release=None,
    ):
        """Initializes the class.

        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int
        :param species: species of the files to dowload, defaults to None
        :type species: str, optional
        :param genome_assembly: genome_assembly of the files to dowload, defaults to None
        :type genome_assembly: str, optional
        :param annotation_release: annotation_release of the files to dowload, defaults to None
        :type annotation_release: str, optional
        :param filters: list of filters classes already initialized, defaults to None
        :type filters: list of classes, optional
        """

        self.species = species
        self.genome_assembly = genome_assembly
        self.annotation_release = annotation_release
        self.filters = filters
        self.probe_length_max = probe_length_max
        self.probe_length_min = probe_length_min

        self.annotation = None  # dataframe containing the gene annotation
        self.file_reference_DB = None
        self.file_oligos_DB = None
        self.oligos_DB = None
        self.file_annotation = None
        self.file_sequence = None
        self.gene_transcript = None

    def read_reference_DB(self, file_reference_DB):
        """Saves the path of a previously generated reference DB in the <self.file_reference_DB> attribute.

        :param file_reference_DB: path of the reference_DB file
        :type file_reference_DB: str
        """
        if os.path.exists(file_reference_DB):
            if data_parser.check_fasta_format(file_reference_DB):
                self.file_reference_DB = file_reference_DB
            else:
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

    def read_oligos_DB(self, file_oligos_DB):
        """Reads a previously generated oligos DB and saves it in the <self.oligos_DB> attribute as a dictionary.
        The order of columns is : gene_id, probe_sequence, 'transcript_id', 'exon_id', 'chromosome', 'start', 'end', 'strand', all the additional info computed by the filtering class.

        :param file_oligos_DB: path of the oligos_DB file
        :type file_oligos_DB: str
        """
        if os.path.exists(file_oligos_DB):
            if data_parser.check_tsv_format(file_oligos_DB):
                self.file_oligos_DB = file_oligos_DB
            else:
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        self.file_oligos_DB = file_oligos_DB

        self.oligos_DB = {}
        handle_probe = open(self.file_oligos_DB, "r")
        # read the header
        line = handle_probe.readline()
        columns = line.split("\t")
        columns[-1] = columns[-1][0:-1]  # delete \n in the last word
        add_features = len(columns) > 8
        # read the first line
        line = handle_probe.readline()
        line = line.split("\t")
        line[-1] = line[-1][0:-1]  # delete \n of the last word
        current_gene = line[0]
        sequence = line[1]
        self.oligos_DB[current_gene] = {}
        self.oligos_DB[current_gene][sequence] = {}
        self.oligos_DB[current_gene][sequence]["transcript_id"] = line[2].split(";")
        self.oligos_DB[current_gene][sequence]["exon_id"] = line[3].split(";")
        self.oligos_DB[current_gene][sequence]["chromosome"] = line[4]
        self.oligos_DB[current_gene][sequence]["start"] = line[5].split(";")
        self.oligos_DB[current_gene][sequence]["end"] = line[6].split(";")
        self.oligos_DB[current_gene][sequence]["strand"] = line[7]
        # retrive the remaining features if they were computed
        if add_features:
            for i, column in enumerate(columns[8:]):
                self.oligos_DB[current_gene][sequence][column] = line[i + 8]
        # read the rest of the file
        for line in handle_probe:
            line = line.split("\t")
            line[-1] = line[-1][0:-1]
            if line[0] != current_gene:
                current_gene = line[0]
                self.oligos_DB[current_gene] = {}
            sequence = line[1]
            # what if we have duplicated sequences?
            self.oligos_DB[current_gene][sequence] = {}
            self.oligos_DB[current_gene][sequence]["transcript_id"] = line[2].split(";")
            self.oligos_DB[current_gene][sequence]["exon_id"] = line[3].split(";")
            self.oligos_DB[current_gene][sequence]["chromosome"] = line[4]
            self.oligos_DB[current_gene][sequence]["start"] = line[5].split(";")
            self.oligos_DB[current_gene][sequence]["end"] = line[6].split(";")
            self.oligos_DB[current_gene][sequence]["strand"] = line[7]
            # retrive the remaining features if they were computed
            if add_features:
                for i, column in enumerate(columns[8:]):
                    self.oligos_DB[current_gene][sequence][column] = line[i + 8]

        handle_probe.close()

    def write_oligos_DB(self):
        """Writes the data structure self.oligos_DB in a tsv file in the <self.file_oligos_DB> path.
        The order of columns is : gene_id, probe_sequence, 'transcript_id', 'exon_id', 'chromosome', 'start', 'end', 'strand', all the additional info computed by the filtering class.
        """

        with open(self.file_oligos_DB, "w") as handle_probe:
            columns = ["gene_id", "probe_sequence"]
            # find all the other names of the columns (depend on the filters applied)
            tmp = list(self.oligos_DB.values())[0]
            tmp = list(list(tmp.values())[0].keys())  # what if the gene has no oligos?
            columns.extend(tmp)
            handle_probe.write("\t".join(columns) + "\n")

            for gene_id, probe in self.oligos_DB.items():
                for probe_sequence, probe_attributes in probe.items():
                    # write the basic information information we compute for each probe
                    output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                        gene_id,
                        probe_sequence,
                        ";".join(probe_attributes["transcript_id"]),
                        ";".join(probe_attributes["exon_id"]),
                        probe_attributes["chromosome"],
                        ";".join(str(s) for s in probe_attributes["start"]),
                        ";".join(str(e) for e in probe_attributes["end"]),
                        probe_attributes["strand"],
                    )
                    # if we computed additional features we write also those
                    if len(columns) > 8:
                        for column in columns[8:]:
                            output += "\t{}".format(probe_attributes[column])
                        # \n at the end f the string
                    output += "\n"
                    handle_probe.write(output)

    def _create_reference_DB(
        self,
        genome=False,
        gene_transcript=True,
        gene_CDS=False,
        dir_output="./output/annotation",
        block_size=None,
    ):
        """Creates a fasta file for each of the region selected (genome, gene_transcript, gene_CDS) which will be used for alignements.
        If not specified the exon juctions size is set to <probe_length_max> + 5.

        :param genome: create the reference file for the whole genome, defaults to False
        :type genome: bool, optional
        :param gene_transcript: create the reference file for the gene transcript, defaults to True
        :type gene_transcript: bool, optional
        :param gene_CDS: create the reference file for the coding region, defaults to False
        :type gene_CDS: bool, optional
        :param dir_output: folder where the fasta file will be saved, defaults to './output/annotation'
        :type dir_output: str, optional
        :param block_size: size of the exon junctions, defaults to None
        :type block_size: int, optional
        """

        def get_files_fasta(genome, gene_transcript, gene_CDS):
            """generates the fasta files that will compose the reference_DB

            :return: list of the fasta files
            :rtype: list of str
            """
            files_fasta = []
            if genome:
                file_genome = os.path.join(dir_output, "genome.fna")
                shutil.copyfile(self.file_sequence, file_genome)
                files_fasta.append(file_genome)
            if gene_transcript:
                (
                    file_gene_transcript_annotation,
                    file_gene_transcript_fasta,
                ) = self.gene_transcript.generate_for_reference(
                    block_size, dir_output
                )  # call the outer class to generate the gene trascript
                os.remove(
                    file_gene_transcript_annotation
                )  # not required anymore in this case
                files_fasta.append(file_gene_transcript_fasta)
            if gene_CDS:
                # generate gene cds
                warnings.warn("Gene CDS not implemented yet")
            return files_fasta

        if block_size is None:  # when not specified define it automatically
            block_size = self.probe_length_max + 5

        files_fasta = get_files_fasta(genome, gene_transcript, gene_CDS)
        data_parser.merge_fasta(files_fasta, self.file_reference_DB)

    def __generate_gene_CDS():
        """'Creates a fasta file containing the whole transcriptome."""

    def _create_oligos_DB(
        self,
        genes=None,
        region="gene_transcript",
        number_batchs=1,
        dir_output="./output/annotation",
        write=True,
    ):
        """creates the DB containing all the oligo sequence extracted form the given <region> and belonging the the specified genes. If no genes are specified then
        will be used all the genes. The DB is a dictionary data structure and can be written in a tsv format by setting <write> = True.

        :param genes: genes for which compute the probes, defaults to None
        :type genes: list of str, optional
        :param region: region ofrm whihc generate the probes, it can be 'genome', 'gene_transcript', 'gene_CDS', defaults to 'gene_transcript'
        :type region: str, optional
        :param number_batchs: probes are computes in batches of genes, defaults to 1
        :type number_batchs: int, optional
        :param dir_output: folder where the file will be written, defaults to './output/annotation'
        :type dir_output: str, optional
        :param write: write the file, defaults to True
        :type write: bool, optional
        """

        def create_target_region(region, genes):
            """cretares the annotation file for the specified region."""

            if region == "genome":
                file_region_annotation = None
                # TODO
            elif region == "gene_transcript":
                file_region_annotation = self.gene_transcript.generate_for_oligos(
                    self.probe_length_max - 1, dir_output, genes
                )
            elif region == "gene_CDS":
                file_region_annotation = None
                # TODO
            else:
                raise ValueError(
                    f"The given region does not exists. You selected {region}"
                )
            return file_region_annotation

        def get_gene_list_from_annotation():
            """generates all the different genes in the annotation file of the specified region."""

            genes = self.annotation.loc[self.annotation["feature"] == "gene"]
            genes = list(genes["gene_id"].unique())
            random.shuffle(genes)
            return genes

        file_region_annotation = create_target_region(region, genes)
        if genes is None:
            genes = get_gene_list_from_annotation()
        self.oligos_DB = self.__generate_oligos(
            file_region_annotation, genes, number_batchs, dir_output
        )
        if write:
            self.write_oligos_DB()

        # clean folder
        os.remove(file_region_annotation)

    def __generate_oligos(
        self,
        file_region_annotation,
        genes,
        number_batchs,
        dir_output="./output/annotation",
    ):
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
            loaded_probes = 0

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

                            if (
                                probe_sequence in gene_probes[gene_id]
                            ):  # already checked if fulfills
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
                                # check if the probe fulfills the filters conditions
                                fulfills, additional_features = self.filter_oligo(
                                    probe_sequence
                                )
                                if fulfills:
                                    loaded_probes += 1
                                    gene_probes[gene_id][probe_sequence] = {
                                        "transcript_id": [transcript_id],
                                        "exon_id": [exon_id],
                                        "chromosome": chrom,
                                        "start": [probe_start],
                                        "end": [probe_end],
                                        "strand": strand,
                                    }
                                    gene_probes[gene_id][probe_sequence].update(
                                        additional_features
                                    )

            print(f"the total number of probes found: {loaded_probes}")
            return gene_probes

        # create index file
        pyfaidx.Fasta(self.file_sequence)

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

        # remove index file
        os.remove("{}.fai".format(self.file_sequence))

        return probes

    def filter_oligo(self, sequence):
        """Applies the used-defined filters and returns the result and the additional computed features

        :param sequence: _description_
        :type sequence: _type_
        :return: _description_
        :rtype: _type_
        """

        fulfills = True
        additional_features = {}
        for filter in self.filters:
            out, feat = filter.apply(sequence)
            if not out:  # stop at the first false we obtain
                return False, {}
            additional_features.update(feat)
        return fulfills, additional_features


class NcbiDB(BaseDB):
    """Class to create reference and oligos DB using gtf and fasta files taken from the NCBI server"""

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        filters,
        species=None,
        genome_assembly=None,
        annotation_release=None,
        dir_output="./output/annotation",
    ):
        """Sets species, genome_assembly, annotation_release to a predefined value if thay are not given in input
            Dowloads the fasta and annotation files from the NCBI server and stores them in the dir_output folder

        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int
        :param species: species of the files to dowload, defaults to None
        :type species: str, optional
        :param genome_assembly: genome_assembly of the files to dowload, defaults to None
        :type genome_assembly: str, optional
        :param annotation_release: annotation_release of the files to dowload, defaults to None
        :type annotation_release: str, optional
        :param dir_output: directory where the files are saved, defaults to './output/annotation'
        :type dir_output: str, optional
        :param filters: list of filters classes already initialized, defaults to None
        :type filters: list of classes, optional
        """
        if species is None:
            species = "human"
            warnings.warn(f"No species defined. Using default species {species}!")

        if genome_assembly is None:
            genome_assembly = "hg38"
            warnings.warn(
                f"No genome assembly defined. Using default assembly {genome_assembly}!"
            )

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        super().__init__(
            probe_length_min,
            probe_length_max,
            filters,
            species,
            genome_assembly,
            annotation_release,
        )

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        ftp = ftp_loader.FTPLoaderNCBI(
            dir_output, species, genome_assembly, annotation_release
        )
        self.file_annotation = ftp.download_files("gtf")
        self.file_sequence = ftp.download_files("fasta")
        self.gene_transcript = GeneTranscript(self.file_sequence, self.file_annotation)

    def create_reference_DB(
        self,
        genome=False,
        gene_transcript=True,
        gene_CDS=False,
        dir_output="./output/annotation",
        block_size=None,
    ):
        """Creates a fasta file for each of the region selected (genome, gene_transcript, gene_CDS) which will be used for alignements.

        :param genome: create the reference file for the whole genome, defaults to False
        :type genome: bool, optional
        :param gene_transcript: create the reference file for the gene transcript, defaults to True
        :type gene_transcript: bool, optional
        :param gene_CDS: create the reference file for the coding region, defaults to False
        :type gene_CDS: bool, optional
        :param dir_output: folder where the fasta file will be saved, defaults to './output/annotation'
        :type dir_output: str, optional
        :param block_size: size of the exon junctions, defaults to None
        :type block_size: int, optional
        """

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_reference_DB = os.path.join(
            dir_output,
            f"reference_DB_{self.species}_{self.genome_assembly}_NCBI_release_{self.annotation_release}_genome_{genome}_gene_transcript_{gene_transcript}",
        )
        self._create_reference_DB(
            genome, gene_transcript, gene_CDS, dir_output, block_size
        )

    def create_oligos_DB(
        self,
        genes=None,
        region="gene_transcript",
        number_batchs=1,
        dir_output="./output/annotation",
        write=True,
    ):
        """creates the DB containing all the oligo sequence extracted form the given <region> and belonging the the specified genes. If no genes are specified then
        will be used all the genes. The DB is a dictionary data structure and can be written in a tsv format by setting <write> = True.

        :param genes: genes for which compute the probes, defaults to None
        :type genes: list of str, optional
        :param region: region ofrm whihc generate the probes, it can be 'genome', 'gene_transcript', 'gene_CDS', defaults to 'gene_transcript'
        :type region: str, optional
        :param number_batchs: probes are computes in batches of genes, defaults to 1
        :type number_batchs: int, optional
        :param dir_output: folder where the file will be written, defaults to './output/annotation'
        :type dir_output: str, optional
        :param write: write the file, defaults to True
        :type write: bool, optional
        """

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_oligos_DB = os.path.join(
            dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_NCBI_release_{self.annotation_release}_{region}",
        )
        self._create_oligos_DB(genes, region, number_batchs, dir_output, write)


class EnsemblDB(BaseDB):
    """Class to create reference and oligos DB using gtf and fasta files taken from the Ensembl server"""

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        filters,
        species=None,
        genome_assembly=None,
        annotation_release=None,
        dir_output="./output/annotation",
    ):
        """Sets species, genome_assembly, annotation_release to a predefined value if thay are not given in input
            Dowloads the fasta and annotation files from the Ensemble server and stores them in the dir_output folder

        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int
        :param species: species of the files to dowload, defaults to None
        :type species: str, optional
        :param genome_assembly: genome_assembly of the files to dowload, defaults to None
        :type genome_assembly: str, optional
        :param annotation_release: annotation_release of the files to dowload, defaults to None
        :type annotation_release: str, optional
        :param dir_output: directory where the files are saved, defaults to './output/annotation'
        :type dir_output: str, optional
        :param filters: list of filters classes already initialized, defaults to None
        :type filters: list of classes, optional
        """
        if species is None:  # change to some standard values for Ensemble
            species = "human"
            warnings.warn(f"No species defined. Using default species {species}!")

        if genome_assembly is None:
            genome_assembly = "hg38"
            warnings.warn(
                f"No genome assembly defined. Using default assembly {genome_assembly}!"
            )

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        super().__init__(
            probe_length_min,
            probe_length_max,
            filters,
            species,
            genome_assembly,
            annotation_release,
        )

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        ftp = ftp_loader.FtpLoaderEnsembl(
            species, dir_output, genome_assembly, annotation_release
        )
        self.file_annotation = ftp.download_files("gtf")
        self.file_sequence = ftp.download_files("fasta")
        self.gene_transcript = GeneTranscript(self.file_sequence, self.file_annotation)

    def create_reference_DB(
        self,
        genome=False,
        gene_transcript=True,
        gene_CDS=False,
        dir_output="./output/annotation",
        block_size=None,
    ):
        """Creates a fasta file for each of the region selected (genome, gene_transcript, gene_CDS) which will be used for alignements.

        :param genome: create the reference file for the whole genome, defaults to False
        :type genome: bool, optional
        :param gene_transcript: create the reference file for the gene transcript, defaults to True
        :type gene_transcript: bool, optional
        :param gene_CDS: create the reference file for the coding region, defaults to False
        :type gene_CDS: bool, optional
        :param dir_output: folder where the fasta file will be saved, defaults to './output/annotation'
        :type dir_output: str, optional
        :param block_size: size of the exon junctions, defaults to None
        :type block_size: int, optional
        """

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_reference_DB = os.path.join(
            dir_output,
            f"reference_DB_{self.species}_{self.genome_assembly}_Ensembl_release_{self.annotation_release}_genome_{genome}_gene_transcript_{gene_transcript}",
        )
        self._create_reference_DB(
            genome, gene_transcript, gene_CDS, dir_output, block_size
        )

    def create_oligos_DB(
        self,
        genes=None,
        region="gene_transcript",
        number_batchs=1,
        dir_output="./output/annotation",
        write=True,
    ):
        """creates the DB containing all the oligo sequence extracted form the given <region> and belonging the the specified genes. If no genes are specified then
        will be used all the genes. The DB is a dictionary data structure and can be written in a tsv format by setting <write> = True.

        :param genes: genes for which compute the probes, defaults to None
        :type genes: list of str, optional
        :param region: region ofrm whihc generate the probes, it can be 'genome', 'gene_transcript', 'gene_CDS', defaults to 'gene_transcript'
        :type region: str, optional
        :param number_batchs: probes are computes in batches of genes, defaults to 1
        :type number_batchs: int, optional
        :param dir_output: folder where the file will be written, defaults to './output/annotation'
        :type dir_output: str, optional
        :param write: write the file, defaults to True
        :type write: bool, optional
        """

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_oligos_DB = os.path.join(
            dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_Ensembl_release_{self.annotation_release}_{region}",
        )
        self._create_oligos_DB(genes, region, number_batchs, dir_output, write)


class CustomDB(BaseDB):
    """Class to create reference and oligos DB where gtf and fasta files are given by the user"""

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        filters,
        species=None,
        genome_assembly=None,
        annotation_release=None,
        file_annotation=None,
        file_sequence=None,
    ):
        """Sets species, genome_assembly, annotation_release to 'unknown' if thay are not given in input
            Saves the path of the user defined annoation and fasta file

        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int
        :param species: species of the files to dowload, defaults to None
        :type species: str, optional
        :param genome_assembly: genome_assembly of the files to dowload, defaults to None
        :type genome_assembly: str, optional
        :param annotation_release: annotation_release of the files to dowload, defaults to None
        :type annotation_release: str, optional
        :param file_annotation: path to the gtf annotation file, defaults to None
        :type file_annotation: str, optional
        :param file_sequence: path to the fasta file, defaults to None
        :type file_sequence: str, optional
        :param filters: list of filters classes already initialized, defaults to None
        :type filters: list of classes, optional
        """
        if species is None:
            species = "unknown"
            warnings.warn("Species not specified.")

        if genome_assembly is None:
            genome_assembly = "unknown"
            warnings.warn(f"Genome assembly not specified.")

        if annotation_release is None:
            annotation_release = "unknown"
            warnings.warn(f"Annotation release not specified.")

        super().__init__(
            probe_length_min,
            probe_length_max,
            filters,
            species,
            genome_assembly,
            annotation_release,
        )

        # check the filles format
        if file_annotation == None:
            raise ValueError("Annotation File not defined!")

        if file_sequence == None:
            raise ValueError("Sequence File not defined!")

        if not data_parser.check_gtf_format(file_annotation):
            raise ValueError("Annotation File has incorrect format!")

        if not data_parser.check_fasta_format(file_sequence):
            raise ValueError("Sequence File has incorrect format!")

        self.file_annotation = file_annotation
        self.file_sequence = file_sequence
        self.gene_transcript = GeneTranscript(self.file_sequence, self.file_annotation)

    def create_reference_DB(
        self,
        genome=False,
        gene_transcript=True,
        gene_CDS=False,
        dir_output="./output/annotation",
        block_size=None,
    ):
        """Creates a fasta file for each of the region selected (genome, gene_transcript, gene_CDS) which will be used for alignements.

        :param genome: create the reference file for the whole genome, defaults to False
        :type genome: bool, optional
        :param gene_transcript: create the reference file for the gene transcript, defaults to True
        :type gene_transcript: bool, optional
        :param gene_CDS: create the reference file for the coding region, defaults to False
        :type gene_CDS: bool, optional
        :param dir_output: folder where the fasta file will be saved, defaults to './output/annotation'
        :type dir_output: str, optional
        :param block_size: size of the exon junctions, defaults to None
        :type block_size: int, optional
        """

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_reference_DB = self.file_reference_DB = os.path.join(
            dir_output,
            f"reference_DB_{self.species}_{self.genome_assembly}_Custom_release_{self.annotation_release}_genome_{genome}_gene_transcript_{gene_transcript}",
        )
        self._create_reference_DB(
            genome, gene_transcript, gene_CDS, dir_output, block_size
        )

    def create_oligos_DB(
        self,
        genes=None,
        region="gene_transcript",
        number_batchs=1,
        dir_output="./output/annotation",
        write=True,
    ):
        """creates the DB containing all the oligo sequence extracted form the given <region> and belonging the the specified genes. If no genes are specified then
        will be used all the genes. The DB is a dictionary data structure and can be written in a tsv format by setting <write> = True.

        :param genes: genes for which compute the probes, defaults to None
        :type genes: list of str, optional
        :param region: region ofrm whihc generate the probes, it can be 'genome', 'gene_transcript', 'gene_CDS', defaults to 'gene_transcript'
        :type region: str, optional
        :param number_batchs: probes are computes in batches of genes, defaults to 1
        :type number_batchs: int, optional
        :param dir_output: folder where the file will be written, defaults to './output/annotation'
        :type dir_output: str, optional
        :param write: write the file, defaults to True
        :type write: bool, optional
        """

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_oligos_DB = os.path.join(
            dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_Custom_release_{self.annotation_release}_{region}",
        )
        self._create_oligos_DB(genes, region, number_batchs, dir_output, write)
