import copy
import os
import shutil
import warnings
from pathlib import Path

import IO._data_parser as data_parser
import IO._ftp_loader as ftp_loader
import pyfaidx
from oligo_transcript_generation._gene_transcript import GeneTranscript
from oligo_transcript_generation._oligos import Oligos


class CustomDB:
    """This class generates all possible guides that can be designed for a given list of genes,
    based on the transcriptome annotation or the gene CDS or the whole genome and the reference fasta file.
    The gtf and a fasta files are passed in input
    """

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        filters,
        species=None,
        genome_assembly=None,
        annotation_release=None,
        annotation_source=None,
        file_annotation=None,
        file_sequence=None,
    ):
        """Sets species, genome_assembly, annotation_release to 'unknown' if thay are not given in input
            Saves the path of the user defined annoation and fasta file and initializes the

        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int
        :param filters: list of filters classes already initialized
        :type filters: list of classes
        :param species: species of the fasta and gtf files, defaults to None
        :type species: str, optional
        :param genome_assembly: genome_assembly of the fasta and gtf files, defaults to None
        :type genome_assembly: str, optional
        :param annotation_release: annotation_release of the fasta and gtf files, defaults to None
        :type annotation_release: str, optional
        :param annotation_source: source of the fasta and gtf files, defaults to None
        :type annotation_source: str, optional
        :param file_annotation: path to the gtf annotation file, defaults to None
        :type file_annotation: str, optional
        :param file_sequence: path to the fasta file, defaults to None
        :type file_sequence: str, optional
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

        if annotation_source is None:
            annotation_source = "Custom"
            warnings.warn(f"Annotation source not specified.")

        # check the files format
        if file_annotation == None:
            raise ValueError("Annotation File not defined!")

        if file_sequence == None:
            raise ValueError("Sequence File not defined!")

        if not data_parser.check_gtf_format(file_annotation):
            raise ValueError("Annotation File has incorrect format!")

        if not data_parser.check_fasta_format(file_sequence):
            raise ValueError("Sequence File has incorrect format!")

        self.species = species
        self.genome_assembly = genome_assembly
        self.annotation_release = annotation_release
        self.annotation_source = annotation_source
        self.probe_length_max = probe_length_max
        self.probe_length_min = probe_length_min

        self.file_reference_DB = None
        self.file_oligos_DB_tsv = None
        self.file_oligos_DB_gtf = None
        self.oligos_DB = None

        self.file_sequence = file_sequence
        # create index file
        pyfaidx.Fasta(self.file_sequence)
        self.gene_transcript = GeneTranscript(self.file_sequence, file_annotation)
        self.oligos = Oligos(
            self.probe_length_min, self.probe_length_max, self.file_sequence, filters
        )

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

    def read_oligos_DB_gtf(self, file_oligos_DB_gtf):
        """Create the oligo db dictionary from a gtf file.

        :param file_oligos_DB_gtf: Path to the file.
        :type file_oligos_DB_gtf: str
        :raises ValueError: When the file given has the wrong format.
        :raises ValueError: When the path given does not exist.
        """

        if os.path.exists(file_oligos_DB_gtf):
            if data_parser.check_gtf_format(file_oligos_DB_gtf):
                self.file_oligos_DB_gtf = file_oligos_DB_gtf
            else:
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        oligos_df = data_parser.read_gtf(self.file_oligos_DB_gtf)
        self.oligos_DB = {}
        # compute the additional columns
        columns_fixed = [
            "transcript_id",
            "exon_id",
            "chromosome",
            "start",
            "end",
            "strand",
            "length",
            "probe_id",
        ]
        columns_df = list(oligos_df.columns)
        # remove the columns we don't need
        columns_df.remove("seqname")
        columns_df.remove("source")
        columns_df.remove("feature")
        columns_df.remove("score")
        columns_df.remove("frame")
        columns_df.remove("gene_id")
        columns_df.remove("probe_sequence")
        additional_columns = [
            column for column in columns_df if (column not in columns_fixed)
        ]
        # crated the oligos db
        current_gene = ""
        current_probe = ""
        for index in range(len(oligos_df)):
            row = oligos_df.iloc[index]
            if row["gene_id"] != current_gene:
                current_gene = row["gene_id"]
                self.oligos_DB[current_gene] = {}
            if (
                row["probe_sequence"] != current_probe
            ):  # easier to compare that the whole seqeunce
                current_probe = row["probe_sequence"]
                self.oligos_DB[current_gene][current_probe] = {}
                # add all the values
                self.oligos_DB[current_gene][current_probe]["transcript_id"] = [
                    row["transcript_id"]
                ]
                self.oligos_DB[current_gene][current_probe]["exon_id"] = [
                    row["exon_id"]
                ]
                self.oligos_DB[current_gene][current_probe]["chromosome"] = row[
                    "chromosome"
                ]
                self.oligos_DB[current_gene][current_probe]["start"] = [
                    int(row["start"])
                ]
                self.oligos_DB[current_gene][current_probe]["end"] = [int(row["end"])]
                self.oligos_DB[current_gene][current_probe]["strand"] = row["strand"]
                self.oligos_DB[current_gene][current_probe]["length"] = int(
                    row["length"]
                )
                self.oligos_DB[current_gene][current_probe]["probe_id"] = row["seqname"]
                for column in additional_columns:
                    self.oligos_DB[current_gene][current_probe][column] = float(
                        row[column]
                    )
            else:
                # append the values saved as a list
                self.oligos_DB[current_gene][current_probe]["transcript_id"].append(
                    row["transcript_id"]
                )
                self.oligos_DB[current_gene][current_probe]["exon_id"].append(
                    row["exon_id"]
                )
                self.oligos_DB[current_gene][current_probe]["start"].append(
                    row["start"]
                )
                self.oligos_DB[current_gene][current_probe]["end"].append(row["end"])

    def read_oligos_DB_tsv(self, file_oligos_DB_tsv):
        """Reads a previously generated oligos DB and saves it in the <self.oligos_DB> attribute as a dictionary.
        The order of columns is : probe_id, probe_sequence, gene_id,  'transcript_id', 'exon_id', 'chromosome', 'start', 'end', 'strand', all the additional info computed by the filtering class.

        :param file_oligos_DB_tsv: path of the oligos_DB file
        :type file_oligos_DB_tsv: str
        """

        def parse_line_tsv(line, current_gene, add_features):
            """Parses the lines of the tsv file of the oligos db and puts teh data in the dictionary.

            :param line: current line of the tsv file
            :type line: str
            :param current_gene: gene until wich we have created the dictionary
            :type current_gene: str
            :return: current gene we are after this iteration
            :rtype: str
            """
            line = line.split("\t")
            line[-1] = line[-1][0:-1]
            if line[2] != current_gene:
                current_gene = line[2]
                self.oligos_DB[current_gene] = {}
            sequence = line[1]
            # what if we have duplicated sequences?
            self.oligos_DB[current_gene][sequence] = {}
            self.oligos_DB[current_gene][sequence]["transcript_id"] = line[3].split(";")
            self.oligos_DB[current_gene][sequence]["exon_id"] = line[4].split(";")
            self.oligos_DB[current_gene][sequence]["chromosome"] = line[5]
            self.oligos_DB[current_gene][sequence]["start"] = list(
                map(int, line[6].split(";"))
            )
            self.oligos_DB[current_gene][sequence]["end"] = list(
                map(int, line[7].split(";"))
            )
            self.oligos_DB[current_gene][sequence]["strand"] = line[8]
            self.oligos_DB[current_gene][sequence]["length"] = int(line[9])
            self.oligos_DB[current_gene][sequence]["probe_id"] = line[0]
            # retrive the remaining features if they were computed
            if add_features:
                for i, column in enumerate(columns[10:]):
                    self.oligos_DB[current_gene][sequence][column] = float(line[i + 10])
            return current_gene

        if os.path.exists(file_oligos_DB_tsv):
            if data_parser.check_tsv_format(file_oligos_DB_tsv):
                self.file_oligos_DB_tsv = file_oligos_DB_tsv
            else:
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        self.file_oligos_DB_tsv = file_oligos_DB_tsv

        self.oligos_DB = {}
        handle_probe = open(self.file_oligos_DB_tsv, "r")
        # read the header
        line = handle_probe.readline()
        columns = line.split("\t")
        columns[-1] = columns[-1][0:-1]  # delete \n in the last word
        add_features = len(columns) > 10
        # read the rest of the file
        current_gene = ""
        for line in handle_probe:
            current_gene = parse_line_tsv(line, current_gene, add_features)

        handle_probe.close()

    def write_oligos_DB_gtf(self):
        """Writes the data structure self.oligos_DB in a gtf file in the <self.file_oligos_DB_gtf> path.
        The additional features are written in the 9th column.
        """
        with open(self.file_oligos_DB_gtf, "w") as handle_gtf:
            # write some comment?
            source = "oligo-designer-toolsuite"
            feature = "Oligonucleotide"
            score = "."
            frame = "."
            for gene_id, probe in self.oligos_DB.items():
                for probe_sequence, probe_attributes in probe.items():
                    for i in range(len(probe_attributes["start"])):
                        # write the annotation file
                        probe_id = probe_attributes["probe_id"]
                        start = str(probe_attributes["start"][i])
                        end = str(probe_attributes["end"][i])
                        strand = probe_attributes["strand"]
                        output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                            probe_id, source, feature, start, end, score, strand, frame
                        )
                        # add all the oher features
                        output += f'gene_id "{gene_id}"; '
                        output += f'probe_sequence "{probe_sequence}"; '
                        probe_attributes_copy = copy.deepcopy(probe_attributes)
                        # delete already written features
                        del probe_attributes_copy["start"]
                        del probe_attributes_copy["end"]
                        del probe_attributes_copy["strand"]
                        del probe_attributes_copy["probe_id"]
                        for key, value in probe_attributes_copy.items():
                            if type(value) == list:
                                output += f'{key} "{value[i]}"; '
                            else:
                                output += f'{key} "{value}"; '
                        output += "\n"
                        handle_gtf.write(output)

    def write_oligos_DB_tsv(self):
        """Writes the data structure self.oligos_DB in a tsv file in the <self.file_oligos_DB_tsv> path.
        The order of columns is : gene_id, probe_sequence, 'transcript_id', 'exon_id', 'chromosome', 'start', 'end', 'strand', 'length', all the additional info computed by the filtering class.
        """

        with open(self.file_oligos_DB_tsv, "w") as handle_probe:
            columns = [
                "probe_id",
                "probe_sequence",
                "gene_id",
                "transcript_id",
                "exon_id",
                "chromosome",
                "start",
                "end",
                "strand",
                "length",
            ]  # keep fixed the structure for these coulums
            # find all the other names of the columns (depend on the filters applied)
            genes = list(self.oligos_DB.keys())
            i = 0
            while self.oligos_DB[genes[i]] == {}:
                i += 1
            tmp = self.oligos_DB[genes[i]]
            tmp = list(list(tmp.values())[0].keys())
            additional_columns = [column for column in tmp if (column not in columns)]
            columns.extend(additional_columns)
            handle_probe.write("\t".join(columns) + "\n")

            for gene_id, probe in self.oligos_DB.items():
                for probe_sequence, probe_attributes in probe.items():
                    # write the basic information information we compute for each probe
                    output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                        probe_attributes["probe_id"],
                        probe_sequence,
                        gene_id,
                        ";".join(probe_attributes["transcript_id"]),
                        ";".join(probe_attributes["exon_id"]),
                        probe_attributes["chromosome"],
                        ";".join(str(s) for s in probe_attributes["start"]),
                        ";".join(str(e) for e in probe_attributes["end"]),
                        probe_attributes["strand"],
                        probe_attributes["length"],
                    )
                    # if we computed additional features we write also those
                    if len(additional_columns) > 0:
                        for column in additional_columns:
                            output += "\t{}".format(probe_attributes[column])
                        # \n at the end f the string
                    output += "\n"
                    handle_probe.write(output)

    def create_reference_DB(
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

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_reference_DB = os.path.join(
            dir_output,
            f"reference_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_genome_{genome}_gene_transcript_{gene_transcript}.fna",
        )

        if block_size is None:  # when not specified define it automatically
            block_size = self.probe_length_max + 5

        files_fasta = get_files_fasta(genome, gene_transcript, gene_CDS)
        data_parser.merge_fasta(files_fasta, self.file_reference_DB)

    def __generate_gene_CDS():
        """'Creates a fasta file containing the whole transcriptome."""
        # should be implemented in an external class
        raise NotImplementedError

    def create_oligos_DB(
        self,
        genes=None,
        region="gene_transcript",
        number_batchs=1,
        dir_output="./output/annotation",
        write_tsv=True,
        write_gtf=True,
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
                warnings.warn("Genome not implemented yet")
                # TODO
            elif region == "gene_transcript":
                file_region_annotation = self.gene_transcript.generate_for_oligos(
                    self.probe_length_max - 1, dir_output, genes
                )
            elif region == "gene_CDS":
                file_region_annotation = None
                warnings.warn("Gene CDS not implemented yet")
                # TODO
            else:
                raise ValueError(
                    f"The given region does not exists. You selected {region}"
                )
            return file_region_annotation

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        self.file_oligos_DB_tsv = os.path.join(
            dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.tsv",
        )
        self.file_oligos_DB_gtf = os.path.join(
            dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.gtf",
        )
        file_region_annotation = create_target_region(region, genes)
        if genes is None:
            genes = self.gene_transcript.get_genes_from_annotation()
        self.oligos_DB = self.oligos.generate(
            file_region_annotation, genes, number_batchs, dir_output
        )
        if write_tsv:
            self.write_oligos_DB_tsv()
        if write_gtf:
            self.write_oligos_DB_gtf()

        # clean folder
        os.remove(file_region_annotation)


class NcbiDB(CustomDB):
    """Class to create reference and oligos DB using gtf and fasta files taken from the NCBI server"""

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        filters,
        species=None,
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

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        genome_assembly = "GRCh38"
        annotation_source = "NCBI"

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        ftp = ftp_loader.FTPLoaderNCBI(dir_output, species, annotation_release)
        file_annotation = ftp.download_files("gtf")
        file_sequence = ftp.download_files("fasta")

        super().__init__(
            probe_length_min,
            probe_length_max,
            filters,
            species,
            genome_assembly,
            annotation_release,
            annotation_source,
            file_annotation,
            file_sequence,
        )


class EnsemblDB(CustomDB):
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
            genome_assembly = "GRCh38"
            warnings.warn(
                f"No genome assembly defined. Using default assembly {genome_assembly}!"
            )

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        annotation_source = "Ensembl"

        Path(dir_output).mkdir(parents=True, exist_ok=True)
        ftp = ftp_loader.FtpLoaderEnsembl(
            dir_output, species, genome_assembly, annotation_release
        )
        file_annotation = ftp.download_files("gtf")
        file_sequence = ftp.download_files("fasta")

        super().__init__(
            probe_length_min,
            probe_length_max,
            filters,
            species,
            genome_assembly,
            annotation_release,
            annotation_source,
            file_annotation,
            file_sequence,
        )
