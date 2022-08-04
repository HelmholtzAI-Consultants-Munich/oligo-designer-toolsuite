import copy
import os
import shutil
import warnings
from pathlib import Path

import pyfaidx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import oligo_designer_toolsuite.IO._data_parser as data_parser
import oligo_designer_toolsuite.IO._ftp_loader as ftp_loader
from oligo_designer_toolsuite.oligo_transcript_generation._gene_transcript import (
    GeneTranscript,
)
from oligo_designer_toolsuite.oligo_transcript_generation._oligos import Oligos


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
        dir_output="output",
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

        self.dir_output = os.path.join(dir_output, "annotation")  # choose a different
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)
        self.species = species
        self.genome_assembly = genome_assembly
        self.annotation_release = annotation_release
        self.annotation_source = annotation_source
        self.probe_length_max = probe_length_max
        self.probe_length_min = probe_length_min

        self.file_reference_DB = None
        self.file_oligos_DB_tsv = None
        self.file_oligos_DB_gtf = None
        self.file_oligos_DB_fasta = None
        self.oligos_DB = None

        self.file_sequence = file_sequence
        # create index file
        pyfaidx.Fasta(self.file_sequence)
        self.gene_transcript = GeneTranscript(self.file_sequence, file_annotation)
        self.oligos = Oligos(
            self.probe_length_min, self.probe_length_max, self.file_sequence, filters
        )

    # move read ad write to data parser and leave 2 methods read and write wich take as input teh type and the folder where to write (optional)
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

    def read_oligos_DB_gtf(self, file_oligos_DB_gtf, file_oligos_DB_fasta):
        """Create the oligo db dictionary from a gtf file.

        :param file_oligos_DB_gtf: Path to the file.
        :type file_oligos_DB_gtf: str
        :param file_oligos_DB_fasta: Path to the file.
        :type file_oligos_DB_fasta: str
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

        if os.path.exists(file_oligos_DB_fasta):
            if data_parser.check_fasta_format(file_oligos_DB_fasta):
                self.file_oligos_DB_fasta = file_oligos_DB_fasta
            else:
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        oligos_df = data_parser.read_gtf(self.file_oligos_DB_gtf)
        oligos_fasta = SeqIO.parse(self.file_oligos_DB_fasta, "fasta")
        self.oligos_DB = {}
        # compute the additional columns
        columns_fixed = [
            "probe_sequence",
            "transcript_id",
            "exon_id",
            "chromosome",
            "start",
            "end",
            "strand",
            "length",
        ]
        columns_df = list(oligos_df.columns)
        # remove the columns we don't need
        columns_df.remove("seqname")
        columns_df.remove("source")
        columns_df.remove("feature")
        columns_df.remove("score")
        columns_df.remove("frame")
        columns_df.remove("gene_id")
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
            if row["seqname"] != current_probe:
                current_probe = row["seqname"]
                self.oligos_DB[current_gene][current_probe] = {}
                oligo_fasta = next(oligos_fasta)
                probe_sequence = oligo_fasta.seq
                assert (
                    oligo_fasta.id == current_probe
                )  # check that the fasta and gtf file are in sync
                self.oligos_DB[current_gene][current_probe] = {}
                # add all the values
                self.oligos_DB[current_gene][current_probe][
                    "probe_sequence"
                ] = probe_sequence
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
            probe_id = line[0]
            # what if we have duplicated sequences?
            self.oligos_DB[current_gene][probe_id] = {}
            self.oligos_DB[current_gene][probe_id]["probe_sequence"] = Seq(line[1])
            self.oligos_DB[current_gene][probe_id]["transcript_id"] = line[3].split(";")
            self.oligos_DB[current_gene][probe_id]["exon_id"] = line[4].split(";")
            self.oligos_DB[current_gene][probe_id]["chromosome"] = line[5]
            self.oligos_DB[current_gene][probe_id]["start"] = list(
                map(int, line[6].split(";"))
            )
            self.oligos_DB[current_gene][probe_id]["end"] = list(
                map(int, line[7].split(";"))
            )
            self.oligos_DB[current_gene][probe_id]["strand"] = line[8]
            self.oligos_DB[current_gene][probe_id]["length"] = int(line[9])
            # retrive the remaining features if they were computed
            if add_features:
                for i, column in enumerate(columns[10:]):
                    self.oligos_DB[current_gene][probe_id][column] = float(line[i + 10])
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
        The additional features are written in the 9th column and the sequence of the probes is written on a separate fasta fila
        with heading the probe_id.
        """
        with open(self.file_oligos_DB_gtf, "w") as handle_gtf:
            with open(self.file_oligos_DB_fasta, "w") as handle_fasta:
                source = "oligo-designer-toolsuite"
                feature = "Oligonucleotide"
                score = "."
                frame = "."
                output_fasta = []
                for gene_id, probe in self.oligos_DB.items():
                    for probe_id, probe_attributes in probe.items():
                        output_fasta.append(
                            SeqRecord(
                                probe_attributes["probe_sequence"], probe_id, "", ""
                            )
                        )  # write the sequence in the fasta file
                        for i in range(len(probe_attributes["start"])):
                            # write the annotation file
                            start = str(probe_attributes["start"][i])
                            end = str(probe_attributes["end"][i])
                            strand = probe_attributes["strand"]
                            output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                                probe_id,
                                source,
                                feature,
                                start,
                                end,
                                score,
                                strand,
                                frame,
                            )
                            # add all the oher features
                            output += f'gene_id "{gene_id}"; '
                            probe_attributes_copy = copy.deepcopy(probe_attributes)
                            # delete already written features
                            del probe_attributes_copy["start"]
                            del probe_attributes_copy["end"]
                            del probe_attributes_copy["strand"]
                            del probe_attributes_copy["probe_sequence"]
                            for key, value in probe_attributes_copy.items():
                                if type(value) == list:
                                    output += f'{key} "{value[i]}"; '
                                else:
                                    output += f'{key} "{value}"; '
                            output += "\n"
                            handle_gtf.write(output)
                SeqIO.write(output_fasta, handle_fasta, "fasta")

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
                for probe_id, probe_attributes in probe.items():
                    # write the basic information information we compute for each probe
                    output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                        probe_id,
                        probe_attributes["probe_sequence"],
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
        region="gene_transcript",
        block_size=None,
    ):
        """Creates a fasta file for each of the region selected (genome, gene_transcript, gene_CDS) which will be used for alignements, default is "gene_transcript".
        If not specified the exon juctions size is set to <probe_length_max> + 5.

        :param region: the region to use for the reference DB. Possible values are "genome", "gene_transcript", "gene_CDS"
        :type region: str
        :param block_size: size of the exon junctions, defaults to None
        :type block_size: int, optional
        """

        def get_files_fasta(region):
            """generates the fasta files that will compose the reference_DB

            :param region: the region to use for the reference DB. Possible values are "genome", "gene_transcript", "gene_CDS".
            :type region: str
            :return: list of the fasta files
            :rtype: list of str
            """
            file_fasta = None
            if region == "genome":
                file_fasta = os.path.join(self.dir_output, "genome.fna")
                shutil.copyfile(self.file_sequence, file_fasta)
            elif region == "gene_transcript":
                (
                    file_gene_transcript_annotation,
                    file_gene_transcript_fasta,
                ) = self.gene_transcript.generate_for_reference(
                    block_size, self.dir_output
                )  # call the outer class to generate the gene trascript
                os.remove(
                    file_gene_transcript_annotation
                )  # not required anymore in this case
                file_fasta = file_gene_transcript_fasta
            elif region == "gene_CDS":
                # generate gene cds
                warnings.warn("Gene CDS not implemented yet")
            else:
                raise ValueError(
                    f"The given region does not exists. You selected {region}"
                )
            return file_fasta

        self.file_reference_DB = os.path.join(
            self.dir_output,
            f"reference_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.fna",
        )

        if block_size is None:  # when not specified define it automatically
            block_size = self.probe_length_max + 5

        self.file_reference_DB = get_files_fasta(region)

    def __generate_gene_CDS():
        """'Creates a fasta file containing the whole transcriptome."""
        # should be implemented in an external class
        raise NotImplementedError

    def create_oligos_DB(  # does not automatically write, we haveto call it speately
        self,
        genes=None,
        region="gene_transcript",
        n_jobs=2,
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
                    self.probe_length_max - 1, self.dir_output, genes
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

        self.file_oligos_DB_tsv = os.path.join(
            self.dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.tsv",
        )
        self.file_oligos_DB_gtf = os.path.join(
            self.dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.gtf",
        )
        self.file_oligos_DB_fasta = os.path.join(
            self.dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.fasta",
        )
        file_region_annotation = create_target_region(region, genes)
        if genes is None:
            genes = self.gene_transcript.get_genes_from_annotation()
        self.oligos_DB = self.oligos.generate(
            file_region_annotation, genes, n_jobs, self.dir_output
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
        dir_output="output",
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
        self.dir_output = os.path.join(dir_output, "annotation")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ftp = ftp_loader.FTPLoaderNCBI(self.dir_output, species, annotation_release)
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
            dir_output,
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
        dir_output="output",
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
        :param dir_output: directory where the files are saved, defaults to 'output'
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
        self.dir_output = os.path.join(dir_output, "annotation")

        Path(self.dir_output).mkdir(parents=True, exist_ok=True)
        ftp = ftp_loader.FtpLoaderEnsembl(
            self.dir_output, species, genome_assembly, annotation_release
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
            dir_output,
        )
