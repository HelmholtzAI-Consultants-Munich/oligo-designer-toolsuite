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
        self.file_oligos_DB = None
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
            f"reference_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_genome_{genome}_gene_transcript_{gene_transcript}",
        )

        if block_size is None:  # when not specified define it automatically
            block_size = self.probe_length_max + 5

        files_fasta = get_files_fasta(genome, gene_transcript, gene_CDS)
        data_parser.merge_fasta(files_fasta, self.file_reference_DB)

    def __generate_gene_CDS():
        """'Creates a fasta file containing the whole transcriptome."""
        # should be implemented in an external class

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
        self.file_oligos_DB = os.path.join(
            dir_output,
            f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}",
        )
        file_region_annotation = create_target_region(region, genes)
        if genes is None:
            genes = self.gene_transcript.get_genes_from_annotation()
        self.oligos_DB = self.oligos.generate(
            file_region_annotation, genes, number_batchs, dir_output
        )
        if write:
            self.write_oligos_DB()

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
