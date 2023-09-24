############################################
# imports
############################################

import os
import re
import gzip
import shutil
import itertools

import pandas as pd

from Bio import SeqIO
from ftplib import FTP, error_perm
from pathlib import Path


############################################
# FTP Classes
############################################


class BaseFtpLoader:
    """
    Base class for downloading annotations from different FTP servers.

    :param dir_output: Path to directory where downloaded target file is saved
    :type dir_output: str
    """

    def __init__(self, dir_output: str) -> None:
        """Constructor method"""
        self.dir_output = dir_output

    def _download(self, ftp_link: str, ftp_directory: str, file_name: str):
        """
        Download file from ftp server.

        :param ftp_link: Link to ftp server.
        :type ftp_link: str
        :param ftp_directory: Directory path of target files on FTP server.
        :type ftp_directory: str
        :param file_name: Name of target file.
        :type file_name: str
        :return: Path to downloaded file.
        :rtype: str
        """
        ftp = FTP(ftp_link)
        ftp.login()  # login to ftp server
        ftp.cwd(ftp_directory)  # move to directory

        files = ftp.nlst()
        file_output = None

        for file in files:
            if re.match(file_name, file):
                file_output = os.path.join(self.dir_output, file)
                ftp.retrbinary("RETR " + file, open(file_output, "wb").write)

        ftp.quit()

        return file_output

    def _decompress_gzip(self, file_gzip: str):
        """
        Decompress zip files.

        :param file_gzip: Path to zipped file.
        :type file_gzip: str
        :return: Path to unzipped file.
        :rtype: str
        """
        file_output = file_gzip.split(".gz")[0]
        with gzip.open(file_gzip, "rb") as f_in:
            with open(file_output, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(file_gzip)

        return file_output

    def _download_and_decompress(self, ftp_link: str, ftp_directory: str, file_name: str):
        """Download genome sequence from ftp server and unzip file.

        :param ftp_link: Link to ftp server.
        :type ftp_link: str
        :param ftp_directory: Directory path of target files on FTP server.
        :type ftp_directory: str
        :param file_name: Name of target file.
        :type file_name: str
        :return: Path to downloaded file.
        :rtype: str
        """
        file_download = self._download(ftp_link, ftp_directory, file_name)
        file_unzipped = self._decompress_gzip(file_download)

        return file_unzipped

    def _check_file_type(self, file_type: str):
        """Check if file type is supported by FTP download.

        :param file_type: Type of target file.
        :type file_type: str {'gff', 'gtf', 'fasta'}
        """
        valid_file_types = ["gff", "gtf", "fasta"]

        if file_type not in valid_file_types:
            if isinstance(file_type, str):
                raise Exception(
                    f"An invalid file type name is used as input. Accepted file types: {valid_file_types}"
                )
            else:
                raise TypeError(
                    f"file_type should be a string. Accepted choices for file_type: {valid_file_types}"
                )


class FtpLoaderEnsembl(BaseFtpLoader):
    """
    Class for downloading annotations from Ensembl, inheriting from BaseFtpLoader.

    :param dir_output: Path to directory where downloaded target file is saved
    :type dir_output: str
    :param species: species name in ensemble download format, e.g. 'homo_sapiens' for human; see http://ftp.ensembl.org/pub/release-108/gtf/ for available species names
    :type species: str
    :param annotation_release: release number of annotation, e.g. 'release-108' or 'current' to use most recent annotation release. Check out release numbers for ensemble at ftp.ensembl.org/pub/
    :type annotation_release: str
    """

    def __init__(self, dir_output: str, species: str, annotation_release: str) -> None:
        """Constructor method"""
        super().__init__(dir_output)
        self.species = species
        self.annotation_release = annotation_release
        self.assembly_name = None
        self.assembly_name_placeholder = "[^\\.]*"

        self.ftp_link = "ftp.ensembl.org"

        self.file_type_folder = {"gff": "gff3", "gtf": "gtf", "fasta": "fasta"}

        self.file_type_ending = {
            "gff": "gff3.gz",
            "gtf": "gtf.gz",
            "fasta": "dna_sm.primary_assembly.fa.gz",  # soft-masked version of the genome
        }

    def download_files(self, file_type: str, sequence_nature: str = "dna"):
        """
        Download 'file_type' file (e.g. GFF, GTF or fasta) from Ensembl and unzip file.

        :param file_type: Type of target file.
        :type file_type: str {'gff', 'gtf', 'fasta'}
        :param sequence_nature: Nature of the fasta sequence.
        :type sequence_nature: str {'dna', 'ncrna'}
        :return: Path to downloaded file, annotation release and genome assembly name.
        :rtype: tuple of str
        """
        # TODO: find a better solution for sequence_nature
        ftp_directory, ftp_file = self._get_params(file_type, sequence_nature)
        dowloaded_file = self._download_and_decompress(self.ftp_link, ftp_directory, ftp_file)

        self.assembly_name = re.search("\\.([^\\.]*)\\.", Path(dowloaded_file).name).group().replace(".", "")

        return dowloaded_file, self.annotation_release, self.assembly_name

    def _get_params(self, file_type: str, sequence_nature: str):
        """
        Get directory and file name for specified file type from Ensembl server.

        :param file_type: Type of target file.
        :type file_type: str {'gff', 'gtf', 'fasta'}
        :param sequence_nature: Nature of the fasta sequence.
        :type sequence_nature: str {'dna', 'ncrna'}
        :return: ftp directory and name of target file on Ensembl server.
        :rtype: tuple of str
        """
        self._check_file_type(file_type)

        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        if self.annotation_release == "current":
            file_readme = self._download(self.ftp_link, "pub/", "current_README")
            with open(file_readme, "r") as handle:
                for line in handle:
                    if line.startswith("Ensembl Release"):
                        self.annotation_release = line.strip().split(" ")[2]
            os.remove(file_readme)

        if file_type.casefold() == "fasta".casefold():
            ftp_directory = f"pub/release-{self.annotation_release}/{self.file_type_folder[file_type]}/{self.species}/{sequence_nature}/"
            if sequence_nature == "dna":
                ftp_file = f"{self.species.capitalize()}.{self.assembly_name_placeholder}.{self.file_type_ending[file_type]}"
            else:
                ftp_file = (
                    f"{self.species.capitalize()}.{self.assembly_name_placeholder}.{sequence_nature}.fa.gz"
                )
        else:
            ftp_directory = (
                f"pub/release-{self.annotation_release}/{self.file_type_folder[file_type]}/{self.species}/"
            )
            ftp_file = f"{self.species.capitalize()}.{self.assembly_name_placeholder}.{self.annotation_release}.{self.file_type_ending[file_type]}"

        return ftp_directory, ftp_file


class FtpLoaderNCBI(BaseFtpLoader):
    """
    Class for downloading annotations from NCBI, inheriting from BaseFtpLoader.

    :param dir_output: Path to directory where downloaded target file is saved
    :type dir_output: str
    :param taxon: taxon of the species, valid taxa are: archaea, bacteria, fungi, invertebrate, mitochondrion, plant, plasmid, plastid, protozoa, vertebrate_mammalian, vertebrate_other, viral
    :type taxon: str {'archaea', 'bacteria', 'fungi', 'invertebrate', 'mitochondrion', 'plant', 'plasmid', 'plastid', 'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral'}
    :param species: species name in NCBI download format, e.g. 'Homo_sapiens' for human; see https://ftp.ncbi.nlm.nih.gov/genomes/refseq/ for available species name
    :type species: str
    :param annotation_release: release number of annotation e.g. '109' or '109.20211119'  or 'current' to use most recent annotation release. Check out release numbers for NCBI at ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/.
    :type annotation_release: str
    """

    def __init__(self, dir_output: str, taxon: str, species: str, annotation_release: str) -> None:
        """Constructor method"""
        super().__init__(dir_output)
        self.taxon = taxon
        self.species = species
        self.annotation_release = annotation_release
        self.assembly_name = None
        self.assembly_accession = None

        self.ftp_link = "ftp.ncbi.nlm.nih.gov"

        self.file_type_ending = {
            "gff": "genomic.gff.gz",
            "gtf": "genomic.gtf.gz",
            "fasta": "genomic.fna.gz",  # soft-masked version of the genome
        }

        self.file_type_function = {
            "gff": self._map_chr_names_gene_annotation,
            "gtf": self._map_chr_names_gene_annotation,
            "fasta": self._map_chr_names_genome_sequence,
        }

    def download_files(self, file_type: str):
        """
        Download 'file_type' file (e.g. GFF, GTF or fasta) from NCBI and unzip file.
        Map chromosome annotation to Ref-Seq accession number.

        :param file_type: Type of target file.
        :type file_type: str {'gff', 'gtf', 'fasta'}
        :return: Path to downloaded file, annotation release and genome assembly name.
        :rtype: tuple of str
        """
        ftp_directory, ftp_file, ftp_file_chr_mapping = self._get_params(file_type)

        mapping = self._download_mapping_chr_names(ftp_directory, ftp_file_chr_mapping)
        dowloaded_file = self._download_and_decompress(self.ftp_link, ftp_directory, ftp_file)

        self.file_type_function[file_type](dowloaded_file, mapping)

        return dowloaded_file, self.annotation_release, self.assembly_name

    def _get_params(self, file_type: str):
        """
        Get directory and file name for specified file type from NCBI server.

        :param file_type: Type of target file.
        :type file_type: str {'gff', 'gtf', 'fasta'}
        :return: ftp directory, name of target file and chromosome mapping on NCBI server.
        :rtype: tuple of str
        """
        self._check_file_type(file_type)

        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ftp_directory = "genomes/refseq/" + self.taxon + "/" + self.species + "/annotation_releases/"

        if self.annotation_release == "current":
            ftp_directory = ftp_directory + "current/"
            ftp = FTP(
                self.ftp_link
            )  # inside current dir there is the directory containing the annotation release
            ftp.login()  # login to ftp server
            ftp.cwd(ftp_directory)  # move to directory
            self.annotation_release = ftp.nlst()[0]
            ftp.quit()

        ftp_directory = ftp_directory + f"{self.annotation_release}/"

        file_readme = self._download(self.ftp_link, ftp_directory, f"README_.*{self.annotation_release}")
        with open(file_readme, "r") as handle:
            for line in handle:
                if line.startswith("ASSEMBLY NAME:"):
                    self.assembly_name = line.strip().split("\t")[1]
                if line.startswith("ASSEMBLY ACCESSION:"):
                    self.assembly_accession = line.strip().split("\t")[1]
                    break
        os.remove(file_readme)

        # we need this check here because NCBI changed the folder structure for releases > 110
        ftp = FTP(self.ftp_link)
        ftp.login()
        try:
            ftp.cwd(ftp_directory + f"{self.assembly_accession}_{self.assembly_name}")  # move to directory
            ftp_directory = ftp_directory + f"{self.assembly_accession}_{self.assembly_name}"
        except error_perm as msg:
            ftp_directory = ftp_directory
        ftp.quit()

        ftp_file = f"{self.assembly_accession}_{self.assembly_name}_{self.file_type_ending[file_type]}"
        ftp_file_chr_mapping = f"{self.assembly_accession}_{self.assembly_name}_assembly_report.txt"

        return ftp_directory, ftp_file, ftp_file_chr_mapping

    def _download_mapping_chr_names(self, ftp_directory: str, ftp_file_chr_mapping: str):
        """
        Download file with mapping of chromosome names between GenBank and Ref-Seq accession number
        from ftp server and create a mapping dictionary.

        :param ftp_directory: Directory path of target files on FTP server.
        :type ftp_directory: str
        :param ftp_file_chr_mapping: Name of file on ftp server that contains chromosome mapping information
        :type ftp_file_chr_mapping: str
        :return: Dictionary with mapping of chromsome names from GenBank to Ref-Seq.
        :rtype: dict
        """
        file_mapping = self._download(self.ftp_link, ftp_directory, ftp_file_chr_mapping)

        # skip comment lines but keep last comment line for header
        with open(file_mapping) as handle:
            *_comments, names = itertools.takewhile(lambda line: line.startswith("#"), handle)
            names = names[1:].split()

        assembly_report = pd.read_table(file_mapping, names=names, sep="\t", comment="#")

        mapping_chromosome = assembly_report[assembly_report["Sequence-Role"] == "assembled-molecule"]
        mapping_chromosome = pd.Series(
            mapping_chromosome["Sequence-Name"].values,
            index=mapping_chromosome["RefSeq-Accn"],
        ).to_dict()

        mapping_scaffolds = assembly_report[assembly_report["Sequence-Role"] != "assembled-molecule"]
        mapping_scaffolds = pd.Series(
            mapping_scaffolds["GenBank-Accn"].values,
            index=mapping_scaffolds["RefSeq-Accn"],
        ).to_dict()

        mapping = mapping_chromosome
        mapping.update(mapping_scaffolds)

        return mapping

    def _map_chr_names_gene_annotation(self, ftp_file: str, mapping: dict):
        """
        Process gene annotation file downloaded from NCBI: map chromosome annotation to Ref-Seq.

        :param ftp_file: Path to gene annotation file (gtf/gff).
        :type ftp_file: str
        :param mapping: Chromosome mapping dictionary (GenBank to Ref-Seq).
        :type mapping: dict
        """
        file_tmp = os.path.join(self.dir_output, "temp.gtf")

        # write comment lines to new file
        with open(file_tmp, "w") as handle_out:
            with open(ftp_file) as handle_in:
                *_comments, names = itertools.takewhile(lambda line: line.startswith("#"), handle_in)
                handle_out.write(names)

            # read gtf file without comment lines
            gene_annotation = pd.read_table(
                ftp_file,
                names=[
                    "seqid",
                    "source",
                    "type",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "phase",
                    "attributes",
                ],
                sep="\t",
                comment="#",
            )

            # replace ncbi with genbank chromosome annotation
            gene_annotation["seqid"] = gene_annotation["seqid"].map(mapping)
            gene_annotation.dropna(inplace=True)  # drop if no mapping exists

            gene_annotation.to_csv(handle_out, sep="\t", header=False, index=False)
        os.replace(file_tmp, ftp_file)

    def _map_chr_names_genome_sequence(self, ftp_file: str, mapping: dict):
        """
        Process genome sequence file downloaded from NCBI: map chromosome annotation to Ref-Seq.

        :param ftp_file: Path to genome sequence file (e.g. fasta).
        :type ftp_file: str
        :param mapping: Chromosome mapping dictionary (GenBank to Ref-Seq).
        :type mapping: dict
        """
        file_tmp = os.path.join(self.dir_output, "temp.fna")

        with open(file_tmp, "w") as handle:
            for chromosome_sequnece in SeqIO.parse(ftp_file, "fasta"):
                accession_number = chromosome_sequnece.id
                if accession_number in mapping:
                    chromosome_sequnece.id = mapping[accession_number]
                    chromosome_sequnece.name = mapping[accession_number]
                    chromosome_sequnece.description = chromosome_sequnece.description.replace(
                        accession_number, mapping[accession_number]
                    )
                    SeqIO.write(chromosome_sequnece, handle, "fasta")
                else:
                    self.logging.info("No mapping for accession number: {}".format(accession_number))

        os.replace(file_tmp, ftp_file)
