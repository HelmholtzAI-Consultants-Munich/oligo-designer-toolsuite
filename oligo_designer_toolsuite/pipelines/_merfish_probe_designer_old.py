import os
import random
import warnings
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path
from itertools import combinations

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_property_filter import (
    GCClampFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SecondaryStructureFilter,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    Blastn,
    ExactMatches,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._sequence_design import generate_random_sequence
from oligo_designer_toolsuite.sequence_design import MerfishSequence
from oligo_designer_toolsuite.sequence_generator import FtpLoaderEnsembl

from ._base_oligo_designer import BaseOligoDesigner
from ._utils import initialize_parameters


import os
from pathlib import Path

import numpy as np
import yaml
from scipy.spatial.distance import hamming

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.utils import generate_codebook








# TODO: make it more specific (consider GC- fraction for example)
def generate_random_sequence(oligo_length, sequence_alphabet=["A", "T", "G", "C"]):
    seq_letters = np.random.choice(sequence_alphabet, oligo_length, replace=True)
    seq = "".join(seq_letters)
    return Seq(seq)


def generate_binary_sequences(n_seq, n_bit=16, n_one=1):
    """
    This generator randomly generates the codebook from user input, and returns a list of encoding sequences used for merfish.

    Params:
    n_bit: Number of bits in an encoding code
    n_one: Number of ones in an encoding code
    n_seq: Number of encoding codes in a list
    """

    if n_seq > 2**n_bit:
        raise Exception("Error: the number of sequences is tool large!")
    if n_one > n_bit:
        raise Exception("Error: the number of ones is tool large!")

    else:
        sequences = set()
        while len(sequences) < n_seq:
            candidate = format(random.getrandbits(n_bit), f"0{n_bit}b")
            if bin(int(candidate, 2)).count("1") == n_one and candidate not in sequences:
                sequences.add(candidate)

        return list(sequences)


def generate_codebook(num_seq: int, encoding_scheme: str):
    """
    This function generates bit-codes

    return:
    A list of num_seq numbers of MHD2/MHD4 encoding codes

    Params:
    num_seq: Numer of Codes
    encoding_scheme: MHD2 or MHD4
    """

    # initialize an empty list to store the codes
    codes = []
    pass_sequences = []
    if encoding_scheme == "MHD4":
        # loop through all possible combinations of 16 bits
        for i in range(2**16):
            # convert the integer i to a binary string of length 16
            bin_str = format(i, "016b")
            # count the number of ones in the binary string
            num_ones = bin_str.count("1")
            # if the number of ones is exactly 4
            if num_ones == 4:
                # append the binary string to the codes list
                codes.append(bin_str)
        # loop through all pairs of codes in the list

        for i in range(len(codes)):
            for j in range(i + 1, len(codes)):
                if len(pass_sequences) == 140 or len(pass_sequences) == num_seq:
                    break

                # compute the Hamming distance between two codes
                ham_dist = sum(a != b for a, b in zip(codes[i], codes[j]))
                # if the Hamming distance is less than 4
                if ham_dist >= 4:
                    if len(pass_sequences) == 0:
                        pass_sequences.append(codes[i])
                    # remove one of the codes from the list
                    else:
                        ham_list = []
                        for pass_seq in pass_sequences:
                            ham_dist = sum(a != b for a, b in zip(pass_seq, codes[j]))
                            ham_list.append(ham_dist)

                        if all(4 <= ham for ham in ham_list):
                            pass_sequences.append(codes[j])
                            break
            if len(pass_sequences) == 140 or len(pass_sequences) == num_seq:
                break
                # break out of the inner loop

    if encoding_scheme == "MHD2":
        # loop through all possible combinations of 14 bits
        for i in range(2**14):
            # convert the integer i to a binary string of length 16
            bin_str = format(i, "016b")
            # count the number of ones in the binary string
            num_ones = bin_str.count("1")
            # if the number of ones is exactly 4
            if num_ones == 4:
                # append the binary string to the codes list
                codes.append(bin_str)
        # loop through all pairs of codes in the list

        for i in range(len(codes)):
            for j in range(i + 1, len(codes)):
                if len(pass_sequences) == 1001 or len(pass_sequences) == num_seq:
                    break
                # compute the Hamming distance between two codes
                ham_dist = sum(a != b for a, b in zip(codes[i], codes[j]))
                # if the Hamming distance is less than 2
                if ham_dist >= 2:
                    if len(pass_sequences) == 0:
                        pass_sequences.append(codes[i])
                    # remove one of the codes from the list
                    else:
                        ham_list = []
                        for pass_seq in pass_sequences:
                            ham_dist = sum(a != b for a, b in zip(pass_seq, codes[j]))
                            ham_list.append(ham_dist)
                        # check passed sequences
                        if all(2 <= ham for ham in ham_list):
                            pass_sequences.append(codes[j])
                            break
                    # break out of the inner loop
            if len(pass_sequences) == 1001 or len(pass_sequences) == num_seq:
                break
    return pass_sequences


class MerfishProbeDesigner(BaseOligoDesigner):
    """_summary_

    Args:
        BaseOligoDesigner (_type_): _description_
    """

    def design_readout_probes(
        self,
        primer_fasta: str,
        file_barcode_25nt: str,
        percent_identity: float = 80,
        coverage: float = 50,
        strand: str = "plus",
        primer_blast_word_size: int = 11,
        blast2_word_size: int = 14,
        num_sequences: int = 1000,
        readout_probe_length: int = 30,
        num_readouts: int = 16,
        use_default_readouts: bool = False,
        n_jobs: int = 1,
    ):
        """_summary_"""

        if use_default_readouts:
            readout_probes = self._get_default_readouts()

        else:
            (
                potential_readouts,
                file_potential_readouts,
            ) = self._create_potential_readouts_db(file_barcode_25nt, num_sequences, readout_probe_length)
            readout_probes = self._filter_readouts(
                potential_readouts,
                file_potential_readouts,
                num_readouts,
                primer_fasta,
                percent_identity,
                coverage,
                strand,
                primer_blast_word_size,
                blast2_word_size,
                n_jobs,
            )

        # Create readout_sequences: reverse complement of the readout probes
        return [str(Seq(seq).reverse_complement()) for seq in readout_probes]

    def _get_default_readouts(self):
        """
        This function returns 16 validated readout probs from the merfish paper on Homo Sapiens
        """

        default_readout_probes = [
            "CGCAACGCTTGGGACGGTTCCAATCGGATC",
            "CGAATGCTCTGGCCTCGAACGAACGATAGC",
            "ACAAATCCGACCAGATCGGACGATCATGGG",
            "CAAGTATGCAGCGCGATTGACCGTCTCGTT",
            "GCGGGAAGCACGTGGATTAGGGCATCGACC",
            "AAGTCGTACGCCGATGCGCAGCAATTCACT",
            "CGAAACATCGGCCACGGTCCCGTTGAACTT",
            "ACGAATCCACCGTCCAGCGCGTCAAACAGA",
            "CGCGAAATCCCCGTAACGAGCGTCCCTTGC",
            "GCATGAGTTGCCTGGCGTTGCGACGACTAA",
            "CCGTCGTCTCCGGTCCACCGTTGCGCTTAC",
            "GGCCAATGGCCCAGGTCCGTCACGCAATTT",
            "TTGATCGAATCGGAGCGTAGCGGAATCTGC",
            "CGCGCGGATCCGCTTGTCGGGAACGGATAC",
            "GCCTCGATTACGACGGATGTAATTCGGCCG",
            "GCCCGTATTCCCGCTTGCGAGTAGGGCAAT",
        ]

        return default_readout_probes

    # TODO: rename vars
    def _filter_readouts(
        self,
        potential_readouts_db: OligoDatabase,
        file_potential_readouts_db: str,
        num_readouts: int,
        primer_fasta: str,
        percent_identity: float,
        coverage: float,
        strand: str,
        primer_blast_word_size: int,
        blast2_word_size: int,
        n_jobs: int,
    ):
        """
        Function to create the readout probes
        param num_readouts: number of readout probes that should be created
        type num_readouts: int
        """
        if primer_fasta is not None:
            # blast each potential readout probe against the previous build primer probes library
            dir_specificity_primers = os.path.join(self.readout_dir_output, "specificity_temporary_primers")
            exact_matches = ExactMatches()
            blast_filter = Blastn(
                dir_specificity=dir_specificity_primers,
                word_size=primer_blast_word_size,
                percent_identity=percent_identity,
                coverage=coverage,
                strand=strand,
            )
            filters = [exact_matches, blast_filter]
            reference_database = ReferenceDatabase(file_fasta=primer_fasta)
            reference_database.load_fasta_into_database()
            specificity_filter = SpecificityFilter(filters=filters)
            potential_readouts_db = specificity_filter.apply(
                oligo_database=potential_readouts_db,
                reference_database=reference_database,
                n_jobs=n_jobs,
            )
            if self.write_intermediate_steps:
                file_database = potential_readouts_db.write_database(
                    filename="readout_database_specificity_filter1.txt"
                )

        # blast each potential readout probe against the previous built readout probes library
        dir_specificity1 = os.path.join(self.readout_dir_output, "specificity_temporary1")
        blast_filter1 = Blastn(
            dir_specificity=dir_specificity1,
            word_size=primer_blast_word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
        )
        # create reference DB with fasta file
        reference_database1 = ReferenceDatabase(file_potential_readouts_db)
        reference_database1.load_fasta_into_database()
        specificity_filter1 = SpecificityFilter(filters=[blast_filter1])
        potential_readouts_db = specificity_filter1.apply(
            oligo_database=potential_readouts_db,
            reference_database=reference_database1,
            n_jobs=n_jobs,
        )

        if self.write_intermediate_steps:
            file_database = potential_readouts_db.write_database(
                filename="readout_database_specificity_filter1.txt"
            )

        # remove probs with significant homology to members of the transcriptome
        dir_specificity2 = os.path.join(
            self.readout_dir_output, "specificity_temporary2"
        )  # folder where the temporary files will be written

        reference_database2 = ReferenceDatabase(
            file_fasta=self.file_transcriptome_reference,
            metadata=self.metadata,
            dir_output=self.readout_dir_output,
        )
        reference_database2.load_fasta_into_database()
        blast_filter2 = Blastn(
            dir_specificity=dir_specificity2,
            word_size=blast2_word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
        )
        # initialize the specificity filter class
        specificity_filter2 = SpecificityFilter(filters=[blast_filter2])
        # filter the database
        potential_readouts_db = specificity_filter2.apply(
            oligo_database=potential_readouts_db,
            reference_database=reference_database2,
            n_jobs=n_jobs,
        )

        file_database = potential_readouts_db.write_database(filename="readout_database_full.txt")

        return potential_readouts_db.to_sequence_list()[:num_readouts]

    def _create_potential_readouts_db(
        self, file_barcode_25nt: str, num_sequences: int, readout_probe_length: int
    ):
        self.readout_dir_output = os.path.join(self.dir_output, "readout_probes")
        Path(self.readout_dir_output).mkdir(parents=True, exist_ok=True)

        barcode_25nt_sequences = [rec.seq for rec in SeqIO.parse(file_barcode_25nt, "fasta")]
        barcode_25nt_sequences = random.sample(barcode_25nt_sequences, num_sequences)

        augmented_sequences = []
        for sequence in barcode_25nt_sequences:
            random_sequence = generate_random_sequence(readout_probe_length - len(sequence))
            append_beginning = random.choice([True, False])
            if append_beginning:
                augmented_sequences.append(random_sequence + sequence)
            else:
                augmented_sequences.append(sequence + random_sequence)

        augmented_barcodes = {
            f"augmented_bc25mer_{i}": sequence for i, sequence in enumerate(augmented_sequences)
        }
        readout_database = OligoDatabase(dir_output=self.readout_dir_output)
        readout_database.create_database_from_sequences(augmented_barcodes)

        ##### save database #####
        file_database = readout_database.write_fasta_from_database(
            filename="merfish_potential_readout_probes.tsv"
        )

        return readout_database, file_database

    def design_primer_probes(
        self,
        probe_25nt_path: str,
        primer_length_min: int = 19,
        primer_length_max: int = 20,
        GC_content_min: float = 50,
        GC_content_max: float = 65,
        n_repeats: int = 2,
        gc_clamp_n: int = 2,
        blast1_word_size: float = 14,
        blast2_word_size: float = 11,
        blast3_word_size: float = 5,
        percent_identity: float = 80,
        coverage: float = 50,
        strand: str = "plus",
        T7promoter: str = "TAATACGACTCACTATAG",
        n_probes: int = 1,
        n_jobs: int = 1,
    ):
        probe_25nt_path = Path(probe_25nt_path)
        probe_20nt_path = os.path.join(probe_25nt_path.parent.absolute(), "primer_probe.fasta")

        # create primer sequences
        self.primer_dir_output = os.path.join(self.dir_output, "primer_probes")
        Path(self.primer_dir_output).mkdir(parents=True, exist_ok=True)

        primer_database, file_database, T7_dict = self._create_potential_primers_db(
            probe_25nt_path,
            probe_20nt_path,
            n_probes,
            primer_length_min,
            primer_length_max,
            T7promoter,
            blast3_word_size,
            n_jobs,
        )
        primer_database, file_database = self._filter_primers_by_property(
            primer_database,
            GC_content_min,
            GC_content_max,
            n_repeats,
            gc_clamp_n,
            n_jobs,
        )
        primer_database, file_database = self._filter_primers_by_specificity(
            primer_database,
            T7_dict,
            blast1_word_size,
            blast2_word_size,
            blast3_word_size,
            percent_identity,
            coverage,
            strand,
            n_jobs,
        )

        # save fasta file of primers
        primer_file_database = os.path.join(self.primer_dir_output, "primers.fna")
        output = []

        keys = list(primer_database.database.keys())[0 : 2 * n_probes]
        for key in keys:
            oligo = list(primer_database.database[key])[0]
            seq = primer_database.database[key][oligo]["sequence"]
            output.append(SeqRecord(seq, key, "", ""))

        with open(primer_file_database, "w") as handle:
            SeqIO.write(output, handle, "fasta")

        primer1_probes_dict = {}
        primer1_genes = list(primer_database.database.keys())[0:n_probes]
        primer1_probe_ids = [list(primer_database.database[gene].keys())[0] for gene in primer1_genes]
        for gene, probe_id in zip(primer1_genes, primer1_probe_ids):
            primer1_probes_dict[gene] = str(primer_database.database[gene][probe_id]["sequence"])

        primer2_probes_dict = {}
        primer2_genes = list(primer_database.database.keys())[n_probes + 1 : (n_probes * 2) + 1]
        primer2_probe_ids = [list(primer_database.database[gene].keys())[0] for gene in primer2_genes]
        for gene, probe_id in zip(primer2_genes, primer2_probe_ids):
            primer2_seq = str(primer_database.database[gene][probe_id]["sequence"].reverse_complement())
            primer2_seq = T7promoter + primer2_seq
            primer2_probes_dict[gene] = primer2_seq
        return (
            list(primer1_probes_dict.values()),
            list(primer2_probes_dict.values()),
            primer_file_database,
        )  # maybe take half of them for primer1 half for primer2?

    def _create_potential_primers_db(
        self,
        probe_25nt_path,
        probe_20nt_path,
        num_seq,
        primer_length_min,
        primer_length_max,
        T7promoter,
        blast3_word_size,
        n_jobs,
    ):
        probe_25nt_dict = {rec.id: rec.seq for rec in SeqIO.parse(probe_25nt_path, "fasta")}

        # select n random sequences
        n = num_seq * 1000  # start with 1000 times the number of required primers
        if n > 240000:
            n = 240000
        sequences = list(probe_25nt_dict.keys())
        selected_sequences = random.sample(sequences, n)
        probe_20nt_dict = {}
        for key in selected_sequences:
            value = probe_25nt_dict[key]
            probe_20nt_dict[key] = value[:-5]

        with open(probe_20nt_path, "w") as handle:
            for name, seq in probe_20nt_dict.items():
                handle.write(">" + name + "\n")
                handle.write(str(seq) + "\n")

        T7_dict = dict()
        T7_dict["T7promoter"] = {"sequence": Seq(T7promoter[-blast3_word_size:])}

        primer_database = OligoDatabase(
            n_jobs=n_jobs,
            dir_output=self.primer_dir_output,
        )
        primer_database.create_database(
            file_fasta=probe_20nt_path,
            oligo_length_min=primer_length_min,
            oligo_length_max=primer_length_max,
        )

        return primer_database, probe_20nt_path, T7_dict

    def _filter_primers_by_property(
        self,
        primer_database,
        GC_content_min,
        GC_content_max,
        n_repeats,
        gc_clamp_n,
        n_jobs,
    ):
        GC_content_filter = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        consecutive_repeats = HomopolymericRunsFilter(base=["A", "C", "T", "G"], n=n_repeats)

        GC_clamp = GCClampFilter(gc_clamp_n)

        filters = [GC_content_filter, consecutive_repeats, GC_clamp]

        property_filter = PropertyFilter(filters=filters)
        # property filter
        primer_database = property_filter.apply(oligo_database=primer_database, n_jobs=n_jobs)

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = primer_database.write_database(filename="primer_database_property_filter.txt")
        else:
            file_database = ""

        return primer_database, file_database

    def _filter_primers_by_specificity(
        self,
        primer_database,
        T7_dict,
        blast1_word_size,
        blast2_word_size,
        blast3_word_size,
        percent_identity,
        coverage,
        strand,
        n_jobs,
    ):
        dir_specificity1 = os.path.join(self.primer_dir_output, "specificity_temporary_1")
        blast_filter1 = Blastn(
            dir_specificity=dir_specificity1,
            word_size=blast1_word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
        )
        reference_database1 = ReferenceDatabase(
            file_fasta=self.file_transcriptome_reference,
            metadata=self.metadata,
            dir_output=self.primer_dir_output,
        )
        # self.reference_database1.load_fasta_into_database()

        # second blast against 3' end of other primers
        dir_specificity2 = os.path.join(self.primer_dir_output, "specificity_temporary_2")
        blast_filter2 = Blastn(
            dir_specificity=dir_specificity2,
            word_size=blast2_word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
        )

        # third blast against 3' end of T7 promoter - Trim T7 to blast word size
        dir_specificity3 = os.path.join(self.primer_dir_output, "specificity_temporary3")
        blast_filter3 = Blastn(
            dir_specificity=dir_specificity3,
            word_size=blast3_word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
        )

        # create reference DB with fasta file
        # TODO: doesn't seem like the right way to do it
        fasta_reference_database3 = blast_filter3._create_fasta_file(T7_dict, dir_specificity3, "T7")
        reference_database3 = ReferenceDatabase(file_fasta=fasta_reference_database3)
        reference_database3.load_fasta_into_database()

        # specifity filter 1
        specificity_filter1 = SpecificityFilter(filters=[blast_filter1])

        primer_database = specificity_filter1.apply(
            oligo_database=primer_database,
            reference_database=reference_database1,
            n_jobs=n_jobs,
        )
        if self.write_intermediate_steps:
            file_database = primer_database.write_database(
                filename="primer_database_specificity_filter_1.txt"
            )

        # specificity filter 2
        # create reference db with trimmed primers

        trimmed_primers = {}

        primer_genes = list(primer_database.database.keys())
        primer_probe_ids = [list(primer_database.database[gene].keys())[0] for gene in primer_genes]
        for gene, probe_id in zip(primer_genes, primer_probe_ids):
            # Get 3' end sequences
            trimmed_primers[gene] = str(
                primer_database.database[gene][probe_id]["sequence"][-blast2_word_size:]
            )

        # create reference DB with fasta file

        fasta_reference_database2 = os.path.join(dir_specificity2, "oligos_primers.fna")

        with open(fasta_reference_database2, "w") as handle:
            for name, seq in trimmed_primers.items():
                handle.write(">" + name + "\n")
                handle.write(str(seq) + "\n")

        reference_database2 = ReferenceDatabase(file_fasta=fasta_reference_database2)
        reference_database2.load_fasta_into_database()
        specificity_filter2 = SpecificityFilter(filters=[blast_filter2])
        primer_database = specificity_filter2.apply(
            oligo_database=primer_database,
            reference_database=reference_database2,
            n_jobs=n_jobs,
        )
        if self.write_intermediate_steps:
            file_database = primer_database.write_database(filename="primer_database_specificity_filter2.txt")

        # specificity filter 3
        specificity_filter3 = SpecificityFilter(filters=[blast_filter3])
        primer_database = specificity_filter3.apply(
            oligo_database=primer_database,
            reference_database=reference_database3,
            n_jobs=n_jobs,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = primer_database.write_database(filename="primer_database_specificity_filter3.txt")
        else:
            file_database = ""

        return primer_database, file_database

    def filter_assembled_probes_by_specificity(
        self,
        assembled_probes: OligoDatabase,
        file_highly_expressed_genes: str,
        blast_ncrna_word_size: int = 14,
        percent_identity: float = 66,
        coverage: float = 50,
        strand: str = "plus",
        blast_highly_expressed_genes_word_size: int = 17,
        n_jobs: int = 1,
    ):
        # blast against ncRNA Ensembl
        loader = FtpLoaderEnsembl(
            dir_output=self.dir_output,
            species=self.region_generator.species,
            annotation_release=self.region_generator.annotation_release,
        )
        file_ncrna, _, _ = loader.download_files(file_type="fasta", sequence_nature="ncrna")

        reference_ncrna = ReferenceDatabase(
            file_fasta=file_ncrna,
            metadata=self.metadata,
            dir_output=self.dir_output,
        )
        reference_ncrna.load_fasta_into_database()
        dir_specificity_ncrna = os.path.join(self.dir_output, "specificity_temporary_ncrna")
        blastn1 = Blastn(
            dir_specificity=dir_specificity_ncrna,
            word_size=blast_ncrna_word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
        )
        specificity_filter_ncrna = SpecificityFilter(filters=[blastn1])
        assembled_probes = specificity_filter_ncrna.apply(
            oligo_database=assembled_probes,
            reference_database=reference_ncrna,
            n_jobs=n_jobs,
        )

        # blast against highly expressed genes
        reference_highly_expressed_genes = ReferenceDatabase(
            file_fasta=self.file_transcriptome_reference,
            metadata=self.metadata,
            dir_output=self.dir_output,
        )
        reference_highly_expressed_genes.load_fasta_into_database()
        if file_highly_expressed_genes is None:
            warnings.warn(
                "No file containing the highly expressed genes was provided, all the genes be used."
            )
            genes = None
        else:
            with open(file_highly_expressed_genes) as handle:
                lines = handle.readlines()
                genes = [line.rstrip() for line in lines]
            reference_highly_expressed_genes.filter_database(genes)

        dir_specificity_highly_expressed_genes = os.path.join(
            self.dir_output, "specificity_temporary_highly_expressed_genes"
        )
        blastn2 = Blastn(
            dir_specificity=dir_specificity_highly_expressed_genes,
            word_size=blast_highly_expressed_genes_word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
        )
        specificity_filter_highly_expressed_genes = SpecificityFilter(filters=[blastn2])
        assembled_probes = specificity_filter_highly_expressed_genes.apply(
            oligo_database=assembled_probes,
            reference_database=reference_highly_expressed_genes,
            n_jobs=n_jobs,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = assembled_probes.write_database(
                filename="assembled_probe_database_specificity_filter.txt"
            )
        else:
            file_database = ""

        return assembled_probes, file_database

    # Target probes
    def filter_target_probes_by_property(
        self,
        probe_database: OligoDatabase,
        GC_content_min: int = 40,
        GC_content_max: int = 60,
        Tm_min: int = 70,
        Tm_max: int = 80,
        internal_secondary_structures_T: float = 76,
        internal_secondary_structures_threshold_deltaG: float = 0,
        max_repeats_nt: int = 6,
        Tm_parameters_probe: dict = {
            "check": True,
            "strict": True,
            "c_seq": None,
            "shift": 0,
            "nn_table": "DNA_NN3",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "dnac1": 50,
            "dnac2": 0,
            "selfcomp": False,
            "dNTPs": 0,
            "saltcorr": 7,
            "Na": 1000,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
        },
        Tm_correction_param: dict = {
            "DMSO": 0,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
            "fmd": 20,
        },
        n_jobs: int = 1,
    ):
        ##### preprocess melting temperature params #####
        Tm_parameters_probe["nn_table"] = getattr(mt, Tm_parameters_probe["nn_table"])
        Tm_parameters_probe["tmm_table"] = getattr(mt, Tm_parameters_probe["tmm_table"])
        Tm_parameters_probe["imm_table"] = getattr(mt, Tm_parameters_probe["imm_table"])
        Tm_parameters_probe["de_table"] = getattr(mt, Tm_parameters_probe["de_table"])

        melting_temperature = MeltingTemperatureNNFilter(
            Tm_min=Tm_min,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_correction_param,
        )
        consecutive_repeats = HomopolymericRunsFilter(base=["A", "C", "T", "G"], n=max_repeats_nt)
        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        secondary_structure = SecondaryStructureFilter(
            T=internal_secondary_structures_T,
            DG_thr=internal_secondary_structures_threshold_deltaG,
        )
        masked_sequences = HardMaskedSequenceFilter()

        # create the list of filters
        filters = [
            masked_sequences,
            consecutive_repeats,
            gc_content,
            melting_temperature,
            secondary_structure,
        ]

        # initialize the property filter class
        property_filter = PropertyFilter(filters=filters)
        # filter the database
        probe_database = property_filter.apply(oligo_database=probe_database, n_jobs=n_jobs)

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = probe_database.write_database(filename="oligo_database_property_filter.txt")
        else:
            file_database = ""

        return probe_database, file_database

    # Done
    def filter_target_probes_by_specificity(
        self,
        probe_database: OligoDatabase,
        word_size: int = 15,
        percent_identity: float = 80,
        coverage: float = 50,
        strand: str = "plus",
        n_jobs: int = 1,
    ):
        # Specificity filters to remove probes with more than 1 RNA species target
        dir_specificity = os.path.join(
            self.dir_output, "specificity_temporary"
        )  # folder where the temporary files will be written

        (
            probe_length_min,
            probe_length_max,
        ) = self._get_probe_length_min_max_from_database(probe_database.database)

        self.file_transcriptome_reference = self.region_generator.generate_transcript_reduced_representation(
            include_exon_junctions=True, exon_junction_size=2 * probe_length_max
        )

        reference_database = ReferenceDatabase(
            file_fasta=self.file_transcriptome_reference,
            metadata=self.metadata,
            dir_output=self.dir_output,
        )

        # intialize the filter classes
        exact_matches = ExactMatches(dir_specificity=dir_specificity)
        blastn = Blastn(
            dir_specificity=dir_specificity,
            word_size=word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
        )
        filters = [exact_matches, blastn]
        # initialize the specificity filter class
        specificity_filter = SpecificityFilter(filters=filters)
        # filter the database
        probe_database = specificity_filter.apply(
            oligo_database=probe_database,
            reference_database=reference_database,
            n_jobs=n_jobs,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = probe_database.write_database(filename="probe_database_specificity_filter.txt")
        else:
            file_database = ""

        return probe_database, file_database

    # Done
    def filter_cross_hybridization_targets(
        self,
        probe_database: OligoDatabase,
        word_size: int = 15,
        percent_identity_ch: int = 66,
        coverage: float = 50,
        n_jobs: int = 1,
    ):
        # TODO: check cross hybridization MT
        # Specificity filter to remove cross hybridization targets
        targets_fasta = probe_database.write_fasta_from_database(filename="target_probes_init")
        reference_database = ReferenceDatabase(file_fasta=targets_fasta)
        reference_database.load_fasta_into_database()
        dir_specificity = os.path.join(
            self.dir_output, "specificity_temporary_cross_hybridization"
        )  # folder where the temporary files will be written
        # intialize the filter classes
        blastn = Blastn(
            dir_specificity=dir_specificity,
            word_size=word_size,
            percent_identity=percent_identity_ch,
            coverage=coverage,
            strand="minus",
        )
        filters = [blastn]
        # initialize the specificity filter class
        specificity_filter = SpecificityFilter(filters=filters)
        # filter the database
        probe_database = specificity_filter.apply(
            oligo_database=probe_database,
            reference_database=reference_database,
            n_jobs=n_jobs,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = probe_database.write_database(filename="merfish_target_probes.txt")
        else:
            file_database = ""

        return probe_database, file_database


def main():
    """
    TODO
    """
    parser = ArgumentParser(
        prog="MERFISH Probe Designer",
        usage="merfish_probe_designer [options]",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )

    config = initialize_parameters(parser, exp_name="merfish")

    # with open(merfish_config_path, "r") as handle:
    #     config = yaml.safe_load(handle)

    dir_output = os.path.abspath(config["output"])
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    ##### Initialize ProbeDesigner Class #####
    probe_designer = MerfishProbeDesigner(dir_output=dir_output, log_name="merfish_probe_designer")

    ##### load annotations #####
    probe_designer.load_annotations(source=config["source"], source_params=config["source_params"])

    ##### read the genes file #####
    if config["file_genes"] is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This chioce can use a lot of resources."
        )
        genes = None
    else:
        with open(config["file_genes"]) as handle:
            lines = handle.readlines()
            genes = [line.rstrip() for line in lines]

    ##### create target probe database #####
    target_probe_database, file_database = probe_designer.create_probe_database(
        genes=genes,
        region=config["region"],
        probe_length_min=config["target_probe"]["probe_length_min"],
        probe_length_max=config["target_probe"]["probe_length_max"],
        min_probes_per_gene=config["min_probes_per_gene"],
        n_jobs=config["n_jobs"],
    )

    ##### filter target probes by property #####
    (
        target_probe_database,
        file_database,
    ) = probe_designer.filter_target_probes_by_property(
        probe_database=target_probe_database,
        GC_content_min=config["targets_setup"]["GC_content_min"],
        GC_content_max=config["targets_setup"]["GC_content_max"],
        Tm_min=config["targets_setup"]["Tm_min"],
        Tm_max=config["targets_setup"]["Tm_max"],
        internal_secondary_structures_T=config["targets_setup"]["internal_secondary_structures_T"],
        internal_secondary_structures_threshold_deltaG=config["targets_setup"][
            "internal_secondary_structures_threshold_deltaG"
        ],
        max_repeats_nt=config["targets_setup"]["max_repeats_nt"],
        Tm_parameters_probe=config["Tm_parameters_probe"],
        Tm_correction_param=config["Tm_correction_parameters"],
        n_jobs=config["n_jobs"],
    )

    ##### filter target probes by specificity #####
    (
        target_probe_database,
        file_database,
    ) = probe_designer.filter_target_probes_by_specificity(
        probe_database=target_probe_database,
        word_size=config["targeting_sequences_setup"]["word_size"],
        percent_identity=config["blast_percent_identity"],
        coverage=config["blast_coverage"],
        strand=config["strand"],
        n_jobs=config["n_jobs"],
    )

    ##### filter cross hybridization targets #####
    (
        target_probe_database,
        file_database,
    ) = probe_designer.filter_cross_hybridization_targets(
        probe_database=target_probe_database,
        word_size=config["targeting_sequences_setup"]["word_size"],
        percent_identity_ch=config["targeting_sequences_setup"]["percent_identity_ch"],
        coverage=50,
        n_jobs=config["n_jobs"],
    )

    ##### design primer probes #####
    primer1, primer2, primer_file_database = probe_designer.design_primer_probes(
        probe_25nt_path=config["primer_probe"]["file_bc25mer"],
        primer_length_min=config["primer_probe"]["probe_length_min"],
        primer_length_max=config["primer_probe"]["probe_length_min"],
        GC_content_min=config["primers_setup"]["GC_content_min"],
        GC_content_max=config["primers_setup"]["GC_content_max"],
        n_repeats=config["primers_setup"]["max_repeats_nt"],
        gc_clamp_n=config["primers_setup"]["GC_clamp_n"],
        blast1_word_size=config["primers_blast_setup"]["blast1_word_size"],
        blast2_word_size=config["primers_blast_setup"]["blast2_word_size"],
        blast3_word_size=config["primers_blast_setup"]["blast3_word_size"],
        percent_identity=config["blast_percent_identity"],
        coverage=config["blast_coverage"],
        strand=config["strand"],
        n_jobs=config["n_jobs"],
    )

    ##### design primer probes #####
    encoding = config["encoding"]
    readout_probes = probe_designer.design_readout_probes(
        primer_file_database,
        config["readout_probe"]["file_bc25mer"],
        percent_identity=config["blast_percent_identity"],
        coverage=config["blast_coverage"],
        strand=config["strand"],
        primer_blast_word_size=config["readout_blast_setup"]["blast1_word_size"],
        blast2_word_size=config["readout_blast_setup"]["blast2_word_size"],
        readout_probe_length=config["readout_probe"]["probe_length_max"],
        num_readouts=config[encoding]["num_bits"],
        use_default_readouts=config["use_default_readouts"],
        n_jobs=config["n_jobs"],
    )

    merfish_sequence_designer = MerfishSequence(
        n_genes=len(genes),
        encoding_scheme=encoding,
        percentage_blank=config["percentage_blank_barcodes"],
        dir_output=dir_output,
    )

    assembled_probes = merfish_sequence_designer.assemble_probes(
        target_probe_database, readout_probes, primer1, primer2
    )

    merfish_sequence_designer.write_final_probes(
        assembled_probes,
        readout_probes,
        readout_probe_length=config["readout_probe"]["probe_length_max"],
        primer_probe_length=config["primer_probe"]["probe_length_max"],
        num_bits=config[encoding]["num_bits"],
    )


if __name__ == "__main__":
    main()


n_genes = 3
n_oligos_per_gene = 1
regions = list(oligo_database.database.keys())
reduced_oligo_db = {
    key: {
        inner_key: oligo_database.database[key][inner_key]
        for inner_key in list(oligo_database.database[key].keys())
    }
    for key in regions[:n_genes]
}

reduced_oligo_db


class MerfishSequence:
    def __init__(
        self,
        n_genes,
        encoding_scheme: str = "MHD4",
        percentage_blank: int = 10,
        dir_output: str = "output",
    ):
        self.n_blanks = round(n_genes * percentage_blank / 100)
        n_codes = n_genes + self.n_blanks
        self.code = generate_codebook(num_seq=n_codes, encoding_scheme=encoding_scheme)
        self.dir_output = dir_output

    def assemble_probes(
        self,
        target_probes: OligoDatabase,
        readout_probes: list[str],
        primer1: list[str],
        primer2: list[str],
    ):
        readout_sequences = np.zeros_like(readout_probes)
        for i, seq in enumerate(readout_probes):
            readout_sequences[i] = str(Seq(seq).reverse_complement())

        # generate codebook for the genes, leave a percentage of blank barcodes
        database = target_probes.database
        genes = list(database.keys())

        # assemble the probes
        for gene_idx, gene in enumerate(genes):
            gene_code = np.asarray(list(self.code[gene_idx]))
            ones = np.where(gene_code == "1")[0]  # find 1s in the code
            oligo_ids = list(database[gene].keys())
            for oligo_id_idx, oligo_id in enumerate(oligo_ids):
                # randomly choose 2 1s from the code to be readout probes for this oligo
                readout_idx = np.random.choice(ones, 2, replace=False)
                readout_seq_1 = readout_sequences[readout_idx[0]]
                readout_seq_2 = readout_sequences[readout_idx[1]]

                # target sequence is the reverse complement of the target probe
                target_mRNA = database[gene][oligo_id]["sequence"]
                target_sequence = str(target_mRNA.reverse_complement()).upper()

                # place readout sequences on either side of the targeting sequence
                encoding_probe = readout_seq_1 + target_sequence + readout_seq_2

                # add primers
                encoding_probe = primer1[0] + encoding_probe + primer2[0]

                # put assembled probe in the database
                target_probes.database[gene][oligo_id]["sequence"] = Seq(encoding_probe)

        # save assembled probes in a file
        # assembled_probes_file_database = target_probes.write_database(filename="merfish_assembled_probes.txt")
        return target_probes

    def write_final_probes(
        self,
        assembled_probes,
        readout_probes,
        readout_probe_length,
        primer_probe_length,
        num_bits,
    ):
        # num_bits = self.config[self.config["encoding"]]["num_bits"]
        # create formatted output file
        database = assembled_probes.database
        genes = list(database.keys())
        n_genes = len(genes)
        yaml_dict = {}
        readout_length = readout_probe_length
        primer1_length = primer_probe_length
        primer2_length = primer1_length + 18  # T7 promoter has 18nt
        for gene_idx, gene in enumerate(genes):
            yaml_dict[gene] = {}
            oligo_ids = list(database[gene].keys())
            for oligo_idx, oligo_id in enumerate(oligo_ids):
                full_sequence = database[gene][oligo_id]["sequence"]
                primer1 = full_sequence[:primer1_length]
                primer2 = full_sequence[-primer2_length:]
                readout_seq_1 = full_sequence[primer1_length : primer1_length + readout_length]
                target_sequence = full_sequence[
                    primer1_length + readout_length : -(primer2_length + readout_length)
                ]
                readout_seq_2 = full_sequence[-(primer2_length + readout_length) : -primer2_length]

                yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"] = {}
                yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"]["id"] = oligo_id
                yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"]["gene"] = gene
                database_region = database[gene]
                for key in [
                    "additional_information_fasta",
                    "chromosome",
                    "start",
                    "end",
                    "strand",
                ]:
                    if len(database_region[oligo_id][key]) == 1:
                        yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][key] = str(
                            database_region[oligo_id][key][0]
                        )
                    else:
                        yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][key] = ",".join(
                            str(database_region[oligo_id][key])
                        )

                yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"].update(
                    {
                        "merfish_probe_full_sequence": str(full_sequence),
                        "readout_probe_1": str(readout_seq_1.reverse_complement() + "/3Cy5Sp"),
                        "readout_probe_2": str(readout_seq_2.reverse_complement() + "/3Cy5Sp"),
                        "primer_sequence_1": str(primer1),
                        "primer_sequence_2": str(primer2),
                        "merfish_barcode_sequence": str(self.code[gene_idx]),
                        "target_sequence": str(target_sequence),
                        "recognised_mRNA_sequence": str(target_sequence.reverse_complement()),
                    }
                )
                for key in ["GC_content", "melting_temperature"]:
                    yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][key] = float(database_region[oligo_id][key])
                for key in ["length"]:
                    yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][key] = int(database_region[oligo_id][key])
        # save
        final_dir_output = os.path.join(self.dir_output, "final_merfish_probes")
        Path(final_dir_output).mkdir(parents=True, exist_ok=True)
        with open(os.path.join(final_dir_output, "merfish_probes.yml"), "w") as outfile:
            yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)

        # Create order file
        yaml_order = {}
        for region in yaml_dict:
            yaml_order[region] = {}
            for oligo_id in yaml_dict[region]:
                yaml_order[region][oligo_id] = {}
                yaml_order[region][oligo_id]["merfish_probe_full_sequence"] = yaml_dict[region][oligo_id][
                    "merfish_probe_full_sequence"
                ]
                yaml_order[region][oligo_id]["readout_probe_1"] = yaml_dict[region][oligo_id][
                    "readout_probe_1"
                ]
                yaml_order[region][oligo_id]["readout_probe_2"] = yaml_dict[region][oligo_id][
                    "readout_probe_2"
                ]
        with open(os.path.join(final_dir_output, "merfish_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_order, outfile, default_flow_style=False, sort_keys=False)

        # Create readout probe file
        yaml_readout = {}
        yaml_readout["Bit"] = "Readout Probe"
        for i in range(num_bits):
            yaml_readout[str(i + 1)] = readout_probes[i] + "/3Cy5Sp"
        with open(os.path.join(final_dir_output, "merfish_readout_probes.yml"), "w") as outfile:
            yaml.dump(yaml_readout, outfile, default_flow_style=False, sort_keys=False)

        # Create codebook file
        yaml_codebook = {}
        for gene_idx, gene in enumerate(genes):
            yaml_codebook[gene] = self.code[gene_idx]

        for i in range(self.n_blanks):
            yaml_codebook[f"blank_barcode_{i + 1}"] = self.code[n_genes + i]
        with open(os.path.join(final_dir_output, "merfish_codebook.yml"), "w") as outfile:
            yaml.dump(yaml_codebook, outfile, default_flow_style=False, sort_keys=False)