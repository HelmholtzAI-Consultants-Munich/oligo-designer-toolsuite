import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import os
import yaml
from oligo_designer_toolsuite.oligo_property_filter import (
    PropertyFilter,
    MaskedSequences,
    GCContent,
    MeltingTemperatureNN,
    PadlockArms,
    ConsecutiveRepeats,
    GCClamp
)
from oligo_designer_toolsuite.database import ReferenceDatabase, OligoDatabase

from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilter,
    Blastn,
)


class PrimerProbes:

    def __init__(
            self,
            num_seq,
            config_path,
            file_transcriptome,
            region_generator

    ):
        self.num_seq = num_seq
        self.length = 25  # need to trim the sequece
        # self.config_file = os.path.join(os.getcwd(), "tutorials", "configs", "merfish_probe_designer_test.yaml")
        with open(config_path, 'r') as yaml_file:
            self.config = yaml.safe_load(yaml_file)
        self.primer_oligo_config = self.config["primer_oligo"]
        self.file_transcriptome = file_transcriptome
        self.region_generator = region_generator
        self.config_param = self.config["primers_setup"].copy()
        self.GC_min = self.config_param["GC_content_min"]
        self.GC_max = self.config_param["GC_content_max"]
        self.GC_content_filter = GCContent(GC_content_min=self.GC_min, GC_content_max=self.GC_max)
        self.repeat_num_max = ConsecutiveRepeats(3)
        self.lib = ['A', 'C', 'G', 'T']
        self.T7promoter = "TAATACGACTCACTATAG"
        # self.blast_filter = blast_filter
        # self.ref = os.path.basename(reference_DB.file_reference_DB)
        self.oligo_25nt_path = os.path.join(os.getcwd(), "oligo_designer_toolsuite", "experiment_specific", "data",
                                            "bc25mer.240k.fasta")
        self.oligo_25nt_dict = SeqIO.to_dict(SeqIO.parse(self.oligo_25nt_path, "fasta"))
        self.oligo_20nt_dict = {k: v[:-5] for k, v in self.oligo_25nt_dict.items()}

        self.oligo_20nt_path = os.path.join(self.primer_oligo_config["oligo_output"], "bc20mer_oligo.fasta")
        with open(self.oligo_20nt_path, "w") as handle:
            for name, seq in self.oligo_20nt_dict.items():
                handle.write(">" + name + "\n")
                handle.write(str(seq) + "\n")
            print("write bc20mer_oligo done!")

        self.oligo_20nt_DB = None

        self.consecutive_repeats = ConsecutiveRepeats(self.config_param["Repeat_AA_max"])
        self.GC_clamp = GCClamp(self.config_param["GC_clamp_n"])

        # create the list of filters
        # self.melting_temperature

        self.filters = [self.GC_content_filter, self.consecutive_repeats, self.GC_clamp]

        # specificty filters

        # first blast against human transcriptome
        self.dir_specificity1 = os.path.join(self.config["dir_output"], "specificity_temporary1")
        self.blast_filter1 = Blastn(
            dir_specificity=self.dir_specificity1,
            word_size=self.config["primers_blast_setup"]['blast1_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        self.reference_database1 = ReferenceDatabase(
            file_fasta=self.file_transcriptome,
            files_source=self.region_generator.files_source,
            species=self.region_generator.species,
            annotation_release=self.region_generator.annotation_release,
            genome_assembly=self.region_generator.genome_assembly,
            dir_output=self.dir_output
        )
        self.reference_database1.load_fasta_into_database()

        # second blast against 3' end of other primers
        self.dir_specificity2 = os.path.join(self.config["dir_output"], "specificity_temporary2")
        self.blast_filter2 = Blastn(
            dir_specificity=self.dir_specificity2,
            word_size=self.config["primers_blast_setup"]['blast2_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )

        # third blast against 3' end of T7 promoter - Trim T7 to blast word size
        self.dir_specificity3 = os.path.join(self.config["dir_output"], "specificity_temporary3")
        self.blast_filter3 = Blastn(
            dir_specificity=self.dir_specificity3,
            word_size=self.config["primers_blast_setup"]['blast3_word_size'],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )

        T7_dict = dict()
        T7_dict['T7promoter'] = {
            "sequence": Seq(self.T7promoter[-self.config["primers_blast_setup"]['blast3_word_size']:])}
        # create reference DB with fasta file
        fasta_reference_database3 = self.blast_filter3._create_fasta_file(T7_dict, self.dir_specificity3, 'T7')
        self.reference_database3 = ReferenceDatabase(
            file_fasta=fasta_reference_database3)
        self.reference_database3.load_fasta_into_database()
        self.oligo_20nt_DB = PrimerProbes.create_primer_DB()

    def create_primer_DB(self):

        oligo_database = OligoDatabase(
            file_fasta=self.oligo_20nt_path,
            oligo_length_min=self.primer_oligo_config["oligo_length_min"],
            oligo_length_max=self.primer_oligo_config["oligo_length_max"],
            n_jobs=1,
            dir_output=self.primer_oligo_config["oligo_DB_output"],
        )
        return oligo_database.create_database()




    def create_primer(self):

        property_filter = PropertyFilter(filters=self.filters,
                                         write_genes_with_insufficient_oligos=self.config["write_removed_genes"])
        # property filter
        oligo_database = property_filter.apply(oligo_database=self.oligo_20nt_DB, n_jobs=self.config["n_jobs"])

        # specifity filter 1
        specificity_filter1 = SpecificityFilter(filters=[self.blast_filter1],
                                                write_regions_with_insufficient_oligos=self.config[
                                                    "write_removed_genes"])
        oligo_database = specificity_filter1.apply(oligo_database=oligo_database,
                                                   reference_database=self.reference_database1,
                                                   n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="primer_database_specificity_filter_1.txt")

        # specificity filter 2
        # create reference db with trimmed primers
        trimmed_primers = {k: v[-self.config["primers_blast_setup"]['blast2_word_size']:] for k, v in
                           oligo_database.items()}
        # create reference DB with fasta file
        fasta_reference_database2 = self.blast_filter2._create_fasta_file(trimmed_primers, self.dir_specificity2,
                                                                          'primers')
        reference_database2 = ReferenceDatabase(
            file_fasta=fasta_reference_database2)
        reference_database2.load_fasta_into_database()
        specificity_filter2 = SpecificityFilter(filters=[self.blast_filter2],
                                                write_regions_with_insufficient_oligos=self.config[
                                                    "write_removed_genes"])
        oligo_database = specificity_filter2.apply(oligo_database=oligo_database,
                                                   reference_database=reference_database2, n_jobs=self.config["n_jobs"])
        if self.config["write_intermediate_steps"]:
            file_database = oligo_database.write_database(filename="primer_database_specificity_filter2.txt")

        # specificity filter 3
        specificity_filter3 = SpecificityFilter(filters=[self.blast_filter3],
                                                write_regions_with_insufficient_oligos=self.config[
                                                    "write_removed_genes"])
        oligo_database = specificity_filter3.apply(oligo_database=oligo_database,
                                                   reference_database=self.reference_database3,
                                                   n_jobs=self.config["n_jobs"])

        primer_file_database = oligo_database.write_database(filename="primer_database.txt")

        primer1_oligos_dict = {}
        primer1_genes = list(oligo_database.database.keys())[0:self.num_seq]
        primer1_oligo_ids = [list(oligo_database.database[gene].keys())[0] for gene in primer1_genes]
        for gene, oligo_id in zip(primer1_genes, primer1_oligo_ids):
            primer1_oligos_dict[gene] = str(oligo_database.database[gene][oligo_id]["sequence"])

        primer2_oligos_dict = {}
        primer2_genes = list(oligo_database.database.keys())[self.num_seq+1: self.num_seq*2]
        primer2_oligo_ids = [list(oligo_database.database[gene].keys())[0] for gene in primer2_genes]
        for gene, oligo_id in zip(primer2_genes, primer2_oligo_ids):
            primer2_seq = str(oligo_database.database[gene][oligo_id]["sequence"])
            primer2_seq = self.T7promoter + primer2_seq[::-1]
            primer2_oligos_dict[gene] = primer2_seq

        return list(primer1_oligos_dict.values()),list(primer2_oligos_dict.values()), primer_file_database # maybe take half of them for primer1 half for primer2?

