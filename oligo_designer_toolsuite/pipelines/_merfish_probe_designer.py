import warnings
from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilter,
    Blastn,
)
from oligo_designer_toolsuite.database import CustomGenomicRegionGenerator, NcbiGenomicRegionGenerator, \
    EnsemblGenomicRegionGenerator
from oligo_designer_toolsuite.database import ReferenceDatabase
from oligo_designer_toolsuite.utils import BaseFtpLoader
from oligo_designer_toolsuite.pipelines import (PrimerProbes,ReadoutProbes,TargetProbes,generate_codebook)
import numpy as np
from Bio.Seq import Seq
import yaml
import os
from pathlib import Path


class MerfishProbeDesigner:
    """
    This class is used to design the final merfish probes.
    param config_file: file that contains the configuration used to make the Merfish probes
    type config_file: str
    """

    def __init__(
            self,
            config_file

    ):
        """Constructor method"""
        # set parameters
        self.config_file = config_file
        with open(self.config_file, 'r') as yaml_file:
            self.config = yaml.safe_load(yaml_file)

        dir_output = os.path.join(os.path.dirname(os.getcwd()), self.config["dir_output"])
        self.dir_output = os.path.join(dir_output, "merfish_sequences")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # Create the transcriptome file
        # If the custom config file is selected
        print("Creating transcriptome")
        if self.config["source"] == "custom":
            self.region_generator = CustomGenomicRegionGenerator(
                annotation_file=self.config["annotation_file"],
                sequence_file=self.config["sequence_file"],
                files_source=self.config["files_source"],
                species=self.config["species"],
                annotation_release=self.config["annotation_release"],
                genome_assembly=self.config["genome_assembly"],
                dir_output=self.dir_output
            )
        # If the Ncbi config file is selected
        elif self.config["source"] == "ncbi":
            self.region_generator = NcbiGenomicRegionGenerator(
                taxon=self.config["taxon"],
                species=self.config["species"],
                annotation_release=self.config["annotation_release"],
                dir_output=self.dir_output
            )
        # If the Ensembl config file is generated
        elif self.config["source"] == "ensembl":
            self.region_generator = EnsemblGenomicRegionGenerator(
                species=self.config["species"],
                annotation_release=self.config["annotation_release"],
                dir_output=self.dir_output
            )

        print("making file_transcriptome")
        self.file_transcriptome = self.region_generator.generate_transcript_reduced_representation(
            include_exon_junctions=True, exon_junction_size=2 * self.config["oligo_length_max"])
        print("Done!")

    def design_merfish_probes(
            self):
        # create target probes
        target_probe_class = TargetProbes(self.config,
                                          self.dir_output,
                                          self.file_transcriptome,
                                          self.region_generator)
        target_probes, file_target_probes = target_probe_class.create_target()
        #calculate number of target probes
        n_probes=0
        genes = list(target_probes.database.keys())
        for gene in genes:
            oligos = list(target_probes.database[gene].keys())
            n_probes+=len(oligos)


        # create primer sequences
        print("Creating Primer Probes")
        primer_probes = PrimerProbes( n_probes,
            self.config,
            self.dir_output,
            self.file_transcriptome,
            self.region_generator)
        primer1, primer2, primer_fasta_file = primer_probes.create_primer()  # return dictionary for primer1 primer2
        print("Creating Primer Probes... Done")

        # create readout probes
        readouts = ReadoutProbes(self.config,
                                self.dir_output,
                                self.file_transcriptome,
                                self.region_generator,
                                primer_fasta_file
                                )
        if (self.config["use_default_readouts"]):
            readout_probes = readouts.get_default_readouts()

        else:
            readout_probes = readouts.create_readouts(self.config["MHD4"]["num_bits"])  # return readout dictionary
        # Create readout_sequences: reverse complement of the readout probes
        readout_sequences = np.zeros_like(readout_probes)
        for i, seq in enumerate(readout_probes):
            readout_sequences[i] = str(Seq(seq).reverse_complement())

        # generate codebook for the genes, leave a percentage of blank barcodes
        database = target_probes.database
        genes = list(database.keys())
        n_genes = len(genes)
        n_blanks = round(n_genes * self.config["percentage_blank"] / 100)
        n_codes = n_genes + n_blanks
        code = generate_codebook(num_seq=n_codes, encoding_scheme=self.config["encoding"])

        # assemble the probes
        primer_idx=0
        for gene_idx, gene in enumerate(genes):
            gene_code = np.asarray(list(code[gene_idx]))
            ones = np.where(gene_code == '1')[0]  # find 1s in the code
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
                encoding_probe = primer1[primer_idx] + encoding_probe + primer2[primer_idx]
                primer_idx+=1

                # put assembled probe in the database
                target_probes.database[gene][oligo_id]["sequence"] = Seq(encoding_probe)

        # save assembled probes in a file
        #assembled_probes_file_database = target_probes.write_database(filename="merfish_assembled_probes.txt")

        # blast against ncRNA
        loader= BaseFtpLoader(dir_output=self.dir_output)
        file_ncRNA= loader. _download_and_decompress(
        ftp_link="ftp.ensembl.org", ftp_directory="/pub/release-109/fasta/homo_sapiens/ncrna/", file_name="Homo_sapiens.GRCh38.ncrna.fa.gz"
        )
        reference_database1 = ReferenceDatabase(
            file_fasta = file_ncRNA,
            files_source = 'Ensembl',
            species = self.region_generator.species,
            annotation_release = self.region_generator.annotation_release,
            genome_assembly = self.region_generator.genome_assembly,
            dir_output=self.dir_output
        )
        reference_database1.load_fasta_into_database()
        dir_specificity1 = os.path.join(self.dir_output, "specificity_temporary1") 
        blastn1 = Blastn(
            dir_specificity=dir_specificity1, 
            word_size=self.config["probe_setup"]["blast1_word_size"],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        specificity_filter1 = SpecificityFilter(filters=[blastn1], write_regions_with_insufficient_oligos=self.config["write_removed_genes"])
        target_probes = specificity_filter1.apply(oligo_database=target_probes, reference_database=reference_database1, n_jobs=self.config["n_jobs"])

        # blast against highly expressed genes
        reference_database2 = ReferenceDatabase(
            file_fasta=self.file_transcriptome,
            files_source=self.region_generator.files_source,
            species=self.region_generator.species,
            annotation_release=self.region_generator.annotation_release,
            genome_assembly=self.region_generator.genome_assembly,
            dir_output=self.dir_output
        )
        reference_database2.load_fasta_into_database()
        if self.config["probe_setup"]["file_highly_expressed_genes"] is None:
            warnings.warn(
                "No file containing the highly expressed genes was provided, all the genes be used."
            )
            genes = None
        else:
            with open(self.config["probe_setup"]["file_highly_expressed_genes"]) as handle:
                lines = handle.readlines()
                genes = [line.rstrip() for line in lines]
            reference_database2.filter_database(genes)

        dir_specificity2 = os.path.join(self.dir_output, "specificity_temporary2")
        blastn2 = Blastn(
            dir_specificity=dir_specificity2,
            word_size=self.config["probe_setup"]["blast2_word_size"],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
        )
        specificity_filter2 = SpecificityFilter(filters=[blastn2], write_regions_with_insufficient_oligos=self.config[
            "write_removed_genes"])
        assembled_probes = specificity_filter2.apply(oligo_database=target_probes,
                                                     reference_database=reference_database2,
                                                     n_jobs=self.config["n_jobs"])

        # save final probes in a file
        #final_probes_file_database = assembled_probes.write_database(filename="merfish_final_probes.txt")

        # create formatted output file
        database = assembled_probes.database
        genes = list(database.keys())
        n_genes=len(genes)
        yaml_dict = {}
        readout_length=self.config["readout_oligo"]["oligo_length_max"]
        primer1_length=self.config["primer_oligo"]["oligo_length_max"]
        primer2_length=primer1_length+18 #T7 promoter has 18nt
        for gene_idx, gene in enumerate(genes):
            yaml_dict[gene] = {}
            oligo_ids = list(database[gene].keys())
            for oligo_idx, oligo_id in enumerate(oligo_ids):

                full_sequence = database[gene][oligo_id]["sequence"]
                primer1= full_sequence[:primer1_length]
                primer2= full_sequence[-primer2_length:]
                readout_seq_1= full_sequence[primer1_length:primer1_length+readout_length]
                target_sequence= full_sequence[primer1_length+readout_length:-(primer2_length+readout_length)]
                readout_seq_2= full_sequence[-(primer2_length+readout_length):-primer2_length]

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
                        yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][
                            key
                        ] = ",".join(str(database_region[oligo_id][key]))

                yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"].update(
                    {
                        "merfish_probe_full_sequence": str(encoding_probe),
                        "readout_probe_1": str(readout_seq_1.reverse_complement()),
                        "readout_probe_2": str(readout_seq_2.reverse_complement()),
                        "primer_sequence_1": str(primer1),
                        "primer_sequence_2": str(primer2),
                        "merfish_barcode_sequence": str(code[gene_idx]),
                        "target_sequence": str(target_sequence),
                        "recognised_mRNA_sequence": str(
                            target_sequence.reverse_complement()
                        ),  
                    }
                )
                for key in [
                    "GC_content",
                    "melting_temperature"
                ]:
                    yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][key] = float(
                        database_region[oligo_id][key]
                    )
                for key in ["length"]:  
                    yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][key] = int(
                        database_region[oligo_id][key]
                    )
        #save 
        final_dir_output = os.path.join(self.dir_output, "final_merfish_probes")
        Path(final_dir_output).mkdir(parents=True, exist_ok=True)
        with open(os.path.join(final_dir_output, "merfish_probes.yml"), "w") as outfile:
            yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)

        #Create order file
        yaml_order = {}
        for region in yaml_dict:
            yaml_order[region] = {}
            for oligo_id in yaml_dict[region]:
                yaml_order[region][oligo_id] = {}
                yaml_order[region][oligo_id]["merfish_probe_full_sequence"] = yaml_dict[
                    region
                ][oligo_id]["merfish_probe_full_sequence"]
                yaml_order[region][oligo_id]["readout_probe_1"] = yaml_dict[
                    region
                ][oligo_id]["readout_probe_1"]
                yaml_order[region][oligo_id]["readout_probe_2"] = yaml_dict[
                    region
                ][oligo_id]["readout_probe_2"]
        with open(
            os.path.join(final_dir_output, "merfish_probes_order.yml"), "w"
        ) as outfile:
            yaml.dump(yaml_order, outfile, default_flow_style=False, sort_keys=False)

        #Create codebook file
        print("Creating Codebook")
        yaml_codebook = {}
        for gene_idx, gene in enumerate(genes):
            yaml_codebook[gene] = code[gene_idx]

        for i in range(n_blanks):
            yaml_codebook[f"blank_barcode_{i + 1}"] = code[n_genes + i]
        with open(os.path.join(final_dir_output, "merfish_codebook.yml"), "w") as outfile:
            yaml.dump(yaml_codebook, outfile, default_flow_style=False, sort_keys=False)
