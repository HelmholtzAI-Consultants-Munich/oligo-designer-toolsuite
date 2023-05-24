import numpy as np
from Bio.Seq import Seq

from pathlib import Path
import os

import yaml

from oligo_designer_toolsuite.utils import generate_codebook
from oligo_designer_toolsuite.database import OligoDatabase


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
                readout_seq_1 = full_sequence[
                    primer1_length : primer1_length + readout_length
                ]
                target_sequence = full_sequence[
                    primer1_length + readout_length : -(primer2_length + readout_length)
                ]
                readout_seq_2 = full_sequence[
                    -(primer2_length + readout_length) : -primer2_length
                ]

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
                        "readout_probe_1": str(
                            readout_seq_1.reverse_complement() + "/3Cy5Sp"
                        ),
                        "readout_probe_2": str(
                            readout_seq_2.reverse_complement() + "/3Cy5Sp"
                        ),
                        "primer_sequence_1": str(primer1),
                        "primer_sequence_2": str(primer2),
                        "merfish_barcode_sequence": str(self.code[gene_idx]),
                        "target_sequence": str(target_sequence),
                        "recognised_mRNA_sequence": str(
                            target_sequence.reverse_complement()
                        ),
                    }
                )
                for key in ["GC_content", "melting_temperature"]:
                    yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][key] = float(
                        database_region[oligo_id][key]
                    )
                for key in ["length"]:
                    yaml_dict[gene][f"{gene}_oligo{oligo_idx+1}"][key] = int(
                        database_region[oligo_id][key]
                    )
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
                yaml_order[region][oligo_id]["merfish_probe_full_sequence"] = yaml_dict[
                    region
                ][oligo_id]["merfish_probe_full_sequence"]
                yaml_order[region][oligo_id]["readout_probe_1"] = yaml_dict[region][
                    oligo_id
                ]["readout_probe_1"]
                yaml_order[region][oligo_id]["readout_probe_2"] = yaml_dict[region][
                    oligo_id
                ]["readout_probe_2"]
        with open(
            os.path.join(final_dir_output, "merfish_probes_order.yml"), "w"
        ) as outfile:
            yaml.dump(yaml_order, outfile, default_flow_style=False, sort_keys=False)

        # Create readout probe file
        yaml_readout = {}
        yaml_readout["Bit"] = "Readout Probe"
        for i in range(num_bits):
            yaml_readout[str(i + 1)] = readout_probes[i] + "/3Cy5Sp"
        with open(
            os.path.join(final_dir_output, "merfish_readout_probes.yml"), "w"
        ) as outfile:
            yaml.dump(yaml_readout, outfile, default_flow_style=False, sort_keys=False)

        # Create codebook file
        yaml_codebook = {}
        for gene_idx, gene in enumerate(genes):
            yaml_codebook[gene] = self.code[gene_idx]

        for i in range(self.n_blanks):
            yaml_codebook[f"blank_barcode_{i + 1}"] = self.code[n_genes + i]
        with open(
            os.path.join(final_dir_output, "merfish_codebook.yml"), "w"
        ) as outfile:
            yaml.dump(yaml_codebook, outfile, default_flow_style=False, sort_keys=False)
