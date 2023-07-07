import numpy as np
from Bio.Seq import Seq

from ..utils import create_seqfish_plus_barcodes
import os
from pathlib import Path
import yaml


class SeqfishProbesCreator:
    """
    This class is used to assemble probes using primary probes, readout probes and barcodes, that were designed for each gene.
    """

    def __init__(self, dir_output: str):
        self.dir_output = dir_output

    def assemble_probes(
        self,
        oligo_database,
        readout_database,
        n_pseudocolors,
        seed,
    ):
        barcodes = create_seqfish_plus_barcodes(
            n_pseudocolors=n_pseudocolors,
            seed=seed,
            num_genes=len(oligo_database.database.keys()),
        )
        readout_sequences = readout_database.to_sequence_list()
        for i in oligo_database.keys():
            barcode = barcodes[i]
            left = readout_sequences[barcode[0]] + readout_sequences[barcode[1]]
            right = readout_sequences[barcode[2]] + readout_sequences[barcode[3]]

            for j in oligo_database[i].keys():
                seq = str(oligo_database[i][j]["sequence"])
                seq = left + " T " + seq  # T is a spacer, to be outputted separately
                seq = seq + " T " + right
                oligo_database[i][j]["sequence"] = Seq(seq)
        return oligo_database

    def write_final_probes(
        self,
        oligo_database,
        readout_database,
        n_pseudocolors,
        seed,
    ):
        barcodes = create_seqfish_plus_barcodes(
            n_pseudocolors=n_pseudocolors,
            seed=seed,
            num_genes=len(oligo_database.database.keys()),
        )
        readout_sequences = readout_database.to_sequence_list()

        regions = oligo_database.database.keys()
        yaml_dict = {}

        for region in regions:
            barcode = barcodes[region]
            left = readout_sequences[barcode[0]] + readout_sequences[barcode[1]]
            right = readout_sequences[barcode[2]] + readout_sequences[barcode[3]]

            yaml_dict[region] = {}
            for i, oligo in enumerate(oligo_database[region].keys()):
                yaml_dict[region][f"{region}_oligo_{i+1}"] = {}
                seq = str(oligo_database[region][oligo]["sequence"])
                yaml_dict[region][f"{region}_oligo_{i+1}"]["id"] = oligo
                yaml_dict[region][f"{region}_oligo_{i+1}"]["region"] = region
                yaml_dict[region][f"{region}_oligo_{i+1}"][
                    "seqfish_plus_full_probe"
                ] = (left + "T" + seq + "T" + right)
                yaml_dict[region][f"{region}_oligo_{i+1}"].update(
                    oligo_database[region][oligo]
                )
                yaml_dict[region][f"{region}_oligo_{i+1}"]["start"] = ",".join(
                    f"{start}"
                    for start in yaml_dict[region][f"{region}_oligo_{i+1}"].pop("start")
                )
                yaml_dict[region][f"{region}_oligo_{i+1}"]["end"] = ",".join(
                    f"{end}"
                    for end in yaml_dict[region][f"{region}_oligo_{i+1}"].pop("end")
                )
                yaml_dict[region][f"{region}_oligo_{i+1}"][
                    "target_sequence"
                ] = yaml_dict[region][f"{region}_oligo_{i+1}"].pop("sequence")

                yaml_dict[region][f"{region}_oligo_{i+1}"].update(
                    {
                        f"readout_sequece_{i+1}": readout_sequences[barcode[i]]
                        for i in range(4)
                    }
                )
        # save
        probes_dir = os.path.join(self.dir_output, "final_seqfish_plus_sequences")
        Path(probes_dir).mkdir(parents=True, exist_ok=True)
        with open(os.path.join(probes_dir, "seqfish_plus_probes.yml"), "w") as outfile:
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
        with open(os.path.join(probes_dir, "merfish_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_order, outfile, default_flow_style=False, sort_keys=False)

        # Create readout probe file
        yaml_readout = {}
        yaml_readout["Bit"] = "Readout Probe"
        for i in range(num_bits):
            yaml_readout[str(i + 1)] = readout_probes[i] + "/3Cy5Sp"
        with open(
            os.path.join(probes_dir, "merfish_readout_probes.yml"), "w"
        ) as outfile:
            yaml.dump(yaml_readout, outfile, default_flow_style=False, sort_keys=False)

        # Create codebook file
        yaml_codebook = {}
        for gene_idx, gene in enumerate(genes):
            yaml_codebook[gene] = self.code[gene_idx]

        for i in range(self.n_blanks):
            yaml_codebook[f"blank_barcode_{i + 1}"] = self.code[n_genes + i]
        with open(os.path.join(probes_dir, "merfish_codebook.yml"), "w") as outfile:
            yaml.dump(yaml_codebook, outfile, default_flow_style=False, sort_keys=False)
