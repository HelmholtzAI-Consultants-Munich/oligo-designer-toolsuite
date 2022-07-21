import logging
import os
import warnings
from pathlib import Path

import pandas as pd


class ProbeFilterBase:
    """This is the base class for all filter classes

    :param no_batches: How many sub-processes to use
    :type no_batches: int
    :param ligation_region: If ligation region present
    :type ligation_region: bool"""

    def __init__(
        self,
        number_batches,
        ligation_region,
        dir_output,
        file_probe_info,
        genes,
        dir_annotations=None,
        number_subbatches=1,
        max_genes_in_batch=300,
    ):
        self.number_batches = (number_batches,)
        self.ligation_region = ligation_region
        self.dir_output = dir_output
        self.number_subbatches = number_subbatches  # if number of genes in batch > max_genes_in_batch then split batch into multiple subbatches
        self.max_genes_in_batch = max_genes_in_batch  # if more than 300 genes in one batch split into subbatches to reduce required memory for loading blast results
        self.file_probe_info = file_probe_info
        self.genes = genes
        self.removed_genes = genes

        # set logger
        self.logging = logging.getLogger("filter_probes")

        # set directory
        if dir_annotations == None:
            self.dir_annotations = os.path.join(self.dir_output, "annotations")
            Path(self.dir_annotations).mkdir(parents=True, exist_ok=True)
        else:
            self.dir_annotations = dir_annotations

    def create_batches(self, number_batches):

        probeinfo_tsv = pd.read_csv(self.file_probe_info, sep="\t", header=0)

        batch_size = int(len(self.genes) / number_batches) + (
            len(self.genes) % number_batches > 0
        )

        for batch_id in range(number_batches):

            file_probe_info_batch = os.path.join(
                self.dir_annotations, "probes_info_batch{}.txt".format(batch_id)
            )
            # file_region_fasta_batch = os.path.join(self.dir_annotations, "region_batch{}.fna".format(batch_id))
            file_probe_sequence_batch = os.path.join(
                self.dir_annotations, "probes_sequence_batch{}.txt".format(batch_id)
            )

            # batch probe info
            genes_batch = self.genes[
                (batch_size * batch_id) : (
                    min(batch_size * (batch_id + 1), len(self.genes) + 1)
                )
            ]

            probeinfo = probeinfo_tsv.loc[
                probeinfo_tsv["gene_id"].isin(genes_batch)
            ].copy()
            probeinfo.reset_index(inplace=True, drop=True)
            probeinfo.to_csv(file_probe_info_batch, sep="\t", header=True, index=False)

            sequences = probeinfo[probeinfo.columns[1]]
            with open(file_probe_sequence_batch, "w") as handle_out:
                for seq in sequences:
                    handle_out.write(seq + "\n")

    def apply_base(self):
        """Apply filter to file_transcriptome_fasta and save output in dir_output"""

        warnings.warn(f"No apply function  for {type(self).__name__}")
