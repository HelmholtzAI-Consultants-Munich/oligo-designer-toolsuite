import os

import pandas as pd


def _write_output(self, probes_info, gene_id, probes_wo_match):
    """Write results of probe design pipeline to file and create one file with suitable probes per gene.
    :param probes_info: Dataframe with probe information, filtered based on sequence properties.
    :type probes_info: pandas.DataFrame
    :param gene_id: Gene ID of processed gene.
    :type gene_id: string
    :param probes_wo_match: List of suitable probes that don't have matches in the transcriptome.
    :type probes_wo_match: list
    """
    file_output = os.path.join(self.dir_probes, "probes_{}.txt".format(gene_id))
    valid_probes = probes_info[probes_info["probe_id"].isin(probes_wo_match)]
    valid_probes.to_csv(file_output, sep="\t", index=False)


def _load_probes_info(self, batch_id):
    """Load filtered probe information from tsv file
    :param batch_id: Batch ID.
    :type batch_id: int
    :return: Dataframe with probe information, filtered based on sequence properties.
    :rtype: pandas.DataFrame
    """
    file_probe_info_batch = os.path.join(
        self.dir_annotations, "probes_info_batch{}.txt".format(batch_id)
    )
    probes_info = pd.read_csv(
        file_probe_info_batch,
        sep="\t",
        dtype={
            "probe_id": str,
            "gene_id": str,
            "probe_sequence": str,
            "transcript_id": str,
            "exon_id": str,
            "chromosome": str,
            "start": str,
            "end": str,
            "strand": str,
            "GC_content": float,
            "melting_temperature": float,
            "melt_temp_arm1": float,
            "melt_temp_arm2": float,
            "dif_melt_temp_arms": float,
            "ligation_site": int,
        },
    )

    return probes_info


def _write_removed_genes(self):
    """Write list of genes for which not enough probes could be designed for."""

    # create file where removed genes are saved
    _, _, probe_files = next(os.walk(self.dir_probes))
    for probe_file in probe_files:
        gene_id = probe_file[len("probes_") : -len(".txt")]
        if gene_id in self.removed_genes:
            self.removed_genes.remove(gene_id)

    with open(self.file_removed_genes, "w") as output:
        for gene_id in self.removed_genes:
            output.write(f"{gene_id}\t0\n")


def mismatch_in_ligation(row):
    # Checks if mismatches are in ligation region, and if True, counts and stores how many mismatches are in ligation region
    mismatches_in_region = []
    for mismatch in row["positions"]:
        if int(mismatch) >= int(row["ligation_region_start"]) and int(mismatch) <= int(
            row["ligation_region_end"]
        ):
            mismatches_in_region.append(True)
    if len(mismatches_in_region) > 0:
        row["mismatch_in_ligation"] = mismatches_in_region[0]
        row["num_mismatches_in_ligation"] = len(mismatches_in_region)
    return row
