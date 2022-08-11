def _write_output(probes_info, gene_id, probes_wo_match):
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
