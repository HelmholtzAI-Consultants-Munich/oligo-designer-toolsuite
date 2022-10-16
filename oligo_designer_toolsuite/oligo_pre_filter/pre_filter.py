class PreFilter:
    def __init__(self, filters) -> None:
        """
        :param filters: list of filters classes already initialized
        :type filters: list of classes
        """

        self.filters = filters

    def apply(self, DB):
        """Filter the database of probes based on teh given filters

        :param DB: database class containig the probes and their features
        :type DB: CustomDB class
        :return: database classs cointainig filtered oligos
        :rtype: CustomDB class
        """

        oligos_DB = DB.oligos_DB
        gene_ids = list(oligos_DB.keys())
        for gene_id in gene_ids:
            probes_id = list(oligos_DB[gene_id].keys())
            for probe_id in probes_id:
                fulfills, additional_features = self.filter_sequence(
                    oligos_DB[gene_id][probe_id]["probe_sequence"]
                )
                if fulfills:
                    oligos_DB[gene_id][probe_id].update(additional_features)
                else:
                    del oligos_DB[gene_id][probe_id]
        DB.oligos_DB = oligos_DB
        return DB

    def filter_sequence(self, sequence):
        """Applies the used-defined filters and returns the result and the additional computed features

        :param sequence: sequence to check
        :type sequence: Bio.Seq
        :return: if the filters are fulfilled and the additional features computed
        :rtype: bool, dict
        """

        fulfills = True
        additional_features = {}
        for filter in self.filters:
            out, feat = filter.apply(sequence)
            if not out:  # stop at the first false we obtain
                return False, {}
            additional_features.update(feat)
        return fulfills, additional_features
