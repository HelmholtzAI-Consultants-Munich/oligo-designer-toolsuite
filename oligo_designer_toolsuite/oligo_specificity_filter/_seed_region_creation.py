from abc import ABC, abstractmethod


class SeedRegionCreationBase(ABC):
    """Base class tfor all the seed region creation classes"""

    def __init__(self):
        pass

    @abstractmethod
    def apply(self, oligo_DB):
        """To each probe in the oligo DB add the  `start` and `end` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param oligo_DB: database containin all teh oligo sequences and their features
        :type oligo_DB: dict
        """
        return oligo_DB


class SeedRegionCreationStandard(SeedRegionCreationBase):
    """Standard seed region creation where the relative coordinates of the seed region are directly passed as arguments. If the end positionis bigger than teh length of the probe the seed region ends where teh probe ends.

    :param start: relative coordinate of the start of the seed region (counting starts from 0)
    :type start: int
    :param end: relative coordinate of the end of the seed region (counting starts from 0)
    :type end: int
    """

    def __init__(self, start, end):
        super().__init__()
        self.start = start
        self.end = end

    def apply(self, oligo_DB):
        """To each probe in the oligo DB add a tuple `(start,end)` with key "seed_region" which indicate the relative coordinates of the start and the end of the seed region.

        :param oligo_DB: database containin all teh oligo sequences and their features
        :type oligo_DB: dict
        """
        for gene in oligo_DB.keys():
            for probe_id in oligo_DB[gene].keys():
                oligo_DB[gene][probe_id]["seed_region_start"] = int(self.start)
                oligo_DB[gene][probe_id]["seed_region_end"] = int(
                    min(self.end, oligo_DB[gene][probe_id]["length"])
                )
        return oligo_DB


class SeedRegionCreationPercentage(SeedRegionCreationBase):
    """The seed region relative coordinates are give as a percentage with respect to the probe length in [0,1].
        For example if is requested a seed region betwee 0.3 and 0.7,  then teh start and end positions woudl be:
        - (3,7) for a probe of 10 basis
        - (6,14) for a probe of 20 bases
    r
        :param start: relative coordinate of the start of the seed region as a pecentage wrt to probe lenght
        :type start: float
        :param end: relative coordinate of the end of the seed region as a pecentage wrt to probe lenght
        :type end: float
    """

    def __init__(self, start, end):
        super().__init__()
        if start < 0 or start > 1:
            raise ValueError("Etart position must be in the interval [0,1]!")
        self.start = start
        if end < 0 or end > 1:
            raise ValueError("End position must be in the interval [0,1]!")
        self.end = end

    def apply(self, oligo_DB):
        """To each probe in the oligo DB add a tuple `(start,end)` with key "seed_region" which indicate the relative coordinates of the start and the end of the seed region.

        :param oligo_DB: database containin all teh oligo sequences and their features
        :type oligo_DB: dict
        """
        for gene in oligo_DB.keys():
            for probe_id in oligo_DB[gene].keys():
                length = oligo_DB[gene][probe_id]["length"]
                oligo_DB[gene][probe_id]["seed_region_start"] = int(
                    round(self.start * length)
                )
                oligo_DB[gene][probe_id]["seed_region_end"] = int(
                    round(self.end * length)
                )
        return oligo_DB


class LigationRegionCreation(SeedRegionCreationBase):
    """Seed regionused in the Paddlock probe designer pipeline, the region extends by `ligation_region_size` number of bases in both directions starting from the ligation site.
    To keep the dimention of the region 2*`ligation_region_size` the left size is expanded by 1 base less.

    **Remark:** it is required to have a ligation site argument for each probe in the datset, which can be computed with the `PadlockArms` filter in the prefilterig step

    :param ligation_region_size: number of basis the regioin has to expand starting from the ligation site.
    :type ligation_region_size: int
    """

    def __init__(self, ligation_region_size):
        super().__init__()
        self.ligation_region_size = ligation_region_size

    def apply(self, oligo_DB):
        """To each probe in the oligo DB add a tuple `(start,end)` with key "seed_region" which indicate the relative coordinates of the start and the end of the seed region.

        :param oligo_DB: database containin all teh oligo sequences and their features
        :type oligo_DB: dict
        """
        # check if the feature ligation site has been already computed
        genes = list(oligo_DB.keys())
        i = 0
        while oligo_DB[genes[i]] == {}:  # make sure we are not in an empty gene
            i += 1
        probe_id = list(oligo_DB[genes[i]].keys())[0]
        if "ligation_site" not in oligo_DB[genes[i]][probe_id].keys():
            raise KeyError("The ligation site has not been computed!")

        for gene in oligo_DB.keys():
            for probe_id in oligo_DB[gene].keys():
                ligation_site = oligo_DB[gene][probe_id]["ligation_site"]
                oligo_DB[gene][probe_id]["seed_region_start"] = int(
                    max(ligation_site - (self.ligation_region_size - 1), 0)
                )
                oligo_DB[gene][probe_id]["seed_region_end"] = int(
                    min(
                        ligation_site + self.ligation_region_size,
                        oligo_DB[gene][probe_id]["length"],
                    )
                )
        return oligo_DB
