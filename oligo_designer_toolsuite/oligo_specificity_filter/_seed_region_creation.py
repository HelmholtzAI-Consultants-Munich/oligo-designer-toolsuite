from abc import ABC, abstractmethod


class SeedRegionCreationBase(ABC):
    """
    Base class for all the seed region creation classes. The ``_create_seed_region`` method contains the functionalities for the specific seed region and the ``apply`` method
    is the method called by the user to generate the seed regions for an oligos_DB. Both need to be implemented for any new ``SeedRegionCreate`` class.
    """

    def __init__(self):
        pass

    @abstractmethod
    def apply(self, oligo_DB: dict):
        """
        Generates the relative coordinates of the seed region for each oligo and stores them in the ``oligo_DB`` database.
        To each oligo in the oligo DB are added the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param oligo_DB: database containing all the oligo sequences and their features
        :type oligo_DB: dict
        :return: Oligo_DB with added start and end position of the seed region for each oligo
        :rtype: dict
        """

    def _update_oligo_DB(self, oligo_DB: dict):
        """To each oligo in the oligo DB add the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param oligo_DB: database containing all the oligo sequences and their features
        :type oligo_DB: dict
        :return: Oligo_DB with added start and end position of the seed region for each oligo
        :rtype: dict
        """
        for gene in oligo_DB.keys():
            for oligo_id in oligo_DB[gene].keys():
                start, end = self._create_seed_region(oligo_DB[gene][oligo_id])
                oligo_DB[gene][oligo_id]["seed_region_start"] = start
                oligo_DB[gene][oligo_id]["seed_region_end"] = end
        return oligo_DB

    @abstractmethod
    def _create_seed_region(self, oligo_DB: dict):
        """
        This method actually defines the start and end point of the seed region according to the given conditions.

        :param oligo_DB: Dictionalry containing all the information of the oligo
        :type oligo_DB: dict
        :return: Start and end position of the seed region of the given oligo
        :rtype: int, int
        """


class SeedRegionCreationStandard(SeedRegionCreationBase):
    """Standard seed region creation where the relative coordinates of the seed region are directly passed as arguments. If the end position is bigger than the length of the oligo the seed region ends where the oligo ends.

    :param start: relative coordinate of the start of the seed region (counting starts from 0)
    :type start: int
    :param end: relative coordinate of the end of the seed region (counting starts from 0)
    :type end: int
    """

    def __init__(self, start: int, end: int):
        super().__init__()
        self.start = start
        self.end = end

    def apply(self, oligo_DB: dict):
        """
        Generates the relative coordinates of the seed region for each oligo and stores them in the ``oligo_DB`` database.
        To each oligo in the oligo DB are added the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param oligo_DB: database containing all the oligo sequences and their features
        :type oligo_DB: dict
        :return: Oligo_DB with added start and end position of the seed region for each oligo
        :rtype: dict
        """
        oligo_DB = super()._update_oligo_DB(oligo_DB)
        return oligo_DB

    def _create_seed_region(self, oligo_DB):
        """
        This method actually defines the start and end point of the seed region according to the given conditions.
        To each oligo in the oligo DB the relative coordinates of the start and the end of the seed region are given as arguments of the class.

        :param oligo_DB: Dictionalry containing all the information of the oligo
        :type oligo_DB: dict
        :return: Start and end position of the seed region of the given oligo
        :rtype: int, int
        """
        start = int(self.start)
        end = int(min(self.end, oligo_DB["length"]))
        return start, end


class SeedRegionCreationPercentage(SeedRegionCreationBase):
    """The seed region relative coordinates are give as a percentage with respect to the oligo length in [0,1].
    For example if is requested a seed region betwee 0.3 and 0.7,  then the start and end positions woudl be:
    - (3,7) for a oligo of 10 basis
    - (6,14) for a oligo of 20 bases

    :param start: relative coordinate of the start of the seed region as a pecentage wrt to oligo lenght
    :type start: float
    :param end: relative coordinate of the end of the seed region as a pecentage wrt to oligo lenght
    :type end: float
    """

    def __init__(self, start: float, end: float):
        super().__init__()
        if start < 0 or start > 1:
            raise ValueError("Start position must be in the interval [0,1]!")
        self.start = start
        if end < 0 or end > 1:
            raise ValueError("End position must be in the interval [0,1]!")
        self.end = end

    def apply(self, oligo_DB: dict):
        """
        Generates the relative coordinates of the seed region for each oligo and stores them in the ``oligo_DB`` database.
        To each oligo in the oligo DB are added the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param oligo_DB: database containing all the oligo sequences and their features
        :type oligo_DB: dict
        :return: Oligo_DB with added start and end position of the seed region for each oligo
        :rtype: dict
        """
        oligo_DB = super()._update_oligo_DB(oligo_DB)
        return oligo_DB

    def _create_seed_region(self, oligo_DB):
        """
        This method actually defines the start and end point of the seed region according to the given conditions.

        :param oligo_DB: Dictionalry containing all the information of the oligo
        :type oligo_DB: dict
        :return: Start and end position of the seed region of the given oligo
        :rtype: int, int
        """

        length = oligo_DB["length"]
        start = int(round(self.start * length))
        end = int(round(self.end * length))
        return start, end


class LigationRegionCreation(SeedRegionCreationBase):
    """Seed region used in the Padlock oligo designer pipeline, the region extends by ``ligation_region_size`` number of bases in both directions starting from the ligation site.
    To keep the dimention of the region 2*``ligation_region_size`` the left size is expanded by 1 base less.

    **Remark:** it is required to have a ligation site argument for each oligo in the datset, which can be computed with the ``PadlockArms`` filter in the prefilterig step

    :param ligation_region_size: number of basis the region has to expand starting from the ligation site.
    :type ligation_region_size: int
    """

    def __init__(self, ligation_region_size: int):
        super().__init__()
        self.ligation_region_size = ligation_region_size

    def apply(self, oligo_DB: dict):
        """
        Generates the relative coordinates of the seed region for each oligo and stores them in the ``oligo_DB`` database.
        To each oligo in the oligo DB are added the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param oligo_DB: database containing all the oligo sequences and their features
        :type oligo_DB: dict
        :return: Oligo_DB with added start and end position of the seed region for each oligo
        :rtype: dict
        """
        # check if the feature ligation site has been already computed
        genes = list(oligo_DB.keys())
        i = 0
        while oligo_DB[genes[i]] == {}:  # make sure we are not in an empty gene
            i += 1
        oligo_id = list(oligo_DB[genes[i]].keys())[0]
        if "ligation_site" not in oligo_DB[genes[i]][oligo_id].keys():
            raise KeyError("The ligation site has not been computed!")

        oligo_DB = super()._update_oligo_DB(oligo_DB)
        return oligo_DB

    def _create_seed_region(self, oligo_DB):
        """
        This method actually defines the start and end point of the seed region according to the given conditions.

        :param oligo_DB: Dictionalry containing all the information of the oligo
        :type oligo_DB: dict
        :return: Start and end position of the seed region of the given oligo
        :rtype: int, int
        """

        ligation_site = oligo_DB["ligation_site"]
        start = int(max(ligation_site - (self.ligation_region_size - 1), 0))
        end = int(
            min(
                ligation_site + self.ligation_region_size,
                oligo_DB["length"],
            )
        )
        return start, end
