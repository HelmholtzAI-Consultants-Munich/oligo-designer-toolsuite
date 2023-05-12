############################################
# imports
############################################

from abc import ABC, abstractmethod

############################################
# Oligo Seed Gregion Creation Classes
############################################

class SeedRegionCreationBase(ABC):
    """
    Base class for all the seed region creation classes. The ``_create_seed_region`` method contains the functionalities for the specific seed region and the ``apply`` method
    is the method called by the user to generate the seed regions for an oligos_DB. Both need to be implemented for any new ``SeedRegionCreate`` class.
    """

    def __init__(self):
        pass

    @abstractmethod
    def apply(self, database: dict):
        """
        Generates the relative coordinates of the seed region for each oligo and stores them in the ``database`` database.
        To each oligo in the oligo DB are added the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param database: database containing all the oligo sequences and their features
        :type database: dict
        :return: database with added start and end position of the seed region for each oligo
        :rtype: dict
        """

    def _update_database(self, database: dict):
        """To each oligo in the oligo DB add the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param database: database containing all the oligo sequences and their features
        :type database: dict
        :return: database with added start and end position of the seed region for each oligo
        :rtype: dict
        """
        for region in database.keys():
            for oligo_id in database[region].keys():
                start, end = self._create_seed_region(database[region][oligo_id])
                database[region][oligo_id]["seed_region_start"] = start
                database[region][oligo_id]["seed_region_end"] = end
        return database

    @abstractmethod
    def _create_seed_region(self, database: dict):
        """
        This method actually defines the start and end point of the seed region according to the given conditions.

        :param database: Dictionalry containing all the information of the oligo
        :type database: dict
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

    def apply(self, database: dict):
        """
        Generates the relative coordinates of the seed region for each oligo and stores them in the ``database`` database.
        To each oligo in the oligo DB are added the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param database: database containing all the oligo sequences and their features
        :type database: dict
        :return: database with added start and end position of the seed region for each oligo
        :rtype: dict
        """
        database = super()._update_database(database)
        return database

    def _create_seed_region(self, database):
        """
        This method actually defines the start and end point of the seed region according to the given conditions.
        To each oligo in the oligo DB the relative coordinates of the start and the end of the seed region are given as arguments of the class.

        :param database: Dictionalry containing all the information of the oligo
        :type database: dict
        :return: Start and end position of the seed region of the given oligo
        :rtype: int, int
        """
        start = int(self.start)
        end = int(min(self.end, database["length"]))
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

    def apply(self, database: dict):
        """
        Generates the relative coordinates of the seed region for each oligo and stores them in the ``database`` database.
        To each oligo in the oligo DB are added the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param database: database containing all the oligo sequences and their features
        :type database: dict
        :return: database with added start and end position of the seed region for each oligo
        :rtype: dict
        """
        database = super()._update_database(database)
        return database

    def _create_seed_region(self, database):
        """
        This method actually defines the start and end point of the seed region according to the given conditions.

        :param database: Dictionalry containing all the information of the oligo
        :type database: dict
        :return: Start and end position of the seed region of the given oligo
        :rtype: int, int
        """

        length = database["length"]
        start = int(round(self.start * length))
        end = int(round(self.end * length))
        return start, end


class LigationRegionCreation(SeedRegionCreationBase):
    """The seed region is extended by ``ligation_region_size`` number of bases in both directions starting from the ligation site.
    To keep the dimension of the region 2*``ligation_region_size`` the left size is expanded by 1 base less.

    **Remark:** it is required to have a ligation site argument for each oligo in the datset, which can be computed with the ``PadlockArms`` filter in the prefilterig step

    :param ligation_region_size: number of basis the region has to expand starting from the ligation site.
    :type ligation_region_size: int
    """

    def __init__(self, ligation_region_size: int):
        super().__init__()
        self.ligation_region_size = ligation_region_size

    def apply(self, database: dict):
        """
        Generates the relative coordinates of the seed region for each oligo and stores them in the ``database`` database.
        To each oligo in the oligo DB are added the  ``start`` and ``end`` values with key "seed_region_start" and "seed_region_end" which indicate the relative coordinates of the start and the end (included) of the seed region.

        :param database: database containing all the oligo sequences and their features
        :type database: dict
        :return: database with added start and end position of the seed region for each oligo
        :rtype: dict
        """
        # check if the feature ligation site has been already computed
        genes = list(database.keys())
        i = 0
        while database[genes[i]] == {}:  # make sure we are not in an empty region
            i += 1
        oligo_id = list(database[genes[i]].keys())[0]
        if "ligation_site" not in database[genes[i]][oligo_id].keys():
            raise KeyError("The ligation site has not been computed!")

        database = super()._update_database(database)
        return database

    def _create_seed_region(self, database):
        """
        This method actually defines the start and end point of the seed region according to the given conditions.

        :param database: Dictionalry containing all the information of the oligo
        :type database: dict
        :return: Start and end position of the seed region of the given oligo
        :rtype: int, int
        """

        ligation_site = database["ligation_site"]
        start = int(max(ligation_site - (self.ligation_region_size - 1), 0))
        end = int(
            min(
                ligation_site + self.ligation_region_size,
                database["length"],
            )
        )
        return start, end
