import numpy as np
from Bio.Seq import Seq


class SeqfishProbesCreator:
    """
    This class is used to assemble probes using primary probes, readout probes and barcodes, that were designed for each gene.
    """

    def __init__(
        self,
    ):
        pass

    def create_probes(self, oligos_DB, readout_probes, barcodes):
        """ "
        Function for assembling probes.
        :param oligos_DB: database of oligos; dictionary, were each key corresponds to one gene
        :type oligos_DB: OligoDatabase
        :param readout_probes: list of readout probes
        :type readout_probes: list of Seq
        :param barcodes: barcodes for each gene from the database (each barcode consists of 4 digits)
        :type barcodes: dict (str : list of 4 int)
        :return: oligo_database, where "sequence" are assembled sequences
        :rtype: OligoDatabase
        """
        for i in oligos_DB.keys():
            barcode = barcodes[i]
            s1 = readout_probes[barcode[0]]
            s2 = readout_probes[barcode[1]]
            s3 = readout_probes[barcode[2]]
            s4 = readout_probes[barcode[3]]
            left = str(s1) + str(s2)
            right = str(s3) + str(s4)
            for j in oligos_DB[i].keys():
                seq = str(oligos_DB[i][j]["sequence"])
                seq = left + seq
                seq = seq + right
                oligos_DB[i][j]["sequence"] = Seq(seq)
        return oligos_DB


# TODO: Use database isntead of dict, documentation
class SeqFishReadoutProbes(ReadoutProbesBase):
    """ "
    Class to create readout probes in SeqFISH+ experiment
    :param length: length of readout probe
    :type length: int
    :param number_probes: number of probes, that should be created (correspond to the num of pseudocolors in the experiment)
    :type number_probes: int
    :param GC_min: minimum GC-content of the probe
    :type GC_min: int
    :param GC_max: maximum GC-content of the probe
    :type GC_max: int
    :param number_consecutive: min num of consecutive nucleotides, that are not allowed in the probe
    :type number_consecutive: int
    :param blast_filter: Blast filter (typically this filter was already created during specificity filtering)
    :type blast_filter: BlastnFilter
    :param reference_DB: reference database, that was created for Blast specificity filtering
    :type reference_DB: ReferenceDatabase
    """

    def __init__(
        self,
        length: int,
        number_probes: int,
        GC_min: int,
        GC_max: int,
        number_consecutive: int,
        random_seed: int,
        dir_output: str,
        blast_filter,
        reference_DB,
        sequence_alphabet: list[str] = ["A", "C", "G", "T"],
    ):
        self.length = length
        self.num_probes = number_probes
        self.property_filters = [
            GCContent(GC_content_min=GC_min, GC_content_max=GC_max),
            ConsecutiveRepeats(num_consecutive=number_consecutive),
        ]
        self.specificity_filters = [blast_filter]

        self.sequence_alphabet = sequence_alphabet
        self.ref = os.path.basename(reference_DB.file_fasta)
        self.seed = random_seed
        self.dir_output = dir_output
        self.readout_database = None

    def create_readout_probes(self, property_filter, specificity_filter):
        """
        Function, that creates readout probes
        :return: list of readout sequences, that fulfil experiment contraints
        :rtype: list of Seq()
        """

        self.readout_database = OligoDatabase(
            file_fasta=None, dir_output=self.dir_output
        )

        self.readout_database.create_random_database(
            self.length * 100, self.num_probes, sequence_alphabet=self.sequence_alphabet
        )

        # property_filter = PropertyFilter(self.property_filters)
        self.readout_database = property_filter.apply(self.readout_database)

        # specificity_filter = SpecificityFilter(self.specificity_filters)
        self.readout_database = specificity_filter.apply(
            self.readout_database, self.ref
        )
