############################################
# imports
############################################

import warnings
from typing import List, Union

from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import Seq, gc_fraction
from seqfold import dg

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.utils import check_if_key_exists, check_if_list

############################################
# Attrubite Calculation Class
############################################


class OligoAttributes:
    """
    Class for calculating and managing attributes related to oligonucleotides.

    This class is designed to work with an OligoDatabase instance to compute oligo attributes such as oligonucleotide length.
    """

    def __init__(self):
        """Constructor for the OligoAttributes class."""

    @staticmethod
    def _calc_oligo_length(sequence: str):
        """Calculate the length of a given oligonucleotide sequence.

        :param sequence: The oligonucleotide sequence to calculate the length for.
        :type sequence: str
        :return: The length of the sequence.
        :rtype: int
        """
        length = len(sequence)

        return length

    def calculate_oligo_length(self, oligo_database: OligoDatabase, region_ids: Union[str, List[str]] = None):
        """Calculate the length for each oligonucleotide in the database.

        :param oligo_database: The database containing oligonucleotide sequences and attributes.
        :type oligo_database: OligoDatabase
        :param region_ids: The region IDs for which to calculate the oligo length.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                # oligo and target have always same length
                length = self._calc_oligo_length(oligo_attributes["oligo"])
                oligo_attributes["length"] = length

        return oligo_database

    @staticmethod
    def _calc_num_targeted_transcripts(transcript_ids: list):
        """Calculate the number of unique transcripts targeted by a given oligonucleotide.

        :param transcript_ids: List of transcript IDs targeted by the oligonucleotide.
        :type transcript_ids: list
        :return: The number of unique targeted transcripts.
        :rtype: int
        """
        num_targeted_transcripts = len(
            set(
                item
                for sublist in (transcript_ids if isinstance(transcript_ids[0], list) else [transcript_ids])
                for item in sublist
            )
        )

        return num_targeted_transcripts

    def calculate_num_targeted_transcripts(
        self, oligo_database: OligoDatabase, region_ids: Union[str, List[str]] = None
    ):
        """Calculate the number of targeted transcripts for each oligonucleotide in the database.

        If the necessary information for number of targeted transcripts calculation is not available,
        the 'num_targeted_transcripts' attribute is set to None.

        :param oligo_database: The database containing oligonucleotide sequences and attributes.
        :type oligo_database: OligoDatabase
        :param region_ids: The region IDs for which to calculate the number of targeted transcripts.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                if ("transcript_id" in oligo_attributes) and (oligo_attributes["transcript_id"]):
                    num_targeted_transcripts = self._calc_num_targeted_transcripts(
                        oligo_attributes["transcript_id"]
                    )
                else:
                    num_targeted_transcripts = None
                oligo_attributes["num_targeted_transcripts"] = num_targeted_transcripts

        return oligo_database

    @staticmethod
    def _calc_isoform_consensus(transcript_ids: list, number_transcripts: list):
        """Calculate the isoform consensus percentage of a given oligonucleotide.

        This function calculates the isoform consensus based on the provided transcript information.
        It computes the percentage of unique transcript IDs of the oligo over the total number of transcripts
        associated with the oligo region. The maximum value for the isoform consensus is 100%, which means
        that the oligo is present in all isoforms (transcripts) of the region.

        :param transcript_id: Transcript IDs targeted by the oligonucleotide.
        :type transcript_id: list
        :param number_transcripts: Total number of transcripts for the oligo region.
        :type number_transcripts: list
        :return: The isoform consensus percentage.
        :rtype: float
        """
        # number transcripts is the number of transcripts of a genomic region
        # hence, all values have to be the same for each transcript coming from the same oligo
        # since only oligos from the same genomic region are merged into one entry
        number_transcripts = int([item for sublist in number_transcripts for item in sublist][0])
        num_targeted_transcripts = len(
            set(
                item
                for sublist in (transcript_ids if isinstance(transcript_ids[0], list) else [transcript_ids])
                for item in sublist
            )
        )
        isoform_consensus = num_targeted_transcripts / number_transcripts * 100

        return isoform_consensus

    def calculate_isoform_consensus(
        self, oligo_database: OligoDatabase, region_ids: Union[str, List[str]] = None
    ):
        """Calculate the isoform consensus percentage for each oligonucleotide in the database.

        If the necessary information for isoform consensus calculation is not available,
        the 'isoform_consensus' attribute is set to None.

        :param oligo_database: The database containing oligonucleotide sequences and attributes.
        :type oligo_database: OligoDatabase
        :param region_ids: The region IDs for which to calculate the isoform consensus.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                if (
                    ("transcript_id" in oligo_attributes)
                    and (oligo_attributes["transcript_id"])
                    and ("number_transcripts" in oligo_attributes)
                    and (oligo_attributes["number_transcripts"])
                ):
                    isoform_consensus = self._calc_isoform_consensus(
                        oligo_attributes["transcript_id"], oligo_attributes["number_transcripts"]
                    )
                else:
                    isoform_consensus = None
                oligo_attributes["isoform_consensus"] = isoform_consensus

        return oligo_database

    @staticmethod
    def _calc_seedregion(sequence: str, start: Union[int, float], end: Union[int, float]):
        """Calculate the start and end positions of a seed region within a given oligonucleotide sequence.

        The seed region is calculated based on start and end parameters. The start and end can be specified as absolute
        positions (int) or as a percentage of the oligo's length (float).

        For example:
        start = 4
        end = 6
            will set the relative start and end positions wrt the oligo sequence of the seed region to 4 and 6, respectively.

        start = 0.4
        end = 0.6
            will set the relative start and end positions wrt the oligo sequence of the seed region to 4 and 6, respectively,
            only if the oligo length = 10.

        :param sequence: The sequence for which the seed region is calculated.
        :type sequence: str
        :param start: The start position of the seed region, as an index or a fraction of the sequence length.
        :type start: Union[int, float]
        :param end: The end position of the seed region, as an index or a fraction of the sequence length.
        :type end: Union[int, float]
        :return: The calculated start and end positions of the seed region.
        :rtype: tuple[int, int]
        :raises ValueError: If start and end types do not match, or if float values are out of the [0,1] range.
        """
        length = len(sequence)

        if isinstance(start, int) and isinstance(end, int):
            seedregion_start = max(0, start)
            seedregion_end = min(length, end)
        elif isinstance(start, float) and isinstance(end, float):
            if (not 0 <= start <= 1) or (not 0 <= end <= 1):
                raise ValueError("Start and end positions must be in the interval [0,1] for float type.")
            seedregion_start = int(round(start * length))
            seedregion_end = int(round(end * length))
        else:
            raise ValueError("Start and end parameters must be both integers or both floats.")

        return seedregion_start, seedregion_end

    def calculate_seedregion(
        self,
        oligo_database: OligoDatabase,
        start: Union[int, float],
        end: Union[int, float],
        region_ids: Union[str, List[str]] = None,
    ):
        """Calculate the seed region start and end positions for each oligonucleotide in the database.

        :param oligo_database: The database containing oligonucleotide sequences and attributes.
        :type oligo_database: OligoDatabase
        :param start: The start position of the seed region, either as an absolute position or a fraction of the sequence length.
        :type start: Union[int, float]
        :param end: The end position of the seed region, either as an absolute position or a fraction of the sequence length.
        :type end: Union[int, float]
        :param region_ids: The region IDs for which to calculate the seed region.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                # oligo and target have always same length
                seedregion_start, seedregion_end = self._calc_seedregion(
                    oligo_attributes["oligo"], start, end
                )
                oligo_attributes["seedregion_start"] = seedregion_start
                oligo_attributes["seedregion_end"] = seedregion_end

        return oligo_database

    @staticmethod
    def _calc_seedregion_ligationsite(sequence: str, ligation_site: int, seedregion_size: int):
        """Calculate the seed region around a specified ligation site within a given oligonucleotide sequence.

        The seed region is calculated based on a specified seed region size. The seed region is defined
        symmetrically around the ligation site, considering the provided size.

        :param sequence: The sequence for which the seed region is calculated.
        :type sequence: str
        :param ligation_site: The position of the ligation site within the sequence.
        :type ligation_site: int
        :param seedregion_size: The total size of the seed region to calculate around the ligation site.
        :type seedregion_size: int
        :return: The start and end positions of the calculated seed region.
        :rtype: tuple[int, int]
        """
        length = len(sequence)

        seedregion_start = int(max(0, ligation_site - (seedregion_size - 1)))
        seedregion_end = int(
            min(
                length,
                ligation_site + seedregion_size,
            )
        )
        return seedregion_start, seedregion_end

    def calculate_seedregion_ligationsite(
        self, oligo_database: OligoDatabase, seedregion_size: int, region_ids: Union[str, List[str]] = None
    ):
        """Calculate the seed region around a specified ligation site for each oligonucleotide in the database.

        If the necessary information for seed region calculation is not available,
        the 'seedregion_start' and 'seedregion_end' attributes are set to None.

        :param oligo_database: The database containing oligonucleotide sequences and attributes.
        :param seedregion_size: The size of the seed region to calculate around the ligation site.
        :type oligo_database: OligoDatabase
        :type seedregion_size: int
        :param region_ids: The region IDs for which to calculate the seed region.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        :raises KeyError: If the ligation site attribute is missing from any oligonucleotide in the database.
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            if not check_if_key_exists(database_region, "ligation_site"):
                warnings.warn(
                    f"The ligation_site attribute has not been computed for {region_id}! Setting to None!"
                )
            for oligo_id, oligo_attributes in database_region.items():
                if ("ligation_site" in oligo_attributes) and (oligo_attributes["ligation_site"]):
                    # oligo and target have always same length
                    seedregion_start, seedregion_end = self._calc_seedregion_ligationsite(
                        oligo_attributes["oligo"], oligo_attributes["ligation_site"], seedregion_size
                    )
                else:
                    seedregion_start = seedregion_end = None

                oligo_attributes["seedregion_start"] = seedregion_start
                oligo_attributes["seedregion_end"] = seedregion_end

        return oligo_database

    @staticmethod
    def _calc_GC_content(sequence: str):
        """Calculate the GC content of a given oligonucleotide sequence, expressed as a percentage.

        :param sequence: The DNA sequence for which the GC content is calculated.
        :type sequence: str
        :return: The GC content percentage of the sequence.
        :rtype: float
        """
        GC_content = round(gc_fraction(sequence) * 100, 2)
        return GC_content

    def calculate_GC_content(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ):
        """Calculate the GC content for each oligonucleotide in the database, dependent on the specified sequence type.

        :param oligo_database: The database containing oligonucleotide sequences and attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to consider for GC content calculation (e.g., 'oligo' or 'target').
        :type sequence_type: _TYPES_SEQ
        :param region_ids: The region IDs for which to calculate the GC content.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                GC_content = self._calc_GC_content(oligo_attributes[sequence_type])
                oligo_attributes["GC_content"] = GC_content

        return oligo_database

    @staticmethod
    def _calc_TmNN(
        sequence: str,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ):
        """Calculate the melting temperature of a given oligonucleotide using nearest-neighbor thermodynamics,
        with optional salt and chemical corrections.

        :param sequence: The DNA sequence for which Tm is calculated.
        :type sequence: str
        :param Tm_parameters: Parameters for the nearest-neighbor thermodynamic model to calculate Tm.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type Tm_parameters: dict
        :param Tm_salt_correction_parameters: Optional parameters for salt correction of Tm calculations.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type Tm_salt_correction_parameters: dict, optional
        :param Tm_chem_correction_parameters: Optional parameters for chemical correction of Tm calculations.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type Tm_chem_correction_parameters: dict, optional
        :return: The calculated melting temperature.
        :rtype: float
        """
        TmNN = mt.Tm_NN(sequence, **Tm_parameters)
        if Tm_salt_correction_parameters is not None:
            TmNN += mt.salt_correction(**Tm_salt_correction_parameters, seq=sequence)
        if Tm_chem_correction_parameters is not None:
            TmNN = mt.chem_correction(TmNN, **Tm_chem_correction_parameters)
        TmNN = round(TmNN, 2)
        return TmNN

    def calculate_TmNN(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
        region_ids: Union[str, List[str]] = None,
    ):
        """Calculate the melting temperature for each oligonucleotide in the database, dependent on the specified sequence type,
        using nearest-neighbor thermodynamics with optional salt and chemical corrections.

        :param oligo_database: The database containing oligonucleotide sequences and attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to consider for Tm calculation (e.g., 'oligo' or 'target').
        :type sequence_type: _TYPES_SEQ
        :param Tm_parameters: Parameters for the nearest-neighbor Tm calculation.
        :type Tm_parameters: dict
        :param Tm_salt_correction_parameters: Optional parameters for salt correction.
        :type Tm_salt_correction_parameters: dict, optional
        :param Tm_chem_correction_parameters: Optional parameters for chemical correction.
        :type Tm_chem_correction_parameters: dict, optional
        :param region_ids: The region IDs for which to calculate the Tm.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                TmNN = self._calc_TmNN(
                    oligo_attributes[sequence_type],
                    Tm_parameters,
                    Tm_salt_correction_parameters,
                    Tm_chem_correction_parameters,
                )
                oligo_attributes["TmNN"] = TmNN

        return oligo_database

    @staticmethod
    def _calc_length_complement(sequence1: str, sequence2: str):
        """Calculate the length of the longest complementary region between two sequences.

        This function compares two sequences and determines the length of the longest complementary region,
        considering all possible alignments. Sequence2 is complemented before comparison to simulate binding.

        :param sequence1: The first DNA sequence.
        :type sequence1: str
        :param sequence2: The second DNA sequence.
        :type sequence2: str
        :return: The length of the longest complementary region.
        :rtype: int
        """

        def _calculate_max_overlap(seq1, seq2):
            len_overlap_sub = 0
            len_overlap = 0
            for c1, c2 in zip(seq1, seq2):
                if c1 != c2:
                    len_overlap_sub = 0
                else:
                    len_overlap_sub += 1
                len_overlap = max(len_overlap, len_overlap_sub)
            return len_overlap

        # since we are comparing strings, we take the complement of sequence 2,
        # which should be the exact same sequence as sequence 1 if they bind
        sequence2 = Seq(sequence2).complement()

        # Initialize max_len_overlap with overlap without shift
        max_len_overlap = _calculate_max_overlap(sequence1, sequence2)
        max_shift = max(len(sequence1), len(sequence2)) - max_len_overlap

        # Check all possible shifts
        for shift in range(-max_shift, max_shift + 1):
            if shift < 0:
                # Shift sequence2 to the left
                shifted_seq2 = sequence2[-shift:]
                shifted_seq1 = sequence1[: len(shifted_seq2)]
            else:
                # Shift sequence2 to the right
                shifted_seq2 = sequence2[:-shift] if shift != 0 else sequence2
                shifted_seq1 = sequence1[shift:]

            len_comp = _calculate_max_overlap(shifted_seq1, shifted_seq2)
            max_len_overlap = max(max_len_overlap, len_comp)

        return max_len_overlap

    def calculate_length_selfcomplement(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ):
        """Calculate the length of the longest self-complementary sequence for each oligonucleotide in the database,
        dependent on the specified sequence type.

        :param oligo_database: Database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence to analyze (e.g., 'oligo' or 'target').
        :type sequence_type: _TYPES_SEQ
        :param region_ids: The region IDs for which to calculate the length of the longest self-complementary sequence.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """

        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                # we want to check if the reverse of our sequence is complementary to itself, e.g.
                # 5' - TAA CAA TAT ATA TTG TTA - 3' and it's reverse
                # 3' - ATT GTT ATA TAT AAC AAT - 5' are complementary to each other
                sequence = oligo_attributes[sequence_type]
                sequence_rev = sequence[::-1]
                len_overlap = self._calc_length_complement(sequence, sequence_rev)
                oligo_attributes["length_selfcomplement"] = len_overlap

        return oligo_database

    def calculate_length_complement(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        comparison_sequence: str,
        comparison_sequence_name: str,
        region_ids: Union[str, List[str]] = None,
    ):
        """Calculate the length of the longest complementary sequence between two oligonucleotides in the database,
        dependent on the specified sequence types.

        :param oligo_database: Database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :param sequence_type1: Type of sequence to analyze (e.g., 'oligo' or 'target') for the first oligonucleotide.
        :type sequence_type1: _TYPES_SEQ
        :param comparison_sequence: The second DNA sequence to analyze for complementary sequences.
        :type comparison_sequence: str
        :param comparison_sequence_name: Name of the second DNA sequence.
        :type comparison_sequence_name: str
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id, database_region in oligo_database.database.items():
            for oligo_id, oligo_attributes in database_region.items():
                len_overlap = self._calc_length_complement(
                    oligo_attributes[sequence_type], comparison_sequence
                )
                oligo_attributes["length_complement_" + comparison_sequence_name] = len_overlap

        return oligo_database

    @staticmethod
    def _calc_DG_secondary_structure(sequence: str, T: float):
        """Calculate the Gibbs free energy (ΔG) of the secondary structure formation at a given temperature (T)
        of a given oligonucleotide sequence.

        :param sequence: DNA sequence to analyze.
        :type sequence: str
        :param T: Temperature in degrees Celsius.
        :type T: float
        :return: ΔG of secondary structure formation.
        :rtype: float
        """
        DG_secondary_structure = dg(sequence, temp=T)
        return DG_secondary_structure

    def calculate_DG_secondary_structure(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        T: float,
        region_ids: Union[str, List[str]] = None,
    ):
        """Calculate the Gibbs free energy (ΔG) of the secondary structure formation at a given temperature (T)
        for each oligonucleotide in the database, dependent on the specified sequence type.

        :param oligo_database: Database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence (e.g., 'oligo' or 'target').
        :type sequence_type: _TYPES_SEQ
        :param T: Temperature in degrees Celsius for ΔG calculation.
        :type T: float
        :param region_ids: The region IDs for which to calculate the deltaG of secondary structure formation.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                DG_secondary_structure = self._calc_DG_secondary_structure(oligo_attributes[sequence_type], T)
                oligo_attributes["DG_secondary_structure"] = DG_secondary_structure

        return oligo_database

    def _calc_padlock_arms(
        self,
        sequence: str,
        arm_length_min: int,
        arm_Tm_dif_max: float,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ):
        """Calculate the optimal ligation site and melting temperatures for padlock probe arms of a given oligonucleotide sequence.

        The method iteratively adjusts the ligation site to find arm lengths and Tm values that meet the criteria
        of minimum arm length and Tm differences within the specified maximum. The process stops once suitable
        arms are found or if no configuration meets the criteria.

        :param sequence: DNA sequence to analyze.
        :type sequence: str
        :param arm_length_min: Minimum length for each arm of the padlock probe.
        :type arm_length_min: int
        :param arm_Tm_dif_max: Maximum allowable Tm difference between arms.
        :type arm_Tm_dif_max: float
        :param arm_Tm_min: Minimum allowable melting temperature for each arm.
        :type arm_Tm_min: float
        :param arm_Tm_max: Maximum allowable melting temperature for each arm.
        :type arm_Tm_max: float
        :param Tm_parameters: Parameters for melting temperature calculation.
        :type Tm_parameters: dict
        :param Tm_salt_correction_parameters: Optional parameters for salt correction in Tm calculation.
        :type Tm_salt_correction_parameters: dict, optional
        :param Tm_chem_correction_parameters: Optional parameters for chemical correction in Tm calculation.
        :type Tm_chem_correction_parameters: dict, optional
        :return: Melting temperatures for both arms and the ligation site, or None if suitable arms cannot be found.
        :rtype: tuple(float, float, int) or tuple(None, None, None)
        """
        len_sequence = len(sequence)
        ligation_site = len_sequence // 2

        arms_long_enough = (ligation_site >= arm_length_min) and (
            (len_sequence - ligation_site) >= arm_length_min
        )
        Tm_found = False
        sign_factor = 1  # switch between positive and negative shift
        shift = 1  # distance of ligation site shift

        while arms_long_enough and not Tm_found:
            Tm_arm1 = self._calc_TmNN(
                sequence[:ligation_site],
                Tm_parameters,
                Tm_salt_correction_parameters,
                Tm_chem_correction_parameters,
            )
            Tm_arm2 = self._calc_TmNN(
                sequence[ligation_site:],
                Tm_parameters,
                Tm_salt_correction_parameters,
                Tm_chem_correction_parameters,
            )
            Tm_dif = round(abs(Tm_arm2 - Tm_arm1), 2)
            Tm_found = (
                (Tm_dif <= arm_Tm_dif_max)
                and (arm_Tm_min <= Tm_arm1 <= arm_Tm_max)
                and (arm_Tm_min <= Tm_arm2 <= arm_Tm_max)
            )
            if not Tm_found:
                ligation_site += sign_factor * shift
                sign_factor *= -1
                shift += 1
                arms_long_enough = (ligation_site >= arm_length_min) and (
                    (len_sequence - ligation_site) >= arm_length_min
                )

        if Tm_found:
            return Tm_arm1, Tm_arm2, ligation_site
        else:
            return None, None, None

    def calculate_padlock_arms(
        self,
        oligo_database: OligoDatabase,
        arm_length_min: int,
        arm_Tm_dif_max: float,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
        region_ids: Union[str, List[str]] = None,
    ):
        """Calculate the optimal ligation site and melting temperatures for padlock probe arms for each
        oligonucleotide in the database, dependent on the specified sequence type.

        :param oligo_database: Database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :param arm_length_min: Minimum length for each arm of the padlock probe.
        :type arm_length_min: int
        :param arm_Tm_dif_max: Maximum allowable Tm difference between arms.
        :type arm_Tm_dif_max: float
        :param arm_Tm_min: Minimum allowable melting temperature for each arm.
        :type arm_Tm_min: float
        :param arm_Tm_max: Maximum allowable melting temperature for each arm.
        :type arm_Tm_max: float
        :param Tm_parameters: Parameters for melting temperature calculation.
        :type Tm_parameters: dict
        :param Tm_salt_correction_parameters: Optional parameters for salt correction in Tm calculation.
        :type Tm_salt_correction_parameters: dict, optional
        :param Tm_chem_correction_parameters: Optional parameters for chemical correction in Tm calculation.
        :type Tm_chem_correction_parameters: dict, optional
        :param region_ids: The region IDs for which to calculate the padlock arms.
        :type region_ids: Union[str, List[str]]
        :return: The database containing the new oligo attribute.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                arm1_Tm, arm2_Tm, ligation_site = self._calc_padlock_arms(
                    oligo_attributes["oligo"],
                    arm_length_min,
                    arm_Tm_dif_max,
                    arm_Tm_min,
                    arm_Tm_max,
                    Tm_parameters,
                    Tm_salt_correction_parameters,
                    Tm_chem_correction_parameters,
                )
                oligo_attributes["arm1_Tm"] = arm1_Tm
                oligo_attributes["arm2_Tm"] = arm2_Tm
                oligo_attributes["ligation_site"] = ligation_site

        return oligo_database

    def _calc_detect_oligo(
        self,
        sequence: str,
        ligation_site: int,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        min_thymines: int,
    ):
        """Calculate detection oligo sequences based on specified constraints.

        This function calculates potential detection oligo sequences around the ligation site, ensuring they meet
        the specified length and thymine content constraints.

        :param sequence: The oligo sequence to be evaluated.
        :type sequence: str
        :param ligation_site: The position of the ligation site in the sequence.
        :type ligation_site: int
        :param detect_oligo_length_min: The minimum length of the detection oligo.
        :type detect_oligo_length_min: int
        :param detect_oligo_length_max: The maximum length of the detection oligo.
        :type detect_oligo_length_max: int
        :param min_thymines: The minimum number of thymine bases required in the detection oligo.
        :type min_thymines: int
        :return: A tuple containing the even-length detection oligo, longer left detection oligo, and longer right detection oligo.
        :rtype: tuple
        """
        # constraint: a difference of max 1 nt for the sequences left and right of the ligation site is allowed
        # e.g. AAA|TTTT or AAAA|TTT hence, the detetcion oligo can only be as long as the shorter arm + 1 nt
        detect_oligo_length = 2 * min(ligation_site, len(sequence) - ligation_site) + 1

        # check if min and max constraints are fulfilled
        if (detect_oligo_length_min > detect_oligo_length) or (detect_oligo_length_max == 0):
            return None, None, None

        detect_oligo_even = None
        detect_oligo_long_left = None
        detect_oligo_long_right = None

        # Different scenarios
        # 1. If the max length constraint is smaller than the length of the oligo
        if detect_oligo_length_max < detect_oligo_length:
            detect_oligo_length_max_half = detect_oligo_length_max // 2
            detect_oligo_even = sequence[
                ligation_site - detect_oligo_length_max_half : ligation_site + detect_oligo_length_max_half
            ]
            # 1.2 if the maximal length is odd -> return three different oligos: even, longer left, longer right
            if detect_oligo_length_max % 2 == 1:
                detect_oligo_long_left = sequence[
                    ligation_site
                    - detect_oligo_length_max_half
                    - 1 : ligation_site
                    + detect_oligo_length_max_half
                ]
                detect_oligo_long_right = sequence[
                    ligation_site
                    - detect_oligo_length_max_half : ligation_site
                    + detect_oligo_length_max_half
                    + 1
                ]
        # 2. If the max length constraint is greater than the length of the oligo
        else:
            if ligation_site == (len(sequence) - ligation_site):
                detect_oligo = sequence
            elif ligation_site > (len(sequence) - ligation_site):
                start_pos = len(sequence) - 2 * min(ligation_site, len(sequence) - ligation_site) - 1
                detect_oligo = sequence[start_pos:]
            else:
                end_pos = 2 * min(ligation_site, len(sequence) - ligation_site) + 1
                detect_oligo = sequence[:end_pos]

            # 2.1 if the length of the oligo is even (only when the ligation site is exactly
            #     in the middle of an even length oligo) -> return only an even length oligo
            if (len(detect_oligo) % 2) == 0:
                detect_oligo_even = detect_oligo
            # 2.2 if the length of the oligo is odd
            else:
                # 2.2.1 if the ligation site is closer to the left -> return two different oligos: even, long right
                if (len(detect_oligo) - ligation_site) > ligation_site:
                    detect_oligo_even = detect_oligo[:-1]
                    detect_oligo_long_right = detect_oligo
                # 2.2.2 if the ligation site is closter to the right -> return two different oligos: even, long left
                else:
                    detect_oligo_even = detect_oligo[1:]
                    detect_oligo_long_left = detect_oligo

        for oligo in (detect_oligo_even, detect_oligo_long_left, detect_oligo_long_right):
            if oligo and oligo.count("T") >= min_thymines:
                return detect_oligo_even, detect_oligo_long_left, detect_oligo_long_right

        return None, None, None

    def calculate_detect_oligo(
        self,
        oligo_database: OligoDatabase,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        min_thymines: int,
        region_ids: Union[str, List[str]] = None,
    ):
        """Calculate and assign detection oligos for each oligo in the database.

        This function iterates through the oligos in the database and calculates potential detection oligo sequences
        based on specified length and thymine content constraints. The results are added to each oligo's attributes.

        :param oligo_database: Database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :param detect_oligo_length_min: The minimum length of the detection oligo.
        :type detect_oligo_length_min: int
        :param detect_oligo_length_max: The maximum length of the detection oligo.
        :type detect_oligo_length_max: int
        :param min_thymines: The minimum number of thymine bases required in the detection oligo.
        :type min_thymines: int
        :param region_ids: The region IDs for which to calculate the detection oligo.
        :type region_ids: Union[str, List[str]]
        :return: The updated oligo database with detection oligo sequences added to each oligo's attributes.
        :rtype: OligoDatabase
        """
        if region_ids is None:
            region_ids = oligo_database.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        for region_id in region_ids:
            database_region = oligo_database.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                if ("ligation_site" in oligo_attributes) and (oligo_attributes["ligation_site"]):
                    (
                        detect_oligo_even,
                        detect_oligo_long_left,
                        detect_oligo_long_right,
                    ) = self._calc_detect_oligo(
                        sequence=oligo_attributes["oligo"],
                        ligation_site=oligo_attributes["ligation_site"],
                        detect_oligo_length_min=detect_oligo_length_min,
                        detect_oligo_length_max=detect_oligo_length_max,
                        min_thymines=min_thymines,
                    )
                else:
                    detect_oligo_even = detect_oligo_long_left = detect_oligo_long_right = None
                oligo_attributes["detect_oligo_even"] = detect_oligo_even
                oligo_attributes["detect_oligo_long_left"] = detect_oligo_long_left
                oligo_attributes["detect_oligo_long_right"] = detect_oligo_long_right

        return oligo_database
