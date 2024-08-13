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
from oligo_designer_toolsuite.utils import check_if_list, flatten_attribute_list

############################################
# Attrubite Calculation Class
############################################


class OligoAttributes:

    def __init__(self) -> None:
        """Constructor for the OligoAttributes class."""

    @staticmethod
    def _calc_oligo_length(sequence: str) -> int:
        return len(sequence)

    def calculate_oligo_length(
        self, oligo_database: OligoDatabase, region_ids: Union[str, List[str]] = None
    ) -> OligoDatabase:

        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                # oligo and target have always same length
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute="oligo", region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                length = self._calc_oligo_length(sequence=sequence)
                new_oligo_attribute[oligo_id] = {
                    "length": length,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_num_targeted_transcripts(transcript_id: list) -> int:
        # make sure that transcript id one level list
        return len(set(check_if_list(flatten_attribute_list(transcript_id))))

    def calculate_num_targeted_transcripts(
        self, oligo_database: OligoDatabase, region_ids: Union[str, List[str]] = None
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                transcript_id = oligo_database.get_oligo_attribute_value(
                    attribute="transcript_id", region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                if transcript_id:
                    num_targeted_transcripts = self._calc_num_targeted_transcripts(
                        transcript_id=transcript_id
                    )
                else:
                    num_targeted_transcripts = None

                new_oligo_attribute[oligo_id] = {
                    "num_targeted_transcripts": num_targeted_transcripts,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_isoform_consensus(transcript_id: list, number_total_transcripts: list) -> float:

        # number transcripts is the number of transcripts of a genomic region
        # hence, all values have to be the same for each transcript coming from the same oligo
        # since only oligos from the same genomic region are merged into one entry
        number_total_transcripts = int(check_if_list(number_total_transcripts)[0])
        num_targeted_transcripts = len(set(check_if_list(transcript_id)))
        isoform_consensus = round(num_targeted_transcripts / number_total_transcripts * 100, 2)

        return isoform_consensus

    def calculate_isoform_consensus(
        self, oligo_database: OligoDatabase, region_ids: Union[str, List[str]] = None
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                number_total_transcripts = oligo_database.get_oligo_attribute_value(
                    attribute="number_total_transcripts", region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                transcript_id = oligo_database.get_oligo_attribute_value(
                    attribute="transcript_id", region_id=region_id, oligo_id=oligo_id, flatten=True
                )

                if transcript_id and number_total_transcripts:
                    isoform_consensus = self._calc_isoform_consensus(
                        transcript_id=transcript_id, number_total_transcripts=number_total_transcripts
                    )
                else:
                    isoform_consensus = None

                new_oligo_attribute[oligo_id] = {
                    "isoform_consensus": isoform_consensus,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_seedregion(sequence: str, start: Union[int, float], end: Union[int, float]) -> tuple[int, int]:

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
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:

        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                seedregion_start, seedregion_end = self._calc_seedregion(
                    sequence=sequence, start=start, end=end
                )
                new_oligo_attribute[oligo_id] = {
                    "seedregion_start": seedregion_start,
                    "seedregion_end": seedregion_end,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_seedregion_ligationsite(
        sequence: str, ligation_site: int, seedregion_size: int
    ) -> tuple[int, int]:
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
        self,
        oligo_database: OligoDatabase,
        seedregion_size: int,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:

        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                ligation_site = oligo_database.get_oligo_attribute_value(
                    attribute="ligation_site", region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                if ligation_site:
                    # oligo and target have always same length
                    seedregion_start, seedregion_end = self._calc_seedregion_ligationsite(
                        sequence=sequence, ligation_site=ligation_site, seedregion_size=seedregion_size
                    )
                else:
                    seedregion_start = seedregion_end = None

                new_oligo_attribute[oligo_id] = {
                    "seedregion_start": seedregion_start,
                    "seedregion_end": seedregion_end,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_GC_content(sequence: str) -> float:
        return round(gc_fraction(sequence) * 100, 2)

    def calculate_GC_content(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                GC_content = self._calc_GC_content(sequence=sequence)
                new_oligo_attribute[oligo_id] = {
                    "GC_content": GC_content,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_TmNN(
        sequence: str,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
    ) -> float:
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
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict = None,
        Tm_chem_correction_parameters: dict = None,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                TmNN = self._calc_TmNN(
                    sequence=sequence,
                    Tm_parameters=Tm_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                )
                new_oligo_attribute[oligo_id] = {
                    "TmNN": TmNN,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_length_complement(sequence1: str, sequence2: str) -> int:
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
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                # we want to check if the reverse of our sequence is complementary to itself, e.g.
                # 5' - TAA CAA TAT ATA TTG TTA - 3' and it's reverse
                # 3' - ATT GTT ATA TAT AAC AAT - 5' are complementary to each other
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                sequence_rev = sequence[::-1]
                len_overlap = self._calc_length_complement(sequence1=sequence, sequence2=sequence_rev)
                new_oligo_attribute[oligo_id] = {
                    "length_selfcomplement": len_overlap,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    def calculate_length_complement(
        self,
        oligo_database: OligoDatabase,
        comparison_sequence: str,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                len_overlap = self._calc_length_complement(sequence1=sequence, sequence2=comparison_sequence)
                new_oligo_attribute[oligo_id] = {
                    f"length_complement_{comparison_sequence}": len_overlap,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_DG_secondary_structure(sequence: str, T: float) -> float:
        return dg(sequence, temp=T)

    def calculate_DG_secondary_structure(
        self,
        oligo_database: OligoDatabase,
        T: float,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                DG_secondary_structure = self._calc_DG_secondary_structure(sequence=sequence, T=T)
                new_oligo_attribute[oligo_id] = {
                    "DG_secondary_structure": DG_secondary_structure,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

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
    ) -> tuple[float, float, int]:
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
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                arm1_Tm, arm2_Tm, ligation_site = self._calc_padlock_arms(
                    sequence=sequence,
                    arm_length_min=arm_length_min,
                    arm_Tm_dif_max=arm_Tm_dif_max,
                    arm_Tm_min=arm_Tm_min,
                    arm_Tm_max=arm_Tm_max,
                    Tm_parameters=Tm_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                )
                new_oligo_attribute[oligo_id] = {
                    "arm1_Tm": arm1_Tm,
                    "arm2_Tm": arm2_Tm,
                    "ligation_site": ligation_site,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    def _calc_detect_oligo(
        self,
        sequence: str,
        ligation_site: int,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        min_thymines: int,
    ) -> tuple[str, str, str]:
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
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                ligation_site = oligo_database.get_oligo_attribute_value(
                    attribute="ligation_site", region_id=region_id, oligo_id=oligo_id, flatten=True
                )

                if ligation_site:
                    (
                        detect_oligo_even,
                        detect_oligo_long_left,
                        detect_oligo_long_right,
                    ) = self._calc_detect_oligo(
                        sequence=sequence,
                        ligation_site=ligation_site,
                        detect_oligo_length_min=detect_oligo_length_min,
                        detect_oligo_length_max=detect_oligo_length_max,
                        min_thymines=min_thymines,
                    )
                else:
                    detect_oligo_even = detect_oligo_long_left = detect_oligo_long_right = None

                new_oligo_attribute[oligo_id] = {
                    "detect_oligo_even": detect_oligo_even,
                    "detect_oligo_long_left": detect_oligo_long_left,
                    "detect_oligo_long_right": detect_oligo_long_right,
                }
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database
