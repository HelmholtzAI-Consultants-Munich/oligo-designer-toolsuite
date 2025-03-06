############################################
# imports
############################################

from typing import List, Tuple, Union

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp, gc_fraction
from seqfold import dg

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.utils import check_if_list, flatten_attribute_list

############################################
# Attrubite Calculation Class
############################################


class OligoAttributes:
    """
    The OligoAttributes class provides a comprehensive set of methods to calculate and update various attributes of oligonucleotides within a given OligoDatabase.

    This class includes functionalities for determining oligo length, GC content, melting temperature (Tm), secondary structure stability (ΔG), isoform consensus, and more.
    These calculations are essential for the design and evaluation of oligos forcertain experimental specifications.
    """

    def __init__(self) -> None:
        """Constructor for the OligoAttributes class."""

    @staticmethod
    def _calc_oligo_length(sequence: str) -> int:
        """Calculate the length of an oligonucleotide sequence.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :return: The length of the sequence.
        :rtype: int
        """
        return len(sequence)

    def calculate_oligo_length(
        self, oligo_database: OligoDatabase, region_ids: Union[str, List[str]] = None
    ) -> OligoDatabase:
        """Calculate and update the length of oligonucleotides for the specified regions of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """

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
    def _calculate_shortened_sequence(sequence: str, sequence_length: int, reverse: bool) -> int:
        """Calculate the shortened sequence of an oligonucleotide sequence.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :param sequence_length: The desired length for the shortened sequence.
        :type sequence_length: int
        :param reverse: If True, the shortened sequence is taken from the end of the sequence, otherwise from the beginning.
        :type reverse: bool
        :return: The shortened sequence.
        :rtype: str
        """
        return sequence[:sequence_length] if not reverse else sequence[-sequence_length:]

    def calculate_shortened_sequence(
        self,
        oligo_database: OligoDatabase,
        sequence_length: int,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
        reverse: bool = False,
    ) -> OligoDatabase:
        """Calculate and update the shortened oligonucleotide sequences for the specified regions of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_length: The desired length for the shortened sequence.
        :type sequence_length: int
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :param reverse: If True, the shortened sequence is taken from the end of the sequence, otherwise from the beginning, defaults to False.
        :type reverse: bool, optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )

                sequence_short = self._calculate_shortened_sequence(
                    sequence=sequence, sequence_length=sequence_length, reverse=reverse
                )
                new_oligo_attribute[oligo_id] = {f"{sequence_type}_short": sequence_short}
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calculate_reverse_complement_sequence(sequence: str) -> int:
        """Calculate the reverse complemented sequence of an oligonucleotide sequence.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :return: The reverse complemented sequence.
        :rtype: str
        """
        return str(Seq(sequence).reverse_complement())

    def calculate_reverse_complement_sequence(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        sequence_type_reverse_complement: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        """Calculate and update the reverse complement of oligonucleotide sequences for the specified regions of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for attribute calculation.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param sequence_type_reverse_complement: The type of sequence of the reverse complement.
        :type sequence_type_reverse_complement: int
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()
        new_oligo_attribute = {}

        for region_id in region_ids:
            for oligo_id in oligo_database.database[region_id].keys():
                sequence = oligo_database.get_oligo_attribute_value(
                    attribute=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
                )

                sequence_rc = self._calculate_reverse_complement_sequence(sequence=sequence)
                new_oligo_attribute[oligo_id] = {sequence_type_reverse_complement: sequence_rc}
        oligo_database.update_oligo_attributes(new_oligo_attribute)

        return oligo_database

    @staticmethod
    def _calc_num_targeted_transcripts(transcript_id: list) -> int:
        """Calculate the number of unique transcripts targeted by an oligonucleotide.

        :param transcript_id: List of transcript IDs associated with the oligonucleotide.
        :type transcript_id: list
        :return: Number of unique targeted transcripts.
        :rtype: int
        """
        # make sure that transcript id one level list
        return len(set(check_if_list(flatten_attribute_list(transcript_id))))

    def calculate_num_targeted_transcripts(
        self, oligo_database: OligoDatabase, region_ids: Union[str, List[str]] = None
    ) -> OligoDatabase:
        """Calculate and update the number of targeted transcripts for each oligonucleotide for the specified regions of the OligoDatabase.

        If the required information for `_calc_num_targeted_transcripts` is not available,
        the 'num_targeted_transcripts' attribute is set to None.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
        """Calculate the isoform consensus for an oligonucleotide, representing the percentage of transcripts
        targeted by the oligo out of the total number of transcripts in a region. The maximum value for the
        isoform consensus is 100%, which means that the oligo targets all isoforms (transcripts) of the region.

        :param transcript_id: List of transcript IDs associated with the oligonucleotide.
        :type transcript_id: list
        :param number_total_transcripts: Total number of transcripts in the genomic region.
        :type number_total_transcripts: list
        :return: Isoform consensus as a percentage.
        :rtype: float
        """
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
        """Calculate and update the isoform consensus for each oligonucleotide for the specified regions of the OligoDatabase.

        If the required information for `_calc_isoform_consensus` is not available,
        the 'isoform_consensus' attribute is set to None.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
    def _calc_seedregion(sequence: str, start: Union[int, float], end: Union[int, float]) -> Tuple[int, int]:
        """Calculate the seed region of a nucleotide sequence based on the provided start and end positions.

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

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :param start: The start position of the seed region. Can be an integer (exact position) or a float (fraction of the sequence length).
        :type start: Union[int, float]
        :param end: The end position of the seed region. Can be an integer (exact position) or a float (fraction of the sequence length).
        :type end: Union[int, float]
        :return: A tuple containing the start and end positions of the seed region.
        :rtype: Tuple[int, int]
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
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        """Calculate and update the seed region for each oligonucleotide for the specified regions of the OligoDatabase based on the given start and end positions.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param start: The start position of the seed region. Can be an integer (exact position) or a float (fraction of the sequence length).
        :type start: Union[int, float]
        :param end: The end position of the seed region. Can be an integer (exact position) or a float (fraction of the sequence length).
        :type end: Union[int, float]
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
    ) -> Tuple[int, int]:
        """Calculate the start and end positions of the seed region around a ligation site for a nucleotide sequence.
        The seed region is defined symmetrically around the ligation site, considering the provided `seedregion_size`.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :param ligation_site: The position of the ligation site within the sequence.
        :type ligation_site: int
        :param seedregion_size: The size of the seed region to be calculated.
        :type seedregion_size: int
        :return: A tuple containing the start and end positions of the seed region.
        :rtype: Tuple[int, int]
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
        self,
        oligo_database: OligoDatabase,
        seedregion_size: int,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        """Calculate and update the seed region around the ligation site for each oligonucleotide for the
        specified regions of the OligoDatabase. If the required information for `_calc_seedregion_ligationsite`
        is not available, the 'seedregion_start' and 'seedregion_end' attributes are set to None.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param seedregion_size: The size of the seed region to be calculated around the ligation site.
        :type seedregion_size: int
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
        """Calculate the GC content of a given nucleotide sequence.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :return: The GC content as a percentage, rounded to two decimal places.
        :rtype: float
        """
        return round(gc_fraction(sequence) * 100, 2)

    def calculate_GC_content(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        """Calculate and update the GC content for each oligonucleotide for the specified regions of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
        """Calculate the melting temperature (Tm) of a nucleotide sequence using nearest-neighbor thermodynamics.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :param Tm_parameters: Parameters for the nearest-neighbor Tm calculation.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type Tm_parameters: dict
        :param Tm_salt_correction_parameters: Optional parameters for salt correction.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type Tm_salt_correction_parameters: dict, optional
        :param Tm_chem_correction_parameters: Optional parameters for chemical correction.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type Tm_chem_correction_parameters: dict, optional
        :return: The calculated melting temperature (Tm) in degrees Celsius, rounded to two decimal places.
        :rtype: float
        """
        TmNN = MeltingTemp.Tm_NN(sequence, **Tm_parameters)
        if Tm_salt_correction_parameters is not None:
            TmNN += MeltingTemp.salt_correction(**Tm_salt_correction_parameters, seq=sequence)
        if Tm_chem_correction_parameters is not None:
            TmNN = MeltingTemp.chem_correction(TmNN, **Tm_chem_correction_parameters)
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
        """Calculate and update the melting temperature (Tm) for each oligonucleotide for the specified regions of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param Tm_parameters: Parameters for the nearest-neighbor Tm calculation.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type Tm_parameters: dict
        :param Tm_salt_correction_parameters: Optional parameters for salt correction.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type Tm_salt_correction_parameters: dict, optional
        :param Tm_chem_correction_parameters: Optional parameters for chemical correction.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type Tm_chem_correction_parameters: dict, optional
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
        """Calculate the maximum length of complementary overlap between two sequences.

        :param sequence1: The first nucleotide sequence.
        :type sequence1: str
        :param sequence2: The second nucleotide sequence.
        :type sequence2: str
        :return: The maximum length of the complementary overlap between the sequences.
        :rtype: int
        """

        def _calculate_max_overlap(seq1: str, seq2: str) -> int:
            """
            Calculate the maximum overlap between two sequences.

            This function compares two sequences and determines the maximum length of consecutive matching characters between them.
            It iterates over the characters of both sequences and tracks the length of overlapping substrings. The longest such overlap is returned.

            :param seq1: The first sequence to compare.
            :type seq1: str
            :param seq2: The second sequence to compare.
            :type seq2: str
            :return: The length of the maximum overlap between the two sequences.
            :rtype: int
            """
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
        """Calculate and update the length of the maximum self-complementary region for each oligonucleotide for the specified regions of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
        """Calculate and update the length of the complementary overlap between a specified comparison sequence and each oligonucleotide for the specified regions of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param comparison_sequence: The sequence to compare against for complementary overlap.
        :type comparison_sequence: str
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
        """Calculate the Gibbs free energy (ΔG) of the secondary structure of a nucleotide sequence at a given temperature.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :param T: The temperature at which the ΔG calculation should be performed.
        :type T: float
        :return: The calculated ΔG value.
        :rtype: float
        """
        return dg(sequence, temp=T)

    def calculate_DG_secondary_structure(
        self,
        oligo_database: OligoDatabase,
        T: float,
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        """Calculate and update the Gibbs free energy (ΔG) for secondary structure formation at a specified temperature for each oligonucleotide for the specified regions of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param T: The temperature in degrees Celsius at which to calculate the ΔG.
        :type T: float
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
    ) -> Tuple[float, float, int]:
        """Calculate the melting temperatures (Tm) of padlock probe arms and determine the ligation site.

        This function evaluates potential padlock probe arms in a given sequence by calculating their melting temperatures (Tm)
        and finding the optimal ligation site that satisfies the specified Tm conditions. It iteratively adjusts the ligation site
        to find arm lengths and Tm values that meet the criteria of minimum arm length and Tm differences within the specified maximum.
        The process stops once suitable arms are found or if no configuration meets the criteria.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :param arm_length_min: The minimum length required for each arm of the padlock probe.
        :type arm_length_min: int
        :param arm_Tm_dif_max: The maximum allowable difference between the Tm of the two arms.
        :type arm_Tm_dif_max: float
        :param arm_Tm_min: The minimum allowable Tm for each arm.
        :type arm_Tm_min: float
        :param arm_Tm_max: The maximum allowable Tm for each arm.
        :type arm_Tm_max: float
        :param Tm_parameters: Parameters for the nearest-neighbor Tm calculation.
        :type Tm_parameters: dict
        :param Tm_salt_correction_parameters: Optional parameters for salt correction.
        :type Tm_salt_correction_parameters: dict, optional
        :param Tm_chem_correction_parameters: Optional parameters for chemical correction.
        :type Tm_chem_correction_parameters: dict, optional
        :return: A tuple containing the Tm of the first arm, the Tm of the second arm, and the ligation site.
                Returns (None, None, None) if no valid ligation site is found.
        :rtype: Tuple[float, float, int]
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
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        """Calculate and update the melting temperatures (Tm) of the padlock probe arms for each oligonucleotide for the specified regions of the OligoDatabase,
        and identify the optimal ligation site based on specified Tm constraints.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param arm_length_min: The minimum length of each arm of the padlock probe.
        :type arm_length_min: int
        :param arm_Tm_dif_max: The maximum allowable difference in Tm between the two arms.
        :type arm_Tm_dif_max: float
        :param arm_Tm_min: The minimum allowable Tm for each arm.
        :type arm_Tm_min: float
        :param arm_Tm_max: The maximum allowable Tm for each arm.
        :type arm_Tm_max: float
        :param Tm_parameters: Parameters for the nearest-neighbor Tm calculation.
        :type Tm_parameters: dict
        :param Tm_salt_correction_parameters: Optional parameters for salt correction.
        :type Tm_salt_correction_parameters: dict, optional
        :param Tm_chem_correction_parameters: Optional parameters for chemical correction.
        :type Tm_chem_correction_parameters: dict, optional
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
    ) -> Tuple[str, str, str]:
        """Calculate potential detection oligos around a ligation site, ensuring they meet specified length and thymine content criteria.

        :param sequence: The nucleotide sequence.
        :type sequence: str
        :param ligation_site: The position of the ligation site within the sequence.
        :type ligation_site: int
        :param detect_oligo_length_min: The minimum allowable length for the detection oligo.
        :type detect_oligo_length_min: int
        :param detect_oligo_length_max: The maximum allowable length for the detection oligo.
        :type detect_oligo_length_max: int
        :param min_thymines: The minimum number of thymine bases required in the detection oligo.
        :type min_thymines: int
        :return: A tuple containing the even-length detection oligo, and possibly longer left and right versions, or None if conditions aren't met.
        :rtype: Tuple[str, str, str]
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
        sequence_type: _TYPES_SEQ = "oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> OligoDatabase:
        """Calculate and update the detection oligo sequences for each oligonucleotide for the specified regions of the OligoDatabase, based on ligation site position and length constraints.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param detect_oligo_length_min: The minimum allowable length for the detection oligo.
        :type detect_oligo_length_min: int
        :param detect_oligo_length_max: The maximum allowable length for the detection oligo.
        :type detect_oligo_length_max: int
        :param min_thymines: The minimum number of thymine bases required in the detection oligo.
        :type min_thymines: int
        :param sequence_type: The type of sequence to be used for attribute calculation, defaults to "oligo".
        :type sequence_type: _TYPES_SEQ["oligo", "target"], optional
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The updated OligoDatabase with the calculated attribute.
        :rtype: OligoDatabase
        """
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
