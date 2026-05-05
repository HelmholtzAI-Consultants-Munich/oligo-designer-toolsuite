############################################
# imports
############################################

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp, gc_fraction
from seqfold import dg

from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.utils import cast_to_list, flatten_property_list

############################################
# Property Calculation Functions
############################################


def calc_oligo_length(sequence: str) -> int:
    """Calculate the length of an oligonucleotide sequence.

    :param sequence: The nucleotide sequence.
    :type sequence: str
    :return: The length of the sequence.
    :rtype: int
    """
    return len(sequence)


def calc_gc_content(sequence: str) -> float:
    """Calculate the GC content of a given nucleotide sequence.

    :param sequence: The nucleotide sequence.
    :type sequence: str
    :return: The GC content as a percentage, rounded to two decimal places.
    :rtype: float
    """
    gc_content: float = round(gc_fraction(sequence) * 100, 2)
    return gc_content


def calc_tm_nn(
    sequence: str,
    Tm_parameters: dict,
    Tm_salt_correction_parameters: dict | None = None,
    Tm_chem_correction_parameters: dict | None = None,
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
    :type Tm_salt_correction_parameters: dict | None, optional
    :param Tm_chem_correction_parameters: Optional parameters for chemical correction.
        For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
        see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
    :type Tm_chem_correction_parameters: dict | None, optional
    :return: The calculated melting temperature (Tm) in degrees Celsius, rounded to two decimal places.
    :rtype: float
    """
    TmNN: float = MeltingTemp.Tm_NN(sequence, **Tm_parameters)
    if Tm_salt_correction_parameters is not None:
        TmNN += MeltingTemp.salt_correction(**Tm_salt_correction_parameters, seq=sequence)
    if Tm_chem_correction_parameters is not None:
        TmNN = MeltingTemp.chem_correction(TmNN, **Tm_chem_correction_parameters)
    TmNN = round(TmNN, 2)
    return TmNN


def calc_dg_secondary_structure(sequence: str, T: float) -> float:
    """Calculate the Gibbs free energy (ΔG) of the secondary structure of a nucleotide sequence at a given temperature.

    :param sequence: The nucleotide sequence.
    :type sequence: str
    :param T: The temperature at which the ΔG calculation should be performed.
    :type T: float
    :return: The calculated ΔG value.
    :rtype: float
    """
    DG_secondary_structure: float = dg(sequence, temp=T)
    return DG_secondary_structure


def calc_length_complement(sequence1: str, sequence2: str) -> int:
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


def calc_length_selfcomplement(sequence: str) -> int:
    """Calculate the maximum length of self-complementary region in a sequence.

    :param sequence: The nucleotide sequence.
    :type sequence: str
    :return: The maximum length of the self-complementary region.
    :rtype: int
    """
    sequence_rev = sequence[::-1]
    return calc_length_complement(sequence1=sequence, sequence2=sequence_rev)


def calculate_shortened_sequence(sequence: str, sequence_length: int, reverse: bool) -> str:
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


def calculate_reverse_complement_sequence(sequence: str) -> str:
    """Calculate the reverse complemented sequence of an oligonucleotide sequence.

    :param sequence: The nucleotide sequence.
    :type sequence: str
    :return: The reverse complemented sequence.
    :rtype: str
    """
    return str(Seq(sequence).reverse_complement())


def calc_split_sequence(sequence: str, split_start_end: list[tuple]) -> list[str]:
    """
    Extracts sub-sequences from a given sequence using a list of (start, end) index pairs.

    :param sequence: The full sequence string to be split.
    :type sequence: str
    :param split_start_end: A list of tuples indicating start and end indices for each subsequence.
    :type split_start_end: list[tuple]
    :return: List of sub-sequences extracted from the input sequence.
    :rtype: list[str]
    """
    split_sequences = []
    for start_end in split_start_end:
        split_sequences.append(sequence[start_end[0] : start_end[1]])

    return split_sequences


def calc_seedregion(sequence: str, start: int | float, end: int | float) -> tuple[int, int]:
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
    :type start: int | float
    :param end: The end position of the seed region. Can be an integer (exact position) or a float (fraction of the sequence length).
    :type end: int | float
    :return: A tuple containing the start and end positions of the seed region.
    :rtype: tuple[int, int]
    """
    length = len(sequence)

    if isinstance(start, int) and isinstance(end, int):
        seedregion_start = max(0, start)
        seedregion_end = min(length, end)
    elif isinstance(start, float) and isinstance(end, float):
        if (not 0 <= start <= 1) or (not 0 <= end <= 1):
            raise ConfigurationError(
                f"Start and end positions must be in the interval [0,1] for float type. "
                f"Received: start={start}, end={end}."
            )
        seedregion_start = int(round(start * length))
        seedregion_end = int(round(end * length))
    else:
        raise ConfigurationError(
            f"Start and end parameters must be both integers or both floats. "
            f"Received: start={type(start).__name__}, end={type(end).__name__}."
        )

    return seedregion_start, seedregion_end


def calculate_seedregion_site(sequence: str, seedregion_site: int, seedregion_size: int) -> tuple[int, int]:
    """Calculate the start and end positions of the seed region around a seed region site for a nucleotide sequence.
    The seed region is defined symmetrically around the seed region site, considering the provided `seedregion_size`.

    :param sequence: The nucleotide sequence.
    :type sequence: str
    :param seedregion_site: The position of the seed region site within the sequence.
    :type seedregion_site: int
    :param seedregion_size: The size of the seed region to be calculated.
    :type seedregion_size: int
    :return: A tuple containing the start and end positions of the seed region.
    :rtype: tuple[int, int]
    """
    length = len(sequence)

    seedregion_start = int(max(0, seedregion_site - (seedregion_size - 1)))
    seedregion_end = int(
        min(
            length,
            seedregion_site + seedregion_size,
        )
    )
    return seedregion_start, seedregion_end


def calc_padlock_arms(
    sequence: str,
    arm_length_min: int,
    arm_Tm_dif_max: float,
    arm_Tm_min: float,
    arm_Tm_max: float,
    Tm_parameters: dict,
    Tm_salt_correction_parameters: dict | None = None,
    Tm_chem_correction_parameters: dict | None = None,
) -> tuple[float | None, float | None, int | None]:
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
    :type Tm_salt_correction_parameters: dict | None, optional
    :param Tm_chem_correction_parameters: Optional parameters for chemical correction.
    :type Tm_chem_correction_parameters: dict | None, optional
    :return: A tuple containing the Tm of the first arm, the Tm of the second arm, and the ligation site.
            Returns (None, None, None) if no valid ligation site is found.
    :rtype: tuple[float | None, float | None, int | None]
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
        Tm_arm1 = calc_tm_nn(
            sequence[:ligation_site],
            Tm_parameters,
            Tm_salt_correction_parameters,
            Tm_chem_correction_parameters,
        )
        Tm_arm2 = calc_tm_nn(
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


def calc_detect_oligo(
    sequence: str,
    ligation_site: int,
    detect_oligo_length_min: int,
    detect_oligo_length_max: int,
    min_thymines: int,
) -> tuple[str | None, str | None, str | None]:
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
    :rtype: tuple[str | None, str | None, str | None]
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


def calc_num_targeted_transcripts(transcript_id: list) -> int:
    """Calculate the number of unique transcripts targeted by an oligonucleotide.

    :param transcript_id: list of transcript IDs associated with the oligonucleotide.
    :type transcript_id: list
    :return: Number of unique targeted transcripts.
    :rtype: int
    """
    # make sure that transcript id one level list
    return len(set(flatten_property_list(transcript_id)))


def calc_isoform_consensus(transcript_id: list, number_total_transcripts: list) -> float:
    """Calculate the isoform consensus for an oligonucleotide, representing the percentage of transcripts
    targeted by the oligo out of the total number of transcripts in a region. The maximum value for the
    isoform consensus is 100%, which means that the oligo targets all isoforms (transcripts) of the region.

    :param transcript_id: list of transcript IDs associated with the oligonucleotide.
    :type transcript_id: list
    :param number_total_transcripts: Total number of transcripts in the genomic region. All values in this list should be the same since they come from the same genomic region.
    :type number_total_transcripts: list
    :return: Isoform consensus as a percentage.
    :rtype: float
    """
    # Extract the total number of transcripts (all values in the list are the same for the same region)
    total_transcripts = int(number_total_transcripts[0])

    # Count unique targeted transcripts
    unique_transcript_ids = set(cast_to_list(transcript_id))
    num_targeted = len(unique_transcript_ids)

    # Calculate percentage: (targeted / total) * 100
    if total_transcripts == 0:
        return 0.0

    isoform_consensus = round((num_targeted / total_transcripts) * 100, 2)
    return isoform_consensus
