############################################
# imports
############################################

from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_calculator import BaseProperty
from oligo_designer_toolsuite.utils import cast_to_int, cast_to_string

from ._property_functions import (
    calc_detect_oligo,
    calc_dg_secondary_structure,
    calc_gc_content,
    calc_length_complement,
    calc_length_selfcomplement,
    calc_oligo_length,
    calc_padlock_arms,
    calc_seedregion,
    calc_split_sequence,
    calc_tm_nn,
    calculate_reverse_complement_sequence,
    calculate_seedregion_site,
    calculate_shortened_sequence,
)

############################################
# Sequence Property Classes
############################################


class LengthProperty(BaseProperty):
    """
    A property class for calculating the length of oligonucleotide sequences.
    """

    def __init__(self) -> None:
        """Constructor for the LengthProperty class."""
        super().__init__()

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the length of the oligonucleotide sequence.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated length property.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        length: int | None = None

        if sequence:
            length = calc_oligo_length(sequence=sequence)

        return {f"length_{sequence_type}": length}


class GCContentProperty(BaseProperty):
    """
    A property class for calculating the GC content of oligonucleotide sequences.
    """

    def __init__(self) -> None:
        """Constructor for the GCContentProperty class."""
        super().__init__()

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the GC content of the oligonucleotide sequence.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated GC content property.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        GC_content: float | None = None

        if sequence:
            GC_content = calc_gc_content(sequence=sequence)

        return {f"GC_content_{sequence_type}": GC_content}


class TmNNProperty(BaseProperty):
    """
    A property class for calculating the melting temperature (Tm) of oligonucleotide sequences using nearest-neighbor thermodynamics.

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
    """

    def __init__(
        self,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict | None = None,
        Tm_chem_correction_parameters: dict | None = None,
    ) -> None:
        """Constructor for the TmNNProperty class."""
        super().__init__()
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the melting temperature (Tm) of the oligonucleotide sequence.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated Tm property.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        TmNN: float | None = None

        if sequence:
            TmNN = calc_tm_nn(
                sequence=sequence,
                Tm_parameters=self.Tm_parameters,
                Tm_salt_correction_parameters=self.Tm_salt_correction_parameters,
                Tm_chem_correction_parameters=self.Tm_chem_correction_parameters,
            )

        return {f"TmNN_{sequence_type}": TmNN}


class DGSecondaryStructureProperty(BaseProperty):
    """
    A property class for calculating the Gibbs free energy (ΔG) of secondary structure formation for oligonucleotide sequences.

    :param T: The temperature at which the ΔG calculation should be performed.
    :type T: float
    """

    def __init__(self, T: float) -> None:
        """Constructor for the DGSecondaryStructureProperty class."""
        super().__init__()
        self.T = T

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the Gibbs free energy (ΔG) of secondary structure formation.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated ΔG property.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        DG_secondary_structure: float | None = None

        if sequence:
            DG_secondary_structure = calc_dg_secondary_structure(sequence=sequence, T=self.T)

        return {f"DG_secondary_structure_{sequence_type}": DG_secondary_structure}


class LengthSelfComplementProperty(BaseProperty):
    """
    A property class for calculating the length of self-complementary regions in oligonucleotide sequences.
    """

    def __init__(self) -> None:
        """Constructor for the LengthSelfComplementProperty class."""
        super().__init__()

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the length of the self-complementary region.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated self-complement length property.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        len_overlap: int | None = None

        if sequence:
            len_overlap = calc_length_selfcomplement(sequence=sequence)

        return {f"length_selfcomplement_{sequence_type}": len_overlap}


class LengthComplementProperty(BaseProperty):
    """
    A property class for calculating the length of complementary overlap between a sequence and a comparison sequence.

    :param comparison_sequence: The sequence to compare against for complementary overlap.
    :type comparison_sequence: str
    """

    def __init__(self, comparison_sequence: str) -> None:
        """Constructor for the LengthComplementProperty class."""
        super().__init__()
        self.comparison_sequence = comparison_sequence

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the length of complementary overlap with the comparison sequence.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated complement length property.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        len_overlap: int | None = None

        if sequence:
            len_overlap = calc_length_complement(sequence1=sequence, sequence2=self.comparison_sequence)

        return {f"length_complement_{sequence_type}_{self.comparison_sequence}": len_overlap}


class ShortenedSequenceProperty(BaseProperty):
    """
    A property class for calculating shortened versions of oligonucleotide sequences.

    :param sequence_length: The desired length for the shortened sequence.
    :type sequence_length: int
    :param reverse: If True, the shortened sequence is taken from the end of the sequence, otherwise from the beginning.
    :type reverse: bool
    """

    def __init__(self, sequence_length: int, reverse: bool = False) -> None:
        """Constructor for the ShortenedSequenceProperty class."""
        super().__init__()
        self.sequence_length = sequence_length
        self.reverse = reverse

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the shortened sequence.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated shortened sequence property.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        sequence_short: str | None = None

        if sequence:
            sequence_short = calculate_shortened_sequence(
                sequence=sequence, sequence_length=self.sequence_length, reverse=self.reverse
            )

        new_sequence_type = f"{sequence_type}_short"
        return {new_sequence_type: sequence_short}


class ReverseComplementSequenceProperty(BaseProperty):
    """
    A property class for calculating the reverse complement of oligonucleotide sequences.

    :param sequence_type_reverse_complement: The property name for storing the reverse complement sequence.
    :type sequence_type_reverse_complement: str
    """

    def __init__(self, sequence_type_reverse_complement: str) -> None:
        """Constructor for the ReverseComplementSequenceProperty class."""
        super().__init__()
        self.sequence_type_reverse_complement = sequence_type_reverse_complement

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the reverse complement sequence.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated reverse complement sequence property.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )
        sequence_rc: str | None = None

        if sequence:
            sequence_rc = calculate_reverse_complement_sequence(sequence=sequence)

        return {self.sequence_type_reverse_complement: sequence_rc}


class SplitSequenceProperty(BaseProperty):
    """
    A property class for splitting oligonucleotide sequences into sub-sequences.

    :param split_start_end: A list of tuples indicating (start, end) indices for sequence splitting.
    :type split_start_end: list[tuple]
    :param split_names: A list of names to assign to each of the split sequence segments.
    :type split_names: list[str]
    """

    def __init__(
        self,
        split_start_end: list[tuple],
        split_names: list[str],
    ) -> None:
        """Constructor for the SplitSequenceProperty class."""
        super().__init__()
        if len(split_start_end) != len(split_names):
            raise ConfigurationError(
                f"Mismatch between split sequence names and positions: {len(split_names)} names given for {len(split_start_end)} split sequences. "
                f"Must provide one name for each split sequence."
            )
        self.split_start_end = split_start_end
        self.split_names = split_names

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate split sequences from the main sequence.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated split sequence properties.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )
        properties: dict[str, str | None] = dict.fromkeys(self.split_names, None)

        if sequence:
            split_sequences = calc_split_sequence(sequence=sequence, split_start_end=self.split_start_end)
            properties = dict(zip(self.split_names, split_sequences))

        return properties


class SeedregionProperty(BaseProperty):
    """
    A property class for calculating seed regions of oligonucleotide sequences.

    :param start: The start position of the seed region. Can be an integer (exact position) or a float (fraction of the sequence length).
    :type start: int | float
    :param end: The end position of the seed region. Can be an integer (exact position) or a float (fraction of the sequence length).
    :type end: int | float
    """

    def __init__(self, start: int | float, end: int | float) -> None:
        """Constructor for the SeedregionProperty class."""
        super().__init__()
        self.start = start
        self.end = end

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the seed region positions.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated seedregion start and end positions.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        seedregion_start: int | None = None
        seedregion_end: int | None = None

        if sequence:
            seedregion_start, seedregion_end = calc_seedregion(
                sequence=sequence, start=self.start, end=self.end
            )

        return {
            "seedregion_start": seedregion_start,
            "seedregion_end": seedregion_end,
        }


class SeedregionSiteProperty(BaseProperty):
    """
    A property class for calculating seed regions around seed region sites.

    :param seedregion_size: The size of the seed region to be calculated around the seed region site.
    :type seedregion_size: int
    :param seedregion_site_name: The property name of the seed region site stored in the OligoDatabase.
    :type seedregion_site_name: str
    """

    def __init__(
        self,
        seedregion_size: int,
        seedregion_site_name: str,
    ) -> None:
        """Constructor for the SeedregionSiteProperty class."""
        super().__init__()
        self.seedregion_size = seedregion_size
        self.seedregion_site_name = seedregion_site_name

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the seed region around the seed region site.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated seedregion start and end positions.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )
        seedregion_site = cast_to_int(
            oligo_database.get_oligo_property_value(
                property=self.seedregion_site_name, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        seedregion_start: int | None = None
        seedregion_end: int | None = None

        if sequence and seedregion_site:
            seedregion_start, seedregion_end = calculate_seedregion_site(
                sequence=sequence, seedregion_site=seedregion_site, seedregion_size=self.seedregion_size
            )

        return {
            "seedregion_start": seedregion_start,
            "seedregion_end": seedregion_end,
        }


class PadlockArmsProperty(BaseProperty):
    """
    A property class for calculating padlock probe arms and ligation sites.

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
    """

    def __init__(
        self,
        arm_length_min: int,
        arm_Tm_dif_max: float,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters: dict,
        Tm_salt_correction_parameters: dict | None = None,
        Tm_chem_correction_parameters: dict | None = None,
    ) -> None:
        """Constructor for the PadlockArmsProperty class."""
        super().__init__()
        self.arm_length_min = arm_length_min
        self.arm_Tm_dif_max = arm_Tm_dif_max
        self.arm_Tm_min = arm_Tm_min
        self.arm_Tm_max = arm_Tm_max
        self.Tm_parameters = Tm_parameters
        self.Tm_salt_correction_parameters = Tm_salt_correction_parameters
        self.Tm_chem_correction_parameters = Tm_chem_correction_parameters

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the padlock probe arms and ligation site.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated arm temperatures and ligation site.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )
        arm1_Tm: float | None = None
        arm2_Tm: float | None = None
        ligation_site: int | None = None

        if sequence:
            arm1_Tm, arm2_Tm, ligation_site = calc_padlock_arms(
                sequence=sequence,
                arm_length_min=self.arm_length_min,
                arm_Tm_dif_max=self.arm_Tm_dif_max,
                arm_Tm_min=self.arm_Tm_min,
                arm_Tm_max=self.arm_Tm_max,
                Tm_parameters=self.Tm_parameters,
                Tm_salt_correction_parameters=self.Tm_salt_correction_parameters,
                Tm_chem_correction_parameters=self.Tm_chem_correction_parameters,
            )
        return {
            "arm1_Tm": arm1_Tm,
            "arm2_Tm": arm2_Tm,
            "ligation_site": ligation_site,
        }


class DetectOligoProperty(BaseProperty):
    """
    A property class for calculating detection oligos around ligation sites.

    :param detect_oligo_length_min: The minimum allowable length for the detection oligo.
    :type detect_oligo_length_min: int
    :param detect_oligo_length_max: The maximum allowable length for the detection oligo.
    :type detect_oligo_length_max: int
    :param min_thymines: The minimum number of thymine bases required in the detection oligo.
    :type min_thymines: int
    :param ligation_site_name: The property name of the ligation site stored in the OligoDatabase.
    :type ligation_site_name: str
    """

    def __init__(
        self,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        min_thymines: int,
        ligation_site_name: str = "ligation_site",
    ) -> None:
        """Constructor for the DetectOligoProperty class."""
        super().__init__()
        self.detect_oligo_length_min = detect_oligo_length_min
        self.detect_oligo_length_max = detect_oligo_length_max
        self.min_thymines = min_thymines
        self.ligation_site_name = ligation_site_name

    def apply(self, oligo_database: OligoDatabase, region_id: str, oligo_id: str, sequence_type: str) -> dict:
        """
        Calculate the detection oligo sequences around the ligation site.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property is calculated.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A dictionary containing the calculated detection oligo sequences.
        :rtype: dict
        """
        sequence = cast_to_string(
            oligo_database.get_oligo_property_value(
                property=sequence_type, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )
        ligation_site = cast_to_int(
            oligo_database.get_oligo_property_value(
                property=self.ligation_site_name, region_id=region_id, oligo_id=oligo_id, flatten=True
            )
        )

        detect_oligo_even: str | None = None
        detect_oligo_long_left: str | None = None
        detect_oligo_long_right: str | None = None

        if sequence and ligation_site:
            (
                detect_oligo_even,
                detect_oligo_long_left,
                detect_oligo_long_right,
            ) = calc_detect_oligo(
                sequence=sequence,
                ligation_site=ligation_site,
                detect_oligo_length_min=self.detect_oligo_length_min,
                detect_oligo_length_max=self.detect_oligo_length_max,
                min_thymines=self.min_thymines,
            )

        return {
            "detect_oligo_even": detect_oligo_even,
            "detect_oligo_long_left": detect_oligo_long_left,
            "detect_oligo_long_right": detect_oligo_long_right,
        }
