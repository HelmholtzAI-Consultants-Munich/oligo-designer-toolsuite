############################################
# imports
############################################

import os
import warnings
from collections import defaultdict
from pathlib import Path
from typing import List, Union, get_args

import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


from .._constants import _TYPES_SEQ, SEPARATOR_OLIGO_ID
from ..utils._checkers import check_if_key_exists, check_if_list, check_tsv_format
from ..utils._database_processor import (
    check_if_region_in_database,
    collapse_info_for_duplicated_sequences,
    filter_dabase_for_region,
    merge_databases,
)
from ..utils._sequence_parser import FastaParser

############################################
# Oligo Database Class
############################################


class OligoDatabase:
    """Class for managing the oligo databases. This class provides functionality for handling oligo databases,
    allowing users to load, manipulate, and save oligo information efficiently. It includes methods for
    loading metadata, loading and managing oligo databases from various sources (e.g. fasta file or saved database),
    and performing operations such as removing regions with insufficient oligos and calculating the number of targeted
    transcripts for each oligo.

    The header of each sequence must start with '>' and contain the following information:
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
    where the region_id is compulsory and the other fileds are opional.

    Input Format (per sequence):
    >region_id::additional information::chromosome:start-end(strand)
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    Moreover, the database can be saved and loaded to/from a tsv file.

    :param min_oligos_per_region: Minimum number of oligos required per region (default is 0).
    :type min_oligos_per_region: int, optional
    :param write_regions_with_insufficient_oligos: Flag to enable writing regions with insufficient oligos to a file (default is True).
    :type write_regions_with_insufficient_oligos: bool, optional
    :param dir_output: Directory path for the output (default is "output").
    :type dir_output: str, optional
    """

    def __init__(
        self,
        min_oligos_per_region: int = 0,
        write_regions_with_insufficient_oligos: bool = True,
        dir_output: str = "output",
    ):
        """Constructor for the OligoDatabase class."""
        self.min_oligos_per_region = min_oligos_per_region
        self.write_regions_with_insufficient_oligos = write_regions_with_insufficient_oligos

        self.metadata = {}

        self.dir_output = os.path.abspath(os.path.join(dir_output, "oligo_database"))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.fasta_parser = FastaParser()

        # Initialize databse object
        self.database = {}
        self.oligosets = {}  # will be used later in the gereration of non overlpping sets

        # Initialize the file for regions with insufficient oligos
        if self.write_regions_with_insufficient_oligos:
            self.file_removed_regions = os.path.join(self.dir_output, "regions_with_insufficient_oligos.txt")
            with open(self.file_removed_regions, "a") as handle:
                handle.write(f"Region\tPipeline step\n")

    ############################################
    # Load Functions
    ############################################

    def load_metadata(self, metadata: Union[str, dict]) -> None:
        """Load metadata into the OligoDatabase object.

        If metadata already exists, a warning is issued about overwriting the existing metadata. The new metadata can be
        provided either as a path to a YAML file or directly as a dictionary.

        :param metadata: Path to a YAML file or a dictionary containing metadata information.
        :type metadata: Union[str, dict]

        :raises ValueError: If metadata has an incorrect format.
        """
        if self.metadata:
            warnings.warn("Metadata not empty! Overwriting metadata with new metadata from file!")

        if type(metadata) is str and os.path.exists(metadata):
            with open(metadata) as handle:
                self.metadata = yaml.safe_load(handle)
        elif type(metadata) is dict:
            self.metadata = metadata
        else:
            raise ValueError("Metadat has icorrect format!")

    def load_database(
        self,
        file_database: str,
        region_ids: Union[str, List[str]] = None,
        database_overwrite: bool = False,
    ) -> None:
        """Load a previously saved oligo database from a TSV file.

        This function loads the oligo database from a tab-separated values (TSV) file. The file must contain
        columns such as 'region_id', 'oligo_id', 'sequence', and additional attributes, like information from the
        fasta headers or information computed by the filtering classes.
        The order of columns in the database file is:

        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+

        The database can be optionally filtered by specifying a list of region IDs.

        :param file_database: Path to the TSV file containing the oligo database.
        :type file_database: str
        :param region_ids: List of region IDs to filter the database. Defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :param database_overwrite: If True, overwrite the existing database. Defaults to False.
        :type database_overwrite: bool, optional

        :raises ValueError: If the database file has an incorrect format or does not exist.
        """
        region_ids = check_if_list(region_ids)

        if database_overwrite:
            warnings.warn("Overwriting database!")

        if os.path.exists(file_database):
            if not check_tsv_format(file_database):
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        file_tsv_content = pd.read_table(file_database, sep="\t")

        file_tsv_content = file_tsv_content.apply(
            lambda col: col.apply(lambda x: x if pd.notna(x) else [None])
        )

        # convert lists represented as string to proper list format in the table with the eval function
        file_tsv_content = file_tsv_content.apply(
            lambda col: col.apply(
                lambda x: (
                    eval(x)
                    if isinstance(x, str) and x.startswith("[") and x.endswith("]")
                    else ([int(x)] if isinstance(x, str) and x.isdigit() else x)
                )
            )
        )

        database_tmp1 = file_tsv_content.to_dict(orient="records")
        database_tmp2 = defaultdict(dict)
        for entry in database_tmp1:
            region_id, oligo_id = entry.pop("region_id"), entry.pop("oligo_id")
            database_tmp2[region_id][oligo_id] = entry

        database_tmp2 = dict(database_tmp2)

        if not database_overwrite and self.database:
            database_tmp2 = merge_databases(self.database, database_tmp2)

        if region_ids:
            database_tmp2 = filter_dabase_for_region(database_tmp2, region_ids)
            check_if_region_in_database(
                database_tmp2,
                region_ids,
                self.write_regions_with_insufficient_oligos,
                self.file_removed_regions,
            )

        self.database = database_tmp2

    def load_sequences_from_fasta(
        self,
        files_fasta: list[str],
        sequence_type: _TYPES_SEQ,
        region_ids: list[str] = None,
        database_overwrite: bool = False,
    ) -> None:
        """Load "oligo" or "target" sequences from one or more FASTA files into the oligo database.

        This function reads sequences from FASTA file(s) and adds them to the oligo database, eitehr as 'oligo' or
        'target' sequence and computes the reverse compliment sequence for the other sequence type. It parses the
        headers of the FASTA entries to extract region information, and assigns unique IDs to the oligos within
        each region.

        :param files_fasta: Path to the FASTA file(s) containing the sequences.
        :type files_fasta: list[str]
        :param sequence_type: Type of sequence to load, either 'target' or 'oligo'.
        :type sequence_type: _TYPES_SEQ
        :param region_ids: List of region IDs to filter the database. Defaults to None.
        :type region_ids: list[str], optional
        :param database_overwrite: If True, overwrite the existing database. Defaults to False.
        :type database_overwrite: bool, optional

        :raises AssertionError: If the provided sequence type is not supported.
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."
        sequence_type_reverse_complement = options[0] if options[0] != sequence_type else options[1]

        if database_overwrite:
            warnings.warn("Overwriting database!")
            self.database = {}

        region_ids = check_if_list(region_ids)
        files_fasta = check_if_list(files_fasta)

        for file in files_fasta:

            self.fasta_parser.check_fasta_format(file)
            fasta_sequences = self.fasta_parser.read_fasta_sequences(file, region_ids)
            region_sequences = {}
            for entry in fasta_sequences:
                region, additional_info, coordinates = self.fasta_parser.parse_fasta_header(entry.id)
                oligo_info = coordinates | additional_info
                if region in region_sequences:
                    if entry.seq in region_sequences[region]:
                        oligo_info_merged = collapse_info_for_duplicated_sequences(
                            region_sequences[region][entry.seq], oligo_info
                        )
                        region_sequences[region][str(entry.seq)] = oligo_info_merged
                    else:
                        region_sequences[region][str(entry.seq)] = oligo_info
                else:
                    region_sequences[region] = {str(entry.seq): oligo_info}

            database_loaded = {region: {} for region in region_sequences.keys()}
            for region, value in region_sequences.items():
                i = 1
                for oligo_sequence, oligo_info in value.items():
                    oligo_id = f"{region}{SEPARATOR_OLIGO_ID}{i}"
                    oligo_sequence_reverse_complement = str(Seq(oligo_sequence).reverse_complement())
                    oligo_seq_info = {
                        sequence_type: oligo_sequence,
                        sequence_type_reverse_complement: oligo_sequence_reverse_complement,
                    } | oligo_info
                    database_loaded[region][oligo_id] = oligo_seq_info
                    i += 1
            if self.database:
                self.database = merge_databases(self.database, database_loaded)
            else:
                self.database = database_loaded

        # add this step to log regions which are not available in database
        if region_ids:
            check_if_region_in_database(
                self.database,
                region_ids,
                self.write_regions_with_insufficient_oligos,
                self.file_removed_regions,
            )

        self.database = self.database

    ############################################
    # Save Functions
    ############################################

    def save_database(
        self,
        region_ids: list[str] = None,
        filename: str = "oligo_database",
    ):
        """Save the oligo database to YAML and TSV files.

        This function saves the oligo database to two files: a YAML file containing metadata and a TSV (tab-separated values)
        file containing the oligo database entries. The files are saved in the specified output directory with the provided
        filenames. The order of the columns in the TSV file is:

        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+

        Additional feat. includes additional information from the header of the fasta file and additional information
        computed by the filtering classes.

        :param region_ids: A list of region IDs to include in the saved database. If None, all regions are included.
        :type region_ids: list[str], optional
        :param filename: The base filename for the output files (without extensions), defaults to "oligo_database".
        :type filename: str, optional
        :return: Paths to the saved YAML and TSV files.
        :rtype: tuple
        """
        if region_ids:
            region_ids = check_if_list(region_ids)
        else:
            region_ids = self.database.keys()

        file_metadata = os.path.join(self.dir_output, filename + ".yaml")

        with open(file_metadata, "w") as handle:
            yaml.safe_dump(self.metadata, handle, sort_keys=True, default_flow_style=False)

        file_database = os.path.join(self.dir_output, filename + ".tsv")
        file_tsv_content = []

        for region_id, oligo_dict in self.database.items():
            if region_id in region_ids:
                for oligo_id, oligo_attributes in oligo_dict.items():
                    entry = {"region_id": region_id, "oligo_id": oligo_id}
                    entry.update(oligo_attributes)
                    file_tsv_content.append(entry)

        file_tsv_content = pd.DataFrame(data=file_tsv_content)
        file_tsv_content.to_csv(file_database, sep="\t", index=False)

        return file_database, file_metadata

    def write_database_to_fasta(
        self,
        filename: str = "oligo_database",
        region_ids: list[str] = None,
        sequence_type: _TYPES_SEQ = "oligo",
        save_description: bool = False,
    ):
        """Write oligo sequences from the database to a FASTA file.

        This function writes the sequences from the oligo database to a FASTA file. Each sequence is represented as a
        SeqRecord in the FASTA file.

        :param filename: The base filename for the output FASTA file (without extension), defaults to "oligo_database".
        :type filename: str, optional
        :param region_ids: A list of region IDs to include in the saved database. If None, all regions are included.
        :type region_ids: list[str], optional
        :param sequence_type: The type of sequence to write (e.g., "oligo" or "target").
        :type sequence_type: str
        :param save_description: Flag to include sequence descriptions in the FASTA file.
        :type save_description: bool, optional
        :return: Path to the generated FASTA file.
        :rtype: str
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        if region_ids:
            region_ids = check_if_list(region_ids)
        else:
            region_ids = self.database.keys()

        file_fasta = os.path.join(self.dir_output, f"{filename}.fna")
        output_fasta = []

        with open(file_fasta, "w") as handle_fasta:
            for region_id, oligo in self.database.items():
                if region_id in region_ids:
                    for oligo_id, oligo_attributes in oligo.items():
                        description = sequence_type if save_description else ""
                        seq_record = SeqRecord(
                            Seq(oligo_attributes[sequence_type]),
                            id=oligo_id,
                            name=oligo_id.split(SEPARATOR_OLIGO_ID)[0],
                            description=description,
                        )
                        output_fasta.append(seq_record)

            SeqIO.write(output_fasta, handle_fasta, "fasta")

        return file_fasta

    def write_oligosets(self, foldername_out: str = "oligo_sets"):
        """Write oligo sets to individual TSV files.

        This function writes the oligo sets to individual TSV files, with each file representing the oligo sets
        for a specific region.

        :param foldername_out: The name of the folder to store the oligo set files, defaults to "oligo_sets".
        :type foldername_out: str, optional
        :return: Path to the folder containing the generated oligo set files.
        :rtype: str
        """
        dir_oligosets = os.path.join(self.dir_output, foldername_out)
        Path(dir_oligosets).mkdir(parents=True, exist_ok=True)

        for region_id in self.oligosets.keys():
            file_oligosets = os.path.join(dir_oligosets, f"{region_id}_oligosets.tsv")
            self.oligosets[region_id].to_csv(file_oligosets, sep="\t", index=False)

        return dir_oligosets

    def remove_regions_with_insufficient_oligos(self, pipeline_step: str):
        """Remove regions with insufficient oligos from the oligo database.

        This function identifies regions in the oligo database that have fewer oligos than the specified
        minimum threshold (`min_oligos_per_region`). It removes those regions from the database and updates
        the associated oligo sets. If the option to write removed regions to a file is enabled, it logs the
        removed regions along with the pipeline step.

        :param pipeline_step: Step in the pipeline that led to the removal of regions.
        :type pipeline_step: str
        """
        regions_to_remove = [
            region for region, oligos in self.database.items() if len(oligos) <= self.min_oligos_per_region
        ]

        for region in regions_to_remove:
            del self.database[region]
            self.oligosets.pop(region, None)

        if self.write_regions_with_insufficient_oligos and regions_to_remove:
            with open(self.file_removed_regions, "a") as handle:
                handle.write("\n".join(f"{region}\t{pipeline_step}" for region in regions_to_remove) + "\n")

    ############################################
    # Getter Functions
    ############################################

    def get_sequence_list(self, sequence_type: _TYPES_SEQ = "oligo"):
        """Retrieve a list of sequences of the specified type (e.g., 'oligo' or 'target') from the oligo database.

        :param sequence_type: Type of sequences to retrieve (default is 'oligo').
        :type sequence_type: _TYPES_SEQ
        :return: List of sequences.
        :rtype: List[str]
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."
        sequences = [
            str(oligo_attributes[sequence_type])
            for region_id, oligo_dict in self.database.items()
            for oligo_id, oligo_attributes in oligo_dict.items()
        ]

        return sequences

    def get_sequence_oligoid_mapping(
        self, sequence_type: _TYPES_SEQ = "oligo", sequence_to_upper: bool = False
    ):
        """Generate a mapping between sequences and their corresponding oligonucleotide IDs, with an option to convert sequences to uppercase.

        Validates the sequence type against predefined options. If `sequence_to_upper` is True, converts all sequences to uppercase before mapping,
        ensuring case-insensitive comparisons. This function supports scenarios where multiple oligos share the same sequence, grouping their IDs in a list.

        :param sequence_type: The type of sequence to use for the mapping, defaulting to "oligo".
        :type sequence_type: _TYPES_SEQ
        :param sequence_to_upper: Flag indicating whether to convert sequences to uppercase.
        :type sequence_to_upper: bool
        :return: A dictionary mapping sequences to lists of corresponding oligonucleotide IDs.
        :rtype: dict
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        sequence_oligoids_mapping = {}

        for region_id, database_region in self.database.items():
            for oligo_id, oligo_attributes in database_region.items():
                seq = oligo_attributes[sequence_type]
                if sequence_to_upper:
                    seq = seq.upper()
                if seq not in sequence_oligoids_mapping:
                    # If the sequence key doesn't exist, create a new entry with the oligo ID in a list
                    sequence_oligoids_mapping[seq] = [oligo_id]
                else:
                    # If the sequence key already exists, append the oligo ID to the existing list
                    sequence_oligoids_mapping[seq].append(oligo_id)

        return sequence_oligoids_mapping

    def get_oligo_attribute(self, attribute: str):
        """Retrieves a specified attribute for all oligos in the database and returns it as a pandas DataFrame.
        This method assumes the presence of an attribute across all oligo records in the database. If the
        attribute is not found, a KeyError is raised.

        :param attribute: The name of the attribute to retrieve for each oligo.
        :type attribute: str
        :return: A pandas DataFrame with two columns: 'oligo_id' and the specified 'attribute', where each row
                corresponds to an oligo and its attribute value.
        :rtype: pd.DataFrame
        :raises KeyError: If the specified attribute has not been computed and added to the database.
        """
        if not check_if_key_exists(self.database, attribute):
            raise KeyError(f"The {attribute} attribute has not been computed!")
        oligo_ids = [
            oligo_id for region_id, oligo_dict in self.database.items() for oligo_id in oligo_dict.keys()
        ]
        attributes = [
            oligo_attributes[attribute]
            for region_id, oligo_dict in self.database.items()
            for oligo_id, oligo_attributes in oligo_dict.items()
        ]
        return pd.DataFrame({"oligo_id": oligo_ids, attribute: attributes})

    ############################################
    # Update Functions
    ############################################

    def update_oligo_attributes(self, new_oligo_attribute: dict):
        """Update attributes of oligonucleotides in the database with new values.

        Iterates through each oligonucleotide in the database and updates its attributes based on the provided new attribute values.
        The update is done in place, modifying the existing attributes of the oligonucleotides.

        Input format of new attribute:
        new_attribute = {"region_x::1": {"attribute_y": 110}}

        :param new_oligo_attribute: A dictionary mapping oligonucleotide IDs to their new attributes.
        :type new_oligo_attribute: dict
        """
        for region_id, database_region in self.database.items():
            for oligo_id, oligo_attributes in database_region.items():
                oligo_attributes.update(new_oligo_attribute[oligo_id])