from typing import Annotated

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, PositiveInt, field_validator, model_validator

from oligo_designer_toolsuite.validation._types import (
    FilesFastaReferenceDatabaseT,
    GCContentMaxT,
    GCContentMinT,
)
from oligo_designer_toolsuite.validation.models._general import (
    BaseProbabilities,
    HomopolymerThresholds,
    OligoPropertyWeights,
)


class InitiatorProbeHCR(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_initiator_table: Annotated[
        str,
        Field(
            description="Path to a CSV/TSV file containing the initiator table. The file must include columns: 'bit', 'initiator_L_sequence', and 'initiator_R_sequence'.",
        ),
    ]
    file_codebook: Annotated[
        str,
        Field(
            description="Path to a CSV/TSV file containing an existing codebook. The codebook must have region IDs as the index, columns named 'bit_*' with 0/1 entries, and exactly one active bit per row.",
        ),
    ]

    @field_validator("file_initiator_table", "file_codebook")
    @classmethod
    def must_be_csv_or_tsv(cls, v: str | None) -> str | None:
        if v is None:
            return v
        if not (v.endswith(".csv") or v.endswith(".tsv")):
            raise ValueError("File must end with .csv or .tsv")
        return v

    @field_validator("file_codebook")
    @classmethod
    def validate_codebook_content(cls, v: str) -> str:
        try:
            codebook = pd.read_csv(v, sep=None, engine="python", index_col="region_id")
        except Exception as e:
            raise ValueError(f"Could not read codebook file '{v}': {e}")
        if len(codebook.columns) == 0:
            raise ValueError(f"Codebook file '{v}' must contain at least one column.")
        non_bit_columns = [col for col in codebook.columns if not str(col).startswith("bit_")]
        if non_bit_columns:
            raise ValueError(
                f"Codebook file '{v}' must have all columns named 'bit_*'. "
                f"Found non-matching columns: {non_bit_columns}"
            )
        codebook = codebook[codebook.notna().all(axis=1)]
        codebook = codebook[codebook.isin([0, 1]).all(axis=1)]
        codebook = codebook[(codebook == 1).sum(axis=1) == 1]
        if len(codebook) == 0:
            raise ValueError(f"Codebook file '{v}' must contain at least one valid one-hot encoded row.")
        return v

    @field_validator("file_initiator_table")
    @classmethod
    def validate_initiator_table_content(cls, v: str) -> str:
        required_cols = {"bit", "initiator_L_sequence", "initiator_R_sequence"}
        try:
            initiator_table = pd.read_csv(v, sep=None, engine="python")
        except Exception as e:
            raise ValueError(f"Could not read initiator table file '{v}': {e}")
        missing = required_cols - set(initiator_table.columns)
        if missing:
            raise ValueError(
                f"Initiator table '{v}' is missing required columns: {sorted(missing)}. "
                f"Required columns are: {sorted(required_cols)}."
            )
        return v

    @model_validator(mode="after")
    def validate_initiator_bits_in_codebook(self) -> "InitiatorProbeHCR":
        try:
            codebook = pd.read_csv(self.file_codebook, sep=None, engine="python", index_col="region_id")
            initiator_table = pd.read_csv(self.file_initiator_table, sep=None, engine="python")
        except Exception:
            return self
        codebook_bits = set(codebook.columns)
        initiator_bits = set(initiator_table["bit"].dropna().unique())
        missing_bits = initiator_bits ^ codebook_bits
        if missing_bits:
            raise ValueError(
                f"Initiator table references bit columns mismatch with codebook bits: {sorted(missing_bits)}."
            )
        return self


class ReadoutProbeCycleHCR(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_readout_probe_table: Annotated[
        str,
        Field(
            description="Path to a CSV/TSV file containing the readout probe table. The file must include columns: 'channel', 'readout_probe_id', 'L/R', and 'readout_probe_sequence'. If a 'bit' column is not present, it will be automatically assigned.",
        ),
    ]
    file_codebook: Annotated[
        str | None,
        Field(
            description="Path to a CSV/TSV file containing an existing codebook, or None to generate a new codebook. If provided, the codebook must have region IDs as the index and bit columns (with 0/1) named 'bit_1', 'bit_2', etc. If None, a codebook will be automatically generated based on the number of regions and available readout probes.",
        ),
    ]

    @field_validator("file_readout_probe_table", "file_codebook")
    @classmethod
    def must_be_csv_or_tsv(cls, v: str | None) -> str | None:
        if v is None:
            return v
        if not (v.endswith(".csv") or v.endswith(".tsv")):
            raise ValueError("File must end with .csv or .tsv")
        return v


class ReadoutProbeFish(BaseModel):
    model_config = ConfigDict(extra="forbid")

    files_fasta_reference_database: Annotated[
        FilesFastaReferenceDatabaseT,
        Field(
            description="List of paths to FASTA files containing reference sequences used for specificity filtering. These files are used to identify off-target binding sites and potential cross-hybridization events."
        ),
    ]
    length: Annotated[PositiveInt, Field(description="Length of readout probes.")]
    base_probabilities: Annotated[
        BaseProbabilities,
        Field(
            description="Probabilities of each nucleotide base in randomly generated readout probe sequences. Keys should be 'A', 'T', 'G', 'C' and values should sum to 1.0."
        ),
    ]
    GC_content_min: GCContentMinT
    GC_content_max: GCContentMaxT
    homopolymeric_base_n: Annotated[
        HomopolymerThresholds,
        Field(
            description="The minimum number of nucleotides to consider it a homopolymeric run per base",
        ),
    ]
    channels_ids: Annotated[
        list[str],
        Field(
            description="List of fluorescence channel identifiers (e.g., ['Cy3', 'Cy5', 'Alexa488']) to which readout probes will be assigned. Readout probes are distributed across channels in a round-robin fashion."
        ),
    ]


class ReadoutProbeMerfish(ReadoutProbeFish):

    set_size: Annotated[PositiveInt, Field(description="total number of readout probes")]
    homogeneous_properties_weights: Annotated[
        OligoPropertyWeights,
        Field(
            description="Specifying weights for property homogeneity in set selection. Readout probes in one set should have similar properties, the weighted sum of variance of the properties is minimized. The values are weights that determine the relative importance of each property in the homogeneity score."
        ),
    ]
    n_bits: Annotated[
        PositiveInt,
        Field(
            description="Number of bits in each barcode in the codebook. This determines the maximum number of unique barcodes that can be generated."
        ),
    ]
    min_hamming_dist: Annotated[
        PositiveInt,
        Field(
            description="Minimum Hamming distance required between any two barcodes in the codebook. Higher values provide better error correction but reduce the number of available barcodes."
        ),
    ]
    hamming_weight: Annotated[
        PositiveInt,
        Field(
            description="Fixed Hamming weight (number of active bits, value 1) for each barcode. All barcodes will have exactly this many bits set to 1."
        ),
    ]


class ReadoutProbeSeqFishPlus(ReadoutProbeFish):

    n_barcode_rounds: Annotated[
        PositiveInt,
        Field(
            description="Number of barcode rounds in the encoding scheme. This determines how many rounds of barcoding are used, and each round requires readout probes."
        ),
    ]
    n_pseudocolors: Annotated[
        PositiveInt,
        Field(
            description="Number of pseudocolors to use in the barcoding scheme. Pseudocolors expand the barcode space by allowing multiple states per round. The final round includes a checksum pseudocolor."
        ),
    ]
