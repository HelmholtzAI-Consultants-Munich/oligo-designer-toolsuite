from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, PositiveInt

from oligo_designer_toolsuite.validation._types import (
    SecondaryStructuresThresholdDeltaGT,
)
from oligo_designer_toolsuite.validation.models._general import (
    BlastnHitParameters,
    BlastnSearchParameters,
    CrossHybridizationProbabilityFilterConfig,
    HybridizationProbabilityFilterConfig,
    OligoSetSelection,
    TmChemCorrectionParameters,
    TmParameters,
    TmSaltCorrectionParameters,
)


class TargetProbeDev(BaseModel):
    model_config = ConfigDict(extra="forbid")

    secondary_structures_threshold_deltaG: SecondaryStructuresThresholdDeltaGT
    specificity_blastn_search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="BLASTN search parameters for specificity filtering. These parameters control how BLASTN searches are performed to identify off-target binding sites."
        ),
    ]
    specificity_blastn_hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits during specificity analysis. Use either coverage or min_alignment_length."
        ),
    ]
    cross_hybridization_blastn_search_parameters: Annotated[
        BlastnSearchParameters,
        Field(description="Parameters for BLASTN searches used in cross-hybridization filtering."),
    ]
    cross_hybridization_blastn_hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits in cross-hybridization searches. Use either coverage or min_alignment_length."
        ),
    ]


class TargetProbeDevCycleHCR(TargetProbeDev):

    Tm_parameters: Annotated[
        TmParameters,
        Field(
            description="Parameters for calculating melting temperature (Tm) using the nearest-neighbor method. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN"
        ),
    ]
    Tm_chem_correction_parameters: Annotated[
        TmChemCorrectionParameters,
        Field(
            description="Optional parameters for chemical correction.  For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction"
        ),
    ]
    Tm_salt_correction_parameters: Annotated[
        TmSaltCorrectionParameters,
        Field(
            description="Optional parameters to account for the effects of salt concentration on melting temperature. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction"
        ),
    ]

    cross_hybridization_blastn_search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="Pydantic model of BLASTN search parameters for cross-hybridization filtering. These parameters control how BLASTN searches are performed to identify potential cross-hybridization between left and right probe pairs within the same set."
        ),
    ]
    cross_hybridization_blastn_hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits in cross-hybridization searches. Use either coverage or min_alignment_length. Probes with cross-hybridization hits meeting these criteria are removed from the larger region."
        ),
    ]


class DeveloperParametersBase(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_set_selection: OligoSetSelection


class DeveloperParametersCycleHCR(DeveloperParametersBase):
    target_probe: TargetProbeDevCycleHCR


class TargetProbeDevMerfish(TargetProbeDev):
    Tm_parameters: Annotated[
        TmParameters,
        Field(
            description="Parameters for calculating melting temperature (Tm) using the nearest-neighbor method. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN"
        ),
    ]
    Tm_chem_correction_parameters: Annotated[
        TmChemCorrectionParameters,
        Field(
            description="Optional parameters for chemical correction.  For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction"
        ),
    ]
    Tm_salt_correction_parameters: Annotated[
        TmSaltCorrectionParameters,
        Field(
            description="Optional parameters to account for the effects of salt concentration on melting temperature. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction"
        ),
    ]


class ReadoutProbeDevFish(BaseModel):
    model_config = ConfigDict(extra="forbid")

    initial_num_sequences: Annotated[
        PositiveInt,
        Field(
            description="Number of random sequences to generate initially before filtering. Higher values provide more candidates but increase computation time. If not enough readout probes can be generated, increase this number"
        ),
    ]
    specificity_blastn_search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="BLASTN search parameters for specificity filtering. These parameters control how BLASTN searches are performed to identify off-target binding sites."
        ),
    ]
    specificity_blastn_hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits during specificity analysis. Use either coverage or min_alignment_length."
        ),
    ]
    cross_hybridization_blastn_search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="Pydantic model of BLASTN search parameters for cross-hybridization filtering. These parameters control how BLASTN searches are performed to identify potential cross-hybridization between left and right probe pairs within the same set."
        ),
    ]
    cross_hybridization_blastn_hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits in cross-hybridization searches. Use either coverage or min_alignment_length. Probes with cross-hybridization hits meeting these criteria are removed from the larger region."
        ),
    ]


class ReadoutProbeDevMerfish(ReadoutProbeDevFish):
    Tm_parameters: Annotated[
        TmParameters,
        Field(
            description="Parameters for calculating melting temperature (Tm) using the nearest-neighbor method. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN"
        ),
    ]
    Tm_chem_correction_parameters: Annotated[
        TmChemCorrectionParameters,
        Field(
            description="Optional parameters for chemical correction.  For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction"
        ),
    ]
    Tm_salt_correction_parameters: Annotated[
        TmSaltCorrectionParameters,
        Field(
            description="Optional parameters to account for the effects of salt concentration on melting temperature. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction"
        ),
    ]
    n_combinations: Annotated[
        PositiveInt,
        Field(
            description="Number of combinations to evaluate during set creation. Higher values may find better sets but increase computation time.",
        ),
    ]


class PrimerDevFish(BaseModel):
    model_config = ConfigDict(extra="forbid")

    initial_num_sequences: Annotated[
        PositiveInt,
        Field(
            description="Number of random sequences to generate initially before filtering. Higher values provide more candidates but increase computation time. If not enough readout probes can be generated, increase this number"
        ),
    ]
    Tm_parameters: Annotated[
        TmParameters,
        Field(
            description="Parameters for calculating melting temperature (Tm) using the nearest-neighbor method. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN"
        ),
    ]
    Tm_chem_correction_parameters: Annotated[
        TmChemCorrectionParameters,
        Field(
            description="Optional parameters for chemical correction.  For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction"
        ),
    ]
    Tm_salt_correction_parameters: Annotated[
        TmSaltCorrectionParameters,
        Field(
            description="Optional parameters to account for the effects of salt concentration on melting temperature. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction"
        ),
    ]
    specificity_reference_blastn_search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="BLASTN search parameters for specificity filtering. These parameters control how BLASTN searches are performed to identify off-target binding sites."
        ),
    ]
    specificity_reference_blastn_hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits during specificity analysis. Use either coverage or min_alignment_length."
        ),
    ]
    specificity_hybridization_probes_blastn_search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="BLASTN search parameters for specificity filtering against the hybridization probes database.."
        ),
    ]
    specificity_hybridization_probes_blastn_hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits during specificity analysis against the hybridization probes database. Use either coverage or min_alignment_length."
        ),
    ]


class PrimerDevSeqFishPlus(PrimerDevFish):
    secondary_structures_threshold_deltaG: SecondaryStructuresThresholdDeltaGT


class DeveloperParametersMerfish(DeveloperParametersBase):
    target_probe: TargetProbeDevMerfish
    readout_probe: ReadoutProbeDevMerfish
    primer: PrimerDevFish


class DeveloperParametersSeqFishPlus(DeveloperParametersBase):
    target_probe: TargetProbeDev
    readout_probe: ReadoutProbeDevFish
    primer: PrimerDevSeqFishPlus


class TargetProbeDevOligoSeq(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: Annotated[
        TmParameters,
        Field(
            description="Parameters for calculating melting temperature (Tm) using the nearest-neighbor method. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN"
        ),
    ]
    Tm_chem_correction_parameters: Annotated[
        TmChemCorrectionParameters,
        Field(
            description="Optional parameters for chemical correction.  For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction"
        ),
    ]
    Tm_salt_correction_parameters: Annotated[
        TmSaltCorrectionParameters,
        Field(
            description="Optional parameters to account for the effects of salt concentration on melting temperature. For more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction"
        ),
    ]
    hybridization_probability: Annotated[
        HybridizationProbabilityFilterConfig,
        Field(description="Parameters for hybridization probability filtering."),
    ]
    cross_hybridization: Annotated[
        CrossHybridizationProbabilityFilterConfig,
        Field(description="Parameters for cross-hybridization filtering."),
    ]


class DeveloperParametersOligoSeq(DeveloperParametersBase):
    target_probe: TargetProbeDevOligoSeq
