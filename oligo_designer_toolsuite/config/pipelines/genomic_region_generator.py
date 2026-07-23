from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, PositiveInt, model_validator
from typing_extensions import Self

from oligo_designer_toolsuite._constants import SUPPORTED_TAXA_SOURCES


class GenomicRegions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    gene: Annotated[bool, Field(description="Generate gene regions.", default=False)]
    intergenic: Annotated[bool, Field(description="Generate intergenic regions.", default=False)]
    exon: Annotated[bool, Field(description="Generate exon regions.", default=True)]
    exon_exon_junction: Annotated[
        bool, Field(description="Generate exon–exon junction regions.", default=False)
    ]
    utr: Annotated[bool, Field(description="Generate UTR regions.", default=False)]
    cds: Annotated[bool, Field(description="Generate coding sequence (CDS) regions.", default=False)]
    intron: Annotated[bool, Field(description="Generate intron regions.", default=False)]


class AnnotationParamsCustom(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_annotation: Annotated[
        str,
        Field(
            description="Required: The path to the annotation file (GTF).",
        ),
    ]
    file_sequence: Annotated[
        str,
        Field(
            description="Required: The path to the corresponding sequence file (FASTA).",
        ),
    ]
    files_source: Annotated[
        str | None,
        Field(
            description="Optional: The original source of the files (e.g., Ensembl, NCBI), defaults to 'custom'.",
            default="custom",
        ),
    ]
    species: Annotated[
        str | None,
        Field(
            description="Optional: The species name related to the annotation and sequence files, defaults to 'unknown'.",
            default="unknown",
        ),
    ]
    annotation_release: Annotated[
        str | None,
        Field(
            description="Optional: The annotation release version, defaults to 'unknown'.", default="unknown"
        ),
    ]
    genome_assembly: Annotated[
        str | None,
        Field(description="Optional: The genome assembly version, defaults to 'unknown'.", default="unknown"),
    ]


class AnnotationParamsEnsembl(BaseModel):
    model_config = ConfigDict(extra="forbid")

    species: str = Field(
        description="The species of provided annotation.",
        default="homo_sapiens",
    )
    annotation_release: str = Field(
        description="The version of the annotation release to use, e.g. 'release-116'. Use 'current' to use the most recent annotation release. Check out release numbers for ensembl at ftp.ensembl.org/pub/",
        default="current",
    )


class AnnotationParamsNcbiSpecies(BaseModel):
    model_config = ConfigDict(extra="forbid")

    taxon: Annotated[
        Literal[
            "archaea",
            "bacteria",
            "fungi",
            "invertebrate",
            "metagenomes",
            "plant",
            "plants",
            "protozoa",
            "unknown",
            "vertebrate_mammalian",
            "vertebrate_other",
            "viral",
        ],
        Field(
            description="Taxon of the species. Supported taxa are: archaea, bacteria, fungi, invertebrate, metagenomes, plant, plants, protozoa, unknown, vertebrate_mammalian, vertebrate_other, viral.",
            default="vertebrate_mammalian",
        ),
    ]
    species: Annotated[
        str,
        Field(description="Species of provided annotation in NCBI download format.", default="Homo_sapiens"),
    ]
    assembly_source: Annotated[
        Literal["auto", "annotation_releases", "latest_assembly_versions", "reference"],
        Field(
            description="Determines how the assembly for a species is selected from the possible sources within the NCBI FTP directory. 'annotation_releases' directory, should exist for all eukaryotic species and contains assemblies annotated with different annotation versions and the annotation version can be specified by 'annotation_release'. 'latest_assembly_versions' directory is available for all species and contains the latest assembly. 'reference' directory contains the reference genome. This is only available for a subset of species. 'auto' automatically selects an assembly source in the following order (if available): 'annotation_releases', 'latest_assembly_version'",
            default="auto",
        ),
    ]
    annotation_release: Annotated[
        str,
        Field(
            description="Release number of provided annotation (e.g., '109' or 'current'); use release number only with assembly_source=annotation_releases; otherwise set to 'current'",
            default="current",
        ),
    ]

    @model_validator(mode="after")
    def _check_supported_taxa_sources_and_annotation_release(self) -> Self:
        available_sources = SUPPORTED_TAXA_SOURCES[self.taxon]
        if self.assembly_source != "auto" and self.assembly_source not in available_sources:
            supported_sources = ", ".join(sorted(available_sources))
            raise ValueError(
                f"assembly_source '{self.assembly_source}' is not available for taxon '{self.taxon}'. "
                f"Supported sources for this taxon are: {supported_sources}."
            )
        if self.annotation_release != "current" and (
            self.assembly_source in {"latest_assembly_versions", "reference"}
            or (self.assembly_source == "auto" and "annotation_releases" not in available_sources)
        ):
            raise ValueError(
                "A numeric annotation_release is only supported with assembly_source='annotation_releases'. "
                "Use annotation_release='current' for latest_assembly_versions/reference."
            )
        return self


class AnnotationNcbiSpecies(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Annotated[
        Literal["ncbi"],
        Field(
            description="Indicate which annotation should be used. Possible values are 'ensembl', 'ncbi' or 'custom'."
        ),
    ]

    mode: Annotated[
        Literal["species"],
        Field(
            description="How do you want to specify which genome you want to use? Possible values are 'species' and 'assembly'. For 'species', parameters needs the arguments 'taxon', 'species', 'assembly_source' and 'annotation_release'.  For 'assembly', parameters needs the arguments 'refseq_assembly_accession' and 'assembly_name'."
        ),
    ]

    parameters: Annotated[
        AnnotationParamsNcbiSpecies,
        Field(
            description="If a custom source, the metadata of the provided genome and annotation. If ensembl or ncbi, the parameters used to retrieve the genomic data."
        ),
    ]


class AnnotationParamsNcbiAssembly(BaseModel):
    model_config = ConfigDict(extra="forbid")

    refseq_assembly_accession: Annotated[
        str,
        Field(
            description="RefSeq assembly accession number (e.g., 'GCF_000001405.38') of the specified assembly.",
            pattern=r"^(GCF)_(\d{9})\.(\d+)$",
            default="GCF_000001405.38",
        ),
    ]
    assembly_name: Annotated[
        str, Field(description="Name (e.g., 'GRCh38.p12') of the specified assembly.", default="GRCh38.p12")
    ]


class AnnotationNcbiAssembly(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Annotated[
        Literal["ncbi"],
        Field(
            description="Indicate which annotation should be used. Possible values are 'ensembl', 'ncbi' or 'custom'."
        ),
    ]

    mode: Annotated[
        Literal["assembly"],
        Field(
            description="How do you want to specify which genome you want to use? Possible values are 'species' and 'assembly'."
        ),
    ]

    parameters: Annotated[
        AnnotationParamsNcbiAssembly,
        Field(
            description="If a custom source, the metadata of the provided genome and annotation. If ensembl or ncbi, the parameters used to retrieve the genomic data."
        ),
    ]


AnnotationNcbi = Annotated[AnnotationNcbiSpecies | AnnotationNcbiAssembly, Field(discriminator="mode")]


class AnnotationCustom(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Annotated[
        Literal["custom"],
        Field(
            description="Indicate which annotation should be used. Possible values are 'ensembl', 'ncbi' or 'custom'."
        ),
    ]
    parameters: Annotated[
        AnnotationParamsCustom,
        Field(
            description="If a custom source, the metadata of the provided genome and annotation. If ensembl or ncbi, the parameters used to retrieve the genomic data."
        ),
    ]


class AnnotationEnsembl(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Annotated[
        Literal["ensembl"],
        Field(
            description="Indicate which annotation should be used. Possible values are 'ensembl', 'ncbi' or 'custom'."
        ),
    ]
    parameters: Annotated[
        AnnotationParamsEnsembl,
        Field(
            description="If a custom source, the metadata of the provided genome and annotation. If ensembl or ncbi, the parameters used to retrieve the genomic data."
        ),
    ]


AnnotationConfigs = Annotated[
    AnnotationCustom | AnnotationEnsembl | AnnotationNcbi, Field(discriminator="source")
]


class GenomicRegionGeneratorConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = 2
    dir_output: str | None = Field(description="Name of the directory where the output files will be written")
    annotation: Annotated[AnnotationConfigs, Field(description="Parameters for genome and gene annotation.")]
    genomic_regions: Annotated[
        GenomicRegions,
        Field(
            description="List of genomic regions that should be generated, set the genomic regions you want to generate to true."
        ),
    ]
    exon_exon_junction_block_size: PositiveInt | None = Field(
        description="If exon_exon_junction is set to true, specify the block size, i.e. +/- 'block_size' bp around the junction. It does not make sense to set the block size larger than the maximum oligo length",
        default=50,
    )

    @model_validator(mode="after")
    def _default_dir_output(self) -> Self:
        # The dir_output won't be rendered in odt-cloud, because there it is not set by the user, so a validator-derived default is fine here.
        if self.dir_output is None:
            self.dir_output = f"output_genomic_region_generator_{self.annotation.source}"
        return self

    @model_validator(mode="after")
    def _check_exon_exon_junction_block_size(self) -> Self:
        if self.genomic_regions.exon_exon_junction and self.exon_exon_junction_block_size is None:
            raise ValueError(
                "exon_exon_junction_block_size must be set if exon_exon_junction is set to True."
            )
        return self
