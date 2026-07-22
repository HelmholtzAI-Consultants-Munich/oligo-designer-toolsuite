"""
Genomic region generator pipeline.

Probe design pipelines need FASTA files of transcript or genomic regions to
build target libraries and to check specificity. This module prepares those
files from an annotation source such as NCBI, Ensembl, or user-provided GTF and
genome FASTA files.

See :class:`GenomicRegionGenerator` for the callable API. See :func:`main` for
the config-driven command-line workflow.
"""

############################################
# imports
############################################

import inspect
import os
from pathlib import Path

import yaml

from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.pipelines._utils import base_log_parameters, base_parser
from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
)
from oligo_designer_toolsuite.utils import configure_root_logger, logger

############################################
# Genomic Region Generator Functions
############################################


class GenomicRegionGenerator:
    """
    Produce per-region FASTA files from a genomic annotation source.

    This class prepares the reference sequences used by the probe design
    pipelines. It loads annotations from NCBI, Ensembl, or custom GTF and genome
    FASTA files, then writes separate FASTA files for the requested region types.

    Those files can be used later as target sequences or as specificity
    references. Common region types include genes, exons, introns, coding
    sequence, UTRs, exon-exon junctions, and intergenic regions.

    Typical workflow:

    1. Load annotations with :py:meth:`load_annotations`.
    2. Generate the selected region FASTA files with
       :py:meth:`generate_genomic_regions`.

    :param dir_output: Directory where downloaded annotations and per-region
        FASTA files are written. Created if it does not exist.
    :type dir_output: str
    """

    def __init__(self, dir_output: str) -> None:
        """Constructor for the GenomicRegionGenerator class."""
        self.dir_output = os.path.abspath(dir_output)
        Path(dir_output).mkdir(parents=True, exist_ok=True)

    def load_annotations(
        self,
        source: str,
        source_params: dict,
    ) -> CustomGenomicRegionGenerator:
        """
        Load a genomic annotation source for region generation.

        This step prepares a region generator from NCBI, Ensembl, or custom
        annotation files. For NCBI and Ensembl, the required GTF and genome FASTA
        files are downloaded when needed. For a custom source, existing local
        files are used.

        The returned generator is ready to slice region sequences. Pass it to
        :py:meth:`generate_genomic_regions` to write the FASTA files.

        :param source: Annotation source. One of ``"ncbi"``, ``"ensembl"``, or
            ``"custom"``.
        :type source: str
        :param source_params: Source-specific settings from the ``source_params``
            section of the config.

            - ``"ncbi"``: species-based download (``taxon``, ``species``,
              ``annotation_release``) or assembly-based download
              (``refseq_assembly_accession`` and/or ``assembly_name``). The
              ``mode`` setting selects which path is used.
            - ``"ensembl"``: ``species`` and ``annotation_release``.
            - ``"custom"``: paths to the annotation GTF and genome FASTA, plus
              metadata such as ``files_source``, ``species``,
              ``annotation_release``, and ``genome_assembly``.

        :type source_params: dict
        :return: Region generator ready to produce per-region FASTA files.
        :rtype: CustomGenomicRegionGenerator
        :raises ConfigurationError: If ``source`` is unknown or a required NCBI
            setting is missing.
        """
        logger.info("Parameters Load Annotations:")
        frame = inspect.currentframe()
        if frame is not None:
            args, _, _, values = inspect.getargvalues(frame)
            parameters = {i: values[i] for i in args}
            base_log_parameters(parameters)

        region_generator: (
            CustomGenomicRegionGenerator | NcbiGenomicRegionGenerator | EnsemblGenomicRegionGenerator | None
        ) = None

        if source == "ncbi":
            if source_params["mode"] is None:
                raise ConfigurationError(
                    "For source='ncbi', source_params parameter 'mode' must be provided."
                )
            region_generator = NcbiGenomicRegionGenerator(
                mode=source_params["mode"],
                taxon=source_params["taxon"],
                species=source_params["species"],
                annotation_release=source_params["annotation_release"],
                assembly_source=source_params["assembly_source"],
                refseq_assembly_accession=source_params["refseq_assembly_accession"],
                assembly_name=source_params["assembly_name"],
                dir_output=self.dir_output,
            )
        elif source == "ensembl":
            region_generator = EnsemblGenomicRegionGenerator(
                species=source_params["species"],
                annotation_release=source_params["annotation_release"],
                dir_output=self.dir_output,
            )
        elif source == "custom":
            region_generator = CustomGenomicRegionGenerator(
                annotation_file=source_params["file_annotation"],
                sequence_file=source_params["file_sequence"],
                files_source=source_params["files_source"],
                species=source_params["species"],
                annotation_release=source_params["annotation_release"],
                genome_assembly=source_params["genome_assembly"],
                dir_output=self.dir_output,
            )
        else:
            raise ConfigurationError(
                f"Source '{source}' is not supported. Supported sources are: 'NCBI', 'Ensembl', or 'custom'."
            )

        logger.info(
            f"The following annotation files are used for GTF annotation of regions: {region_generator.annotation_file} and for fasta sequence file: {region_generator.sequence_file} ."
        )
        logger.info(
            f"The annotations are from {region_generator.files_source} source, for the species: {region_generator.species}, release number: {region_generator.annotation_release} and genome assembly: {region_generator.genome_assembly}"
        )
        return region_generator

    def generate_genomic_regions(
        self,
        region_generator: CustomGenomicRegionGenerator,
        genomic_regions: dict,
        block_size: int = 50,
    ) -> list:
        """
        Write FASTA files for the requested genomic region types.

        Each enabled region type is written as its own FASTA file under the
        output directory. These files are the sequences later used by probe
        design pipelines as targets or specificity references.

        Supported region types and typical uses:

        - ``gene``: full gene sequences.
        - ``exon``: individual exons.
        - ``intron``: intronic sequences.
        - ``cds``: coding sequence only.
        - ``utr``: untranslated regions.
        - ``exon_exon_junction``: short fragments spanning splice sites, useful
          for isoform-specific probes.
        - ``intergenic``: sequences between genes, often used as an off-target
          reference.

        :param region_generator: Region generator returned by
            :py:meth:`load_annotations`.
        :type region_generator: CustomGenomicRegionGenerator
        :param genomic_regions: Map of region type to an enabled flag. Region
            types with a falsy flag are skipped. Recognised types are ``"gene"``,
            ``"intergenic"``, ``"exon"``, ``"intron"``, ``"cds"``, ``"utr"``, and
            ``"exon_exon_junction"``.
        :type genomic_regions: dict
        :param block_size: Half-length of each exon-exon junction fragment
            centred on the splice site, in bases. Used only for
            ``"exon_exon_junction"``. Defaults to 50.
        :type block_size: int
        :return: Absolute paths to the generated FASTA files, one per enabled
            region type.
        :rtype: list
        :raises ConfigurationError: If ``genomic_regions`` contains an
            unrecognised region type.
        """
        files_fasta = []
        for genomic_region, flag in genomic_regions.items():
            if flag:
                if genomic_region == "gene":
                    file_fasta = region_generator.get_sequence_gene()
                elif genomic_region == "intergenic":
                    file_fasta = region_generator.get_sequence_intergenic()
                elif genomic_region == "exon":
                    file_fasta = region_generator.get_sequence_exon()
                elif genomic_region == "intron":
                    file_fasta = region_generator.get_sequence_intron()
                elif genomic_region == "cds":
                    file_fasta = region_generator.get_sequence_CDS()
                elif genomic_region == "utr":
                    file_fasta = region_generator.get_sequence_UTR()
                elif genomic_region == "exon_exon_junction":
                    file_fasta = region_generator.get_sequence_exon_exon_junction(block_size=block_size)
                else:
                    raise ConfigurationError(
                        f"Genomic region type '{genomic_region}' is not supported. "
                        f"Supported types are: 'gene', 'intergenic', 'exon', 'intron', 'cds', 'utr', 'exon_exon_junction'."
                    )

                files_fasta.append(file_fasta)
                logger.info(f"The genomic region '{genomic_region}' was stored in :{file_fasta}.")

        return files_fasta


############################################
# Genomic Region Generator Pipeline
############################################


def main() -> None:
    """
    Run the genomic region generator pipeline from the command line.

    Parses the required ``-c``/``--config`` argument, loads the YAML
    configuration file, and configures the library logger to write under the
    configured output directory. It then loads the annotation source and writes
    the selected region FASTA files with :class:`GenomicRegionGenerator`.

    The config should follow the YAML files under ``data/configs/``, for example
    ``data/configs/genomic_region_generator_ncbi.yaml``.

    Top-level config sections:

    - ``dir_output``: directory for downloaded annotations and produced FASTA
      files.
    - ``source``: annotation source (``"ncbi"``, ``"ensembl"``, or ``"custom"``).
    - ``source_params``: source-specific settings used by
      :py:meth:`GenomicRegionGenerator.load_annotations`.
    - ``genomic_regions``: map of region type to an enabled flag.
    - ``exon_exon_junction_block_size``: half-length of exon-exon junction
      fragments.

    Files written under ``dir_output``:

    - one FASTA file for each enabled genomic region type.
    - downloaded annotation and genome files when NCBI or Ensembl is used.

    See :class:`GenomicRegionGenerator` for the pipeline description.

    :return: None
    :rtype: None
    """
    print("--------------START PIPELINE--------------")
    args = base_parser()

    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    pipeline = GenomicRegionGenerator(dir_output=config["dir_output"])

    configure_root_logger(
        dir_output=pipeline.dir_output,
        pipeline_name="genomic_region_generation",
        include_console=True,
    )

    region_generator = pipeline.load_annotations(
        source=config["source"],
        source_params=config["source_params"],
    )

    files_fasta = pipeline.generate_genomic_regions(
        region_generator=region_generator,
        genomic_regions=config["genomic_regions"],
        block_size=config["exon_exon_junction_block_size"],
    )

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
