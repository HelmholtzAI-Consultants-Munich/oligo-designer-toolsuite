############################################
# imports
############################################

import inspect
import logging
import os
from pathlib import Path

import yaml
from pydantic import ValidationError

from oligo_designer_toolsuite.pipelines._utils import (
    base_log_parameters,
    base_parser,
    setup_logging,
    write_config_to_yaml,
)
from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
)
from oligo_designer_toolsuite.validation.models._general import GenomicRegions, SourceConfigs
from oligo_designer_toolsuite.validation.models.config_pipelines import GenomicRegionGeneratorConfig

############################################
# Genomic Region Generator Functions
############################################


class GenomicRegionGenerator:
    """
    A class to generate genomic regions and manage annotations. This class allows loading of annotations from different
    sources (NCBI, Ensembl, or custom files), and generates genomic regions such as genes, intergenic regions, exons, etc.

    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(self, dir_output: str) -> None:
        """Constructor for the GenomicRegionGenerator class."""
        # create the output folder
        self.dir_output = os.path.abspath(dir_output)
        Path(dir_output).mkdir(parents=True, exist_ok=True)

        # setup logger
        setup_logging(
            dir_output=self.dir_output,
            pipeline_name="genomic_region_generation",
            include_console=True,
        )

    def load_annotations(
        self,
        source_params: SourceConfigs,
    ) -> CustomGenomicRegionGenerator:
        """
        Loads annotations from the specified source (NCBI, Ensembl, or custom files).

        :param source_params: Parameters required for loading the annotations depending on the source.
            If source is 'ncbi', it should contain 'taxon', 'species', 'annotation_release' and can optionally
            contain 'assembly_source', 'refseq_assembly_accession', and 'assembly_name'.
            If source is 'ensembl', it should contain 'species' and 'annotation_release'.
            If source is 'custom', it should contain 'file_annotation', 'file_sequence', 'files_source',
            'species', 'annotation_release', and 'genome_assembly'.
        :type source_params: SourceConfigs
        :return: An instance of the corresponding region generator class based on the source.
        :rtype: CustomGenomicRegionGenerator
        """
        ##### log parameters #####
        logging.info("Parameters Load Annotations:")
        frame = inspect.currentframe()
        if frame is not None:
            args, _, _, values = inspect.getargvalues(frame)
            parameters = {i: values[i] for i in args}
            base_log_parameters(parameters)

        region_generator: (
            CustomGenomicRegionGenerator | NcbiGenomicRegionGenerator | EnsemblGenomicRegionGenerator | None
        ) = None

        ##### loading annotations from different sources #####
        if source_params.source == "ncbi":
            # dowload the fasta files from the NCBI server
            region_generator = NcbiGenomicRegionGenerator(
                source_params=source_params,
                dir_output=self.dir_output,
            )
        elif source_params.source == "ensembl":
            # dowload the fasta files from the NCBI server
            region_generator = EnsemblGenomicRegionGenerator(
                source_params=source_params,
                dir_output=self.dir_output,
            )
        elif source_params.source == "custom":
            # use already dowloaded files
            region_generator = CustomGenomicRegionGenerator(
                source_params=source_params,
                dir_output=self.dir_output,
            )

        ##### save annotation information #####
        logging.info(
            f"The following annotation files are used for GTF annotation of regions: {region_generator.annotation_file} and for fasta sequence file: {region_generator.sequence_file} ."
        )
        logging.info(
            f"The annotations are from {region_generator.files_source} source, for the species: {region_generator.species}, release number: {region_generator.annotation_release} and genome assembly: {region_generator.genome_assembly}"
        )
        return region_generator

    def generate_genomic_regions(
        self,
        region_generator: CustomGenomicRegionGenerator,
        genomic_regions: GenomicRegions,
        block_size: int = 50,
    ) -> list:
        """
        Generates the specified genomic regions (e.g., genes, intergenic, exons, etc.) using the provided region generator.

        :param region_generator: An instance of CustomGenomicRegionGenerator that contains methods for generating
            genomic regions.
        :type region_generator: CustomGenomicRegionGenerator
        :param genomic_regions: A validated pydantic model where keys are genomic region types (e.g., 'gene', 'intergenic', etc.)
            and values are flags indicating whether to generate that region.
        :type genomic_regions: GenomicRegions
        :param block_size: Size of the block for genomic regions like exon-exon junctions. Defaults to 50.
        :type block_size: int
        :return: A list of file paths to the generated genomic regions.
        :rtype: list
        """
        files_fasta = []
        # loop not parallizeable due to file access restrictions
        for genomic_region, flag in genomic_regions:
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

                files_fasta.append(file_fasta)
                logging.info(f"The genomic region '{genomic_region}' was stored in :{file_fasta}.")

        return files_fasta


############################################
# Genomic Region Generator Pipeline
############################################


def main() -> None:
    """
    Main function to execute the genomic region generation pipeline.

    The pipeline reads a configuration file, initializes a `GenomicRegionGenerator`,
    loads annotations from the specified source, and generates genomic regions based on the provided configuration.

    :param args: Command-line arguments parsed using the base parser. The arguments include:
        - config: Path to the configuration YAML file containing parameters for the pipeline.
    :type args: dict
    """
    logging.info("--------------START PIPELINE--------------")
    args = base_parser()

    # read the config file
    with open(args["config"], "r") as handle:
        config_raw = yaml.safe_load(handle)

    try:
        config = GenomicRegionGeneratorConfig.model_validate(config_raw)
    except ValidationError as e:
        logging.error("Invalid configuration file:\n%s", e)
        raise

    # write used config
    write_config_to_yaml(config=config, dir_output=config.dir_output)

    pipeline = GenomicRegionGenerator(dir_output=config.dir_output)

    # generate the genomic regions
    # as the pydantic model is chosen depending on the `source` field,
    # the `source_params` contain the information from which source the
    # genomic data comes from and `source` is not needed anymore
    region_generator = pipeline.load_annotations(
        source_params=config.source,
    )

    files_fasta = pipeline.generate_genomic_regions(
        region_generator=region_generator,
        genomic_regions=config.genomic_regions,
        block_size=config.exon_exon_junction_block_size,
    )

    logging.info("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
