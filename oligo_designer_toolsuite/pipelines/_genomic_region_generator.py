############################################
# imports
############################################

import inspect
import logging
import os
from pathlib import Path

import yaml

from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.pipelines._utils import base_log_parameters, base_parser, setup_logging
from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
)

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
        source: str,
        mode: str | None,
        source_params: dict,
    ) -> CustomGenomicRegionGenerator:
        """
        Loads annotations from the specified source (NCBI, Ensembl, or custom files).

        :param source: The source of the annotations. Options: 'ncbi', 'ensembl', 'custom'.
        :type source: str
        :param mode: How should be specified which genome to use? Options: 'species', 'assembly'
        :type mode: str
        :param source_params: Parameters required for loading the annotations depending on the source.
            If source is 'ncbi', it should contain 'taxon', 'species', 'annotation_release' and can optionally
            contain 'assembly_source', 'refseq_assembly_accession', and 'assembly_name'.
            If source is 'ensembl', it should contain 'species' and 'annotation_release'.
            If source is 'custom', it should contain 'file_annotation', 'file_sequence', 'files_source',
            'species', 'annotation_release', and 'genome_assembly'.
        :type source_params: dict
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
        if source == "ncbi":
            if mode is None:
                raise ConfigurationError(
                    "For source='ncbi', top-level config parameter 'mode' must be provided."
                )
            # dowload the fasta files formthe NCBI server
            region_generator = NcbiGenomicRegionGenerator(
                mode=mode,
                taxon=source_params["taxon"],
                species=source_params["species"],
                annotation_release=source_params["annotation_release"],
                assembly_source=source_params["assembly_source"],
                refseq_assembly_accession=source_params["refseq_assembly_accession"],
                assembly_name=source_params["assembly_name"],
                dir_output=self.dir_output,
            )
        elif source == "ensembl":
            # dowload the fasta files formthe NCBI server
            region_generator = EnsemblGenomicRegionGenerator(
                species=source_params["species"],
                annotation_release=source_params["annotation_release"],
                dir_output=self.dir_output,
            )
        elif source == "custom":
            # use already dowloaded files
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
        genomic_regions: dict,
        block_size: int = 50,
    ) -> list:
        """
        Generates the specified genomic regions (e.g., genes, intergenic, exons, etc.) using the provided region generator.

        :param region_generator: An instance of CustomGenomicRegionGenerator that contains methods for generating
            genomic regions.
        :type region_generator: CustomGenomicRegionGenerator
        :param genomic_regions: A dictionary where keys are genomic region types (e.g., 'gene', 'intergenic', etc.)
            and values are flags indicating whether to generate that region.
        :type genomic_regions: dict
        :param block_size: Size of the block for genomic regions like exon-exon junctions. Defaults to 50.
        :type block_size: int
        :return: A list of file paths to the generated genomic regions.
        :rtype: list
        """
        files_fasta = []
        # loop not parallizeable due to file access restrictions
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
        config = yaml.safe_load(handle)

    pipeline = GenomicRegionGenerator(dir_output=config["dir_output"])

    # generate the genomic regions
    region_generator = pipeline.load_annotations(
        source=config["source"],
        mode=config["mode"],
        source_params=config["source_params"],
    )

    files_fasta = pipeline.generate_genomic_regions(
        region_generator=region_generator,
        genomic_regions=config["genomic_regions"],
        block_size=config["exon_exon_junction_block_size"],
    )

    logging.info("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
