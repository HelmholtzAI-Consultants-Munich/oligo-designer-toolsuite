############################################
# imports
############################################

import os
import yaml
import logging
import inspect
from pathlib import Path
from datetime import datetime

from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
)
from oligo_designer_toolsuite.pipelines._utils import log_parameters, base_parser

############################################
# Genomic Region Generator Functions
############################################


def load_annotations(
    source: str,
    source_params: dict,
    dir_output,
):
    """
    Loads genomic annotations based on the specified source which can be 'ncbi', 'ensembl', or 'custom'.
    This function configures and initiates the download or loading process for genomic annotations, depending on the source type specified.
    It also logs all operations and configurations used during the process.

    :param source: The source from which to load annotations. Supported options are 'ncbi', 'ensembl', or 'custom'.
    :type source: str
    :param source_params: Configuration parameters for the source. These parameters vary based on the selected source and may include taxon, species, annotation release, file paths, etc.
    :type source_params: dict
    :param dir_output: Directory where the output and logs will be stored.
    :type dir_output: str
    """
    ##### log parameters #####
    logging.info("Parameters Load Annotations:")
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    parameters = {i: values[i] for i in args}
    log_parameters(parameters)

    ##### loading annotations from different sources #####
    if source == "ncbi":
        # dowload the fasta files formthe NCBI server
        region_generator = NcbiGenomicRegionGenerator(
            taxon=source_params["taxon"],
            species=source_params["species"],
            annotation_release=source_params["annotation_release"],
            dir_output=dir_output,
        )
    elif source == "ensembl":
        # dowload the fasta files formthe NCBI server
        region_generator = EnsemblGenomicRegionGenerator(
            species=source_params["species"],
            annotation_release=source_params["annotation_release"],
            dir_output=dir_output,
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
            dir_output=dir_output,
        )
    else:
        raise ValueError(f"Source {source} not supported!")

    ##### save annotation information #####
    logging.info(
        f"The following annotation files are used for GTF annotation of regions: {region_generator.annotation_file} and for fasta sequence file: {region_generator.sequence_file} ."
    )
    logging.info(
        f"The annotations are from {region_generator.files_source} source, for the species: {region_generator.species}, release number: {region_generator.annotation_release} and genome assembly: {region_generator.genome_assembly}"
    )
    return region_generator


def generate_genomic_regions(
    region_generator: CustomGenomicRegionGenerator,
    genomic_regions: dict,
    block_size: int,
):
    """
    Generates specified genomic regions based on the provided flags and returns the paths to the generated fasta files.

    This function iterates through a dictionary of genomic regions and their associated flags, generating sequences for
    those flagged as true. It handles different types of regions such as genes, introns, and more, using the given block size
    for specific regions.

    :param region_generator: The generator used to create sequences for specified genomic regions.
    :type region_generator: CustomGenomicRegionGenerator
    :param genomic_regions: A dictionary specifying which genomic regions to generate, with boolean flags.
    :type genomic_regions: dict
    :param block_size: The block size to use for regions that require specific processing.
    :type block_size: int
    :return: A list containing paths to the fasta files generated for the specified regions.
    :rtype: list
    """
    fasta_files = []
    # loop not parallizeable due to file access restrictions
    for genomic_region, flag in genomic_regions.items():
        if flag:
            if genomic_region == "gene":
                fasta_file = region_generator.get_sequence_gene()
            elif genomic_region == "intergenic":
                fasta_file = region_generator.get_sequence_intergenic()
            elif genomic_region == "exon":
                fasta_file = region_generator.get_sequence_exon()
            elif genomic_region == "intron":
                fasta_file = region_generator.get_sequence_intron()
            elif genomic_region == "cds":
                fasta_file = region_generator.get_sequence_CDS()
            elif genomic_region == "utr":
                fasta_file = region_generator.get_sequence_UTR()
            elif genomic_region == "exon_exon_junction":
                fasta_file = region_generator.get_sequence_exon_exon_junction(block_size=block_size)
            else:
                raise Exception(f"Region generator: {genomic_region} is not implemented.")

            fasta_files.append(fasta_file)
            logging.info(f"The genomic region '{genomic_region}' was stored in :{fasta_file}.")

    return fasta_files


############################################
# Genomic Region Generator Pipeline
############################################


def main():
    """Main function to load sequence and annotation files from various sources and generate genomic regions of interest.
     Configurations are loaded from a YAML file specified at runtime. This function orchestrates the creation of genomic
     regions based on annotations from NCBI, Ensembl, or custom data sources. It handles the generation, logging,
     and storage of the regions. Each source requires specific parameters outlined in the function's detailed docstrings.

    If "ncbi" is choosen the source_params need to contain the keys:
    - "taxon": taxon of the species, valid taxa are: archaea, bacteria, fungi, invertebrate, mitochondrion, plant, plasmid, plastid, protozoa, vertebrate_mammalian, vertebrate_other, viral
    - "species": species name in NCBI download format, e.g. 'Homo_sapiens' for human; see [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/) for available species name
    - "annotation_release": release number (e.g. 109 or 109.20211119 for ncbi) of annotation or 'current' to use most recent annotation release. Check out release numbers for NCBI at ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/

    If "ensembl" is choosen the source_params need to contain the keys:
    - "species": species name in ensembl download format, e.g. 'homo_sapiens' for human; see http://ftp.ensembl.org/pub/release-108/gtf/ for available species names
    - "annotation_release": release number of annotation, e.g. 'release-108' or 'current' to use most recent annotation release. Check out release numbers for ensemble at ftp.ensembl.org/pub/

    If "custom" is choosen the source_params need to contain the keys:
    - "file_annotation": GTF file with gene annotation
    - "file_sequence": FASTA file with genome sequence
    - "files_source": original source of the genomic files -> optional, i.e. can be assigned None
    - "species": species of provided annotation, leave empty if unknown -> optional, i.e. can be assigned None
    - "annotation_release": release number of provided annotation, leave empty if unknown -> optional, i.e. can be assigned None
    - "genome_assembly": genome assembly of provided annotation, leave empty if unknown -> optional, i.e. can be assigned None
    """
    args = base_parser()

    # read the config file
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    # create the output folder
    dir_output = os.path.abspath(config["dir_output"])
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    ##### setup logger #####
    timestamp = datetime.now()
    file_logger = os.path.join(
        dir_output,
        f"log_genomic_region_generation_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
    )
    logging.getLogger("log_name")
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.NOTSET,
        handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
    )
    logging.captureWarnings(True)

    # generate the genomic regions
    region_generator = load_annotations(
        source=config["source"],
        source_params=config["source_params"],
        dir_output=dir_output,
    )

    fasta_files = generate_genomic_regions(
        region_generator=region_generator,
        genomic_regions=config["genomic_regions"],
        block_size=config["exon_exon_junction_block_size"],
    )


if __name__ == "__main__":
    main()
