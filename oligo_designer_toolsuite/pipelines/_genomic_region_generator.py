import logging
import inspect
import os
import yaml
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import datetime

from joblib import Parallel, delayed

from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
)
from oligo_designer_toolsuite.pipelines._utils import log_parameters


def load_annotations(
    source: str,
    source_params: dict,
    dir_output,
):
    """Load annotations for specified source. Source can be either "ncbi", "ensemble" or "custom".

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

    :param source: Indicate from where the annotation files will be loaded. Options:  'ncbi', 'ensembl', 'custom'.
    :type source: str
    :param source_params: Parameters for loading annotations. See above for details.
    :type source_params: dict
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

def generate_single_genomic_region(region_generator: CustomGenomicRegionGenerator, genomic_region: str, block_size: int):
    if genomic_region == "gene":
        return ["gene", region_generator.get_sequence_gene()]
    elif genomic_region == "intergenic":
        return ["intergenic", region_generator.get_sequence_intergenic()]
    elif genomic_region == "exon":
        return ["exon", region_generator.get_sequence_exon()]
    elif genomic_region == "intron":
        return ["intron", region_generator.get_sequence_intron()]
    elif genomic_region == "cds":
        return ["cds", region_generator.get_sequence_CDS()]
    elif genomic_region == "utr":
        return ["utr", region_generator.get_sequence_UTR()]
    elif genomic_region == "exon_exon_junction":
        return [
            "exon_exon_junction", 
            region_generator.get_sequence_exon_exon_junction(block_size=block_size)
        ]
    else:
        raise Exception(
            f"Region generator: {genomic_region} is not implemented."
        )

def generate_genomic_regions(region_generator: CustomGenomicRegionGenerator, genomic_regions: dict, block_size: int = None, n_jobs: int = 1):
    fasta_files = Parallel(n_jobs=n_jobs)(
        delayed(generate_single_genomic_region)(region_generator, genomic_region, block_size) for genomic_region, flag in genomic_regions.items() if flag
    )      
    return fasta_files


def main():
    parser = ArgumentParser(
        prog="Oligo Seq Designer",
        usage="oligo_seq [options]",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to the config yaml file, str",
        default=None,
        type=str,
        metavar="",
    )
    args = parser.parse_args()
    args = vars(args)
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

    ##### log parameters #####
    logging.info("Parameters Init:")
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    parameters = {i: values[i] for i in args}
    log_parameters(parameters)

    # generate the genomic regions
    region_generator = load_annotations(
        source=config["source"],
        source_params=config["source_params"],
        dir_output=dir_output,
    )
    import time
    start = time.time()
    genomic_regions = generate_genomic_regions(region_generator=region_generator, genomic_regions=config["genomic_regions"], block_size=config["exon_exon_junction_block_size"], n_jobs=config["n_jobs"])
    for genomic_region in genomic_regions:
        logging.info(f"The genomic region {genomic_region[0]} was stored in :{genomic_region[1]}.")
    logging.info(f"time required: {time.time() -start}")
if __name__ == "__main__":
    main()
