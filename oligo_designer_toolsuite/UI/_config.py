import os
from pathlib import Path
from distutils.dir_util import copy_tree

import yaml
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap


def generate_config(exp_name: str, directory: str = "output", source: str = "custom"):
    """
    Generate a configuration based on the provided experimental name and source.

    :param exp_name: Name of the experiment for which the configuration is to be generated, can be "MERFISH", "SCRINSHOT" and "SeqFISHPlus" (case insensitive).
    :type  exp_name: str
    :param directory: Directory where the user's configuration will be saved. Defaults to "output".
    :type  directory: str
    :param source: Source of the data, can be "ensembl", "ncbi", or "custom". Defaults to "custom".
    :type  source: str
    :return: Configuration dictionary merged from basic and advanced configurations.
    :rtype: dict
    """
    Path(directory).mkdir(parents=True, exist_ok=True)

    config_parent_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "data",
        "configs",
        exp_name.lower(),
    )

    config_user_dir = os.path.join(directory, "config")
    copy_tree(config_parent_dir, config_user_dir)

    adapt_config_to_source(config_user_dir, source)
    config = load_files_into_dict(config_parent_dir)
    return config


def adapt_config_to_source(config_user_dir, source):
    """
    Adapt the configuration based on the provided source.

    :param config_user_dir: Directory where the user's configuration is saved.
    :type  config_user_dir: str
    :param source: Source of the data. Can be "ensembl", "ncbi", or "custom".
    :type  source: str
    :raises ValueError: If the provided source is not supported.
    """
    file_path = os.path.join(config_user_dir, "basic.yaml")

    yaml = YAML()
    yaml.preserve_quotes = True

    sources_config = {
        "ensembl": {
            "params": {"species": "homo_sapiens", "annotation_release": "current"},
            "comments": {
                "species": "required: species name in ensemble download format, e.g. 'homo_sapiens' for human; see http://ftp.ensembl.org/pub/release-108/gtf/ for available species names",
                "annotation_release": "required: release number of annotation, e.g. 'release-108' or 'current' to use most recent annotation release. Check out release numbers for ensemble at ftp.ensembl.org/pub/",
            },
        },
        "ncbi": {
            "params": {
                "taxon": "vertebrate_mammalian",
                "species": "Homo_sapiens",
                "annotation_release": "110",
            },
            "comments": {
                "taxon": "required: taxon of the species, valid taxa are: archaea, bacteria, fungi, invertebrate, mitochondrion, plant, plasmid, plastid, protozoa, vertebrate_mammalian, vertebrate_other, viral",
                "species": "required: species name in NCBI download format, e.g. 'Homo_sapiens' for human; see https://ftp.ncbi.nlm.nih.gov/genomes/refseq/ for available species name",
                "annotation_release": "required: release number of annotation e.g. '109' or '109.20211119'  or 'current' to use most recent annotation release. Check out release numbers for NCBI at ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/",
            },
        },
        "custom": {
            "params": {
                "file_annotation": "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf",
                "file_sequence": "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna",
                "files_source": "NCBI",
                "species": "Homo_sapiens",
                "annotation_release": "110",
                "genome_assembly": "GRCh38",
            },
            "comments": {
                "file_annotation": "required: GTF file with gene annotation",
                "file_sequence": "required: FASTA file with genome sequence",
                "files_source": "optional: original source of the genomic files",
                "species": "optional: species of provided annotation, leave empty if unknown",
                "annotation_release": "optional: release number of provided annotation, leave empty if unknown",
                "genome_assembly": "optional: genome assembly of provided annotation, leave empty if unknown",
            },
        },
    }

    if source not in sources_config:
        raise ValueError(f"Unsupported source: {source}")

    params = sources_config[source]["params"]
    comments = sources_config[source]["comments"]

    with open(file_path, "r") as f:
        config = yaml.load(f)

    # Update the configuration
    config["source"] = source
    config["source_params"] = CommentedMap(params.items())
    config["file_genes"] = config["file_genes"].format(source=source)

    # Insert comments
    for key, comment in comments.items():
        config["source_params"].yaml_add_eol_comment(comment, key)

    # Save the updated data back to the YAML file
    with open(file_path, "w") as f:
        yaml.dump(config, f)


def merge_dicts(dict1, dict2):
    """Recursive function to merge two dictionaries."""
    merged = dict1.copy()  # Start with dict1's keys and values
    for key, value in dict2.items():
        if key in dict1 and isinstance(dict1[key], dict) and isinstance(value, dict):
            merged[key] = merge_dicts(dict1[key], value)  # Recurse into subdictionaries
        else:
            merged[key] = value  # Overwrite or add new key-value pairs
    return merged


def load_files_into_dict(config_parent_dir):
    """
    Helper function to load configuration from basic and advanced YAML files into a dictionary.

    :param config_parent_dir: Directory from where the configuration files are to be read.
    :type  config_parent_dir: str
    :return: Merged configuration dictionary from basic and advanced configurations.
    :rtype: dict
    """
    basic_file_path = os.path.join(config_parent_dir, "basic.yaml")
    advanced_file_path = os.path.join(config_parent_dir, "advanced.yaml")

    with open(basic_file_path, "r") as f:
        basic_config = yaml.safe_load(f)
    with open(advanced_file_path, "r") as f:
        advanced_config = yaml.safe_load(f)

    merged_config = merge_dicts(basic_config, advanced_config)

    return merged_config
