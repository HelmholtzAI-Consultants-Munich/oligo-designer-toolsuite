import os
from pathlib import Path
import warnings
import yaml
from argparse import ArgumentParser


def generate_config_file(
    exp_name: str, directory: str = "output", source: str = "custom"
):
    directory = os.path.join(directory, "config")
    Path(directory).mkdir(parents=True, exist_ok=True)

    config_parent_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "data",
        "configs",
    )

    if source == "custom":
        # function generating the config file
        # config_file = generate_custom_config(directory)
        config_file = os.path.join(
            config_parent_dir, f"{exp_name}_probe_designer_custom.yaml"
        )
    elif source == "ncbi":
        # function generating the config file
        # config_file = generate_ncbi_config(directory)
        config_file = os.path.join(
            config_parent_dir, f"{exp_name}_probe_designer_ncbi.yaml"
        )
    elif source == "ensembl":
        # function generating the config file
        # config_file = generate_ensembl_config(directory)
        config_file = os.path.join(
            config_parent_dir, f"{exp_name}_probe_designer_ensembl.yaml"
        )
    else:
        config_file = ""
        raise ValueError(f"No config file found for source {source}'.")

    warnings.warn(f"Using default config: {config_file}.")
    return config_file


def initialize_parameters(parser: ArgumentParser, exp_name: str):
    # Common arguments
    parser.add_argument(
        "-o",
        "--output",
        help="path to the output folder, str",
        required=True,
        type=str,
        metavar="",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to the config yaml file, str",
        default=None,
        type=str,
        metavar="",
    )
    parser.add_argument(
        "-n",
        "--n_jobs",
        help="number of cores used, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-mog",
        "--min_probes_per_gene",
        help="genes with less that this number of probes are removed, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-wg",
        "--write_removed_genes",
        help="write in a file the removed genes, bool",
        type=bool,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-ws",
        "--write_intermediate_steps",
        help="write the probe sequences after each step of the pipeline, bool",
        type=bool,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-s",
        "--source",
        help="how to obtain the genomic files: download them from a server [ncbi, ensembl] or provide the files directly [custom]",
        choices=["ncbi", "ensembl", "custom"],
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-tx",
        "--taxon",
        help="[source: ncbi] taxon of the species, [archaea, bacteria, fungi, invertebrate, mitochondrion, plant, plasmid, plastid, protozoa, vertebrate_mammalian, vertebrate_other, viral]",
        choices=[
            "archaea",
            "bacteria",
            "fungi",
            "invertebrate",
            "mitochondrion",
            "plant",
            "plasmid",
            "plastid",
            "protozoa",
            "vertebrate_mammalian",
            "vertebrate_other",
            "viral",
        ],
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-sp",
        "--species",
        help="species name, for valid NCBI species name see https://ftp.ncbi.nlm.nih.gov/genomes/refseq/, for valid Ensembl species name see http://ftp.ensembl.org/pub/release-108/gtf/, str",
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--annotation_release",
        help="release number of annotation, e.g. 'release-108' (Ensembl) or '109' (NCBI) or 'current' to use most recent annotation release, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--genome_assembly",
        help="[source: custom] genome assembly of provided annotation, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-fa",
        "--file_annotation",
        help="[source: custom] path to GTF file with gene annotation, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-fs",
        "--file_sequence",
        help="[source: custom] path to FASTA file with genome sequence, str",
        type=str,
        default=None,
        metavar="",
    )

    parser.add_argument(
        "-fsrc",
        "--files_source",
        help="[source: custom] original source of the genomic files, e.g. NCBI, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-fg",
        "--file_genes",
        help="path to file with a list of the genes that are used to generate the probes sequences, if empty all the genes are used, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-olm",
        "--probe_length_min",
        help="minimum length of probes, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-olM",
        "--probe_length_max",
        help="max length of probes, int",
        type=int,
        default=None,
        metavar="",
    )

    # TODO: adapt here + primers + readouts
    parser.add_argument(
        "--target_GC_content_min",
        help="minimum GC content of target probes, [0, 100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--target_GC_content_max",
        help="maximum GC content of target probes, [0, 100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--target_Tm_min",
        help="minimum melting temperature of target probes, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--target_Tm_max",
        help="maximum melting temperature of target probes, float",
        type=float,
        default=None,
        metavar="",
    )

    # TODO: differentiate all word sizes
    parser.add_argument(
        "--blast_word_size",
        help="word size for the blastn seed, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--blast_coverage",
        help="minimum coverage between probes and target sequence for blastn, [0,100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--blast_percent_identity",
        help="maximum similarity between probes and target sequences for blastn, [0,100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--strand",
        help="strand of the query sequence to search, str",
        type=str,
        default=None,
        metavar="",
    )

    ### Specific parts:
    if exp_name == "scrinshot":
        parser.add_argument(
            "--min_arm_length",
            help="minimum length of each arm, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--max_arm_Tm_dif",
            help="max melting temperature difference of both arms, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--arm_Tm_min",
            help="minimum melting temperature of each arm, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--arm_Tm_max",
            help="max melting temperature of each arm, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--ligation_region_size",
            help="size of the seed region around the ligation site for bowtie seed region filter, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--Tm_opt",
            help="optimal melting temperature of probes, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--GC_content_opt",
            help="optimal GC content of probes, [0, 100]",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--Tm_weight",
            help="weight of the Tm of the probe in the efficiency score, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--GC_weight",
            help="weight of the GC content of the probe in the efficiency score, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--probeset_size_opt",
            help="ideal number of probes per probeset, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--probeset_size_min",
            help="minimum number of probes per probeset, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--n_sets",
            help="maximum number of sets per gene, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--detect_oligo_length_min",
            help="minimum number of probes per probeset, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--detect_oligo_length_max",
            help="maximum length of detection probe, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--detect_oligo_Tm_opt",
            help="optimal melting temperature of detection probe, float",
            type=float,
            default=None,
            metavar="",
        )

    if exp_name == "merfish":
        # target
        parser.add_argument(
            "--target_internal_secondary_structures_T",
            help="temperature at which secondary structures are evaluated, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--target_internal_secondary_structures_threshold_deltaG",
            help="delta G threshold of the secondary structures, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--target_max_repeats_nt",
            help="maximum number of consecutive repeats in target probes, int",
            type=int,
            default=None,
            metavar="",
        )

        # primers
        parser.add_argument(
            "--primer_GC_content_min",
            help="minimum GC content of primer probes, [0, 100]",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--primer_GC_content_max",
            help="maximum GC content of primer probes, [0, 100]",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--primer_Tm_min",
            help="minimum melting temperature of primer probes, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--primer_Tm_max",
            help="maximum melting temperature of primer probes, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--primer_max_repeats_nt",
            help="maximum number of consecutive repeats in primer probes, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--primer_GC_clamp_n",
            help="number of terminal bases to check for a G or a C, int",
            type=int,
            default=None,
            metavar="",
        )

        # readouts
        parser.add_argument(
            "--readout_GC_content_min",
            help="minimum GC content of readout probes, [0, 100]",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--readout_GC_content_max",
            help="maximum GC content of readout probes, [0, 100]",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--readout_Tm_min",
            help="minimum melting temperature of readout probes, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--readout_Tm_max",
            help="maximum melting temperature of readout probes, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--readout_max_repeats_nt",
            help="maximum number of consecutive repeats in readout probes, int",
            type=int,
            default=None,
            metavar="",
        )
        # blast
        # target
        parser.add_argument(
            "--target_percent_identity_ch",
            help="percent_identity for cross hybridization, float",
            type=float,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--target_blast_word_size",
            help="word size for the blastn seed for target probes, int",
            type=int,
            default=None,
            metavar="",
        )

        # TODO: better naming
        # probe
        parser.add_argument(
            "--probe_blast1_word_size",
            help="word size for the blastn seed for probes, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--probe_blast2_word_size",
            help="word size for the blastn seed for probes, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--file_highly_expressed_genes",
            help="text file containing highly expressed genes, str",
            type=str,
            default=None,
            metavar="",
        )

        # readout
        parser.add_argument(
            "--readout_blast1_word_size",
            help="word size for the blastn seed for readout probes, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--readout_blast2_word_size",
            help="word size for the blastn seed for readout probes, int",
            type=int,
            default=None,
            metavar="",
        )
        # primers
        parser.add_argument(
            "--primer_blast1_word_size",
            help="word size for the blastn seed for primer probes, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--primer_blast2_word_size",
            help="word size for the blastn seed for primer probes, int",
            type=int,
            default=None,
            metavar="",
        )
        parser.add_argument(
            "--primer_blast3_word_size",
            help="word size for the blastn seed for primer probes, int",
            type=int,
            default=None,
            metavar="",
        )

        # TODO: better description
        parser.add_argument(
            "--encoding",
            help="Encoding, str",
            type=str,
            default="MHD4",
            metavar="",
        )

    args = parser.parse_args()
    args = vars(args)
    if args["config"] is None:
        warnings.warn(f"No config file defined. Creating default config.")
        if args["source"] is None:
            warnings.warn(f"No source was defined. Using default source: custom")
            args["source"] = "custom"
        args["config"] = generate_config_file(exp_name, args["output"], args["source"])

    # read the config file
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    # update the config file values with the given one in the command line
    for param, value in args.items():
        if value is not None and param != "config":
            # overwrite the config file
            if param in [
                "file_annotation",
                "file_sequence",
                "files_source",
                "taxon",
                "species",
                "annotation_release",
                "genome_assembly",
            ]:
                # source specific parameters
                config["source_params"][param] = value
            elif param in [
                "target_GC_content_min",
                "target_GC_content_max",
                "target_Tm_min",
                "target_Tm_max",
                "target_internal_secondary_structures_T",
                "target_internal_secondary_structures_threshold_deltaG",
                "target_max_repeats_nt",
            ]:
                config["targets_setup"][param.split("target_")[1]] = value
            elif param in [
                "readout_GC_content_min",
                "readout_GC_content_max",
                "readout_Tm_min",
                "readout_Tm_max",
                "readout_max_repeats_nt",
            ]:
                config["readout_setup"][param.split("readout_")[1]] = value
            elif param in [
                "primer_GC_content_min",
                "primer_GC_content_max",
                "primer_GC_clamp_n",
                "primer_Tm_min",
                "primer_Tm_max",
                "primer_max_repeats_nt",
            ]:
                config["primers_setup"][param.split("primers_")[1]] = value
            elif param == "encoding":
                if value == "MHD4":
                    config["MHD4"] = {
                        "num_seq": 100,
                        "num_bits": 16,
                        "num_ones": 4,
                        "hamming_distance": 4,
                    }
                elif value == "MHD2":
                    config["MHD2"] = {
                        "num_seq": 100,
                        "num_bits": 12,
                        "num_ones": 4,
                        "hamming_distance": 2,
                    }
            else:
                config[param] = value

    return config
