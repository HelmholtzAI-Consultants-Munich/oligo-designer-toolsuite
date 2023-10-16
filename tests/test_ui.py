############################################
# imports
############################################

import os
import yaml
from argparse import ArgumentParser

from oligo_designer_toolsuite.UI import (
    initialize_parameters,
    extract_arguments_to_dict,
)

from oligo_designer_toolsuite.UI._config import load_files_into_dict

############################################
# Global Parameters
############################################

config_dir = "tests/data/UI"


############################################
# Tests
############################################


def test_initialize_parameters():
    parser = ArgumentParser()
    initialize_parameters(parser, config_dir)

    # Check if arguments were added to the parser correctly
    args = parser.parse_args(
        [
            "--n-jobs",
            "1",
            "--dir-output",
            "output",
            "--min-probes-per-gene",
            "0",
            "--write-removed-genes",
            "True",
            "--write-intermediate-steps",
            "True",
            "--source",
            "custom",
            "--source-params--files-source",
            "NCBI",
            "--source-params--species",
            "Homo_sapiens",
            "--source-params--annotation-release",
            "110",
            "--source-params--genome-assembly",
            "GRCh38",
            "--source-params--file-annotation",
            "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf",
            "--source-params--file-sequence",
            "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna",
            "--region",
            "transcript",
            "--Tm-chem-correction-param-detection-oligo--DMSO",
            "0",
            "--Tm-chem-correction-param-detection-oligo--DMSOfactor",
            "0.75",
            "--Tm-chem-correction-param-detection-oligo--fmdfactor",
            "0.65",
            "--Tm-chem-correction-param-detection-oligo--fmd",
            "30",
        ]
    )
    assert args.n_jobs == 1, "n_jobs argument initialization failed"
    assert (
        args.write_removed_genes == True
    ), "write_removed_genes argument initialization failed"
    assert (
        args.write_intermediate_steps == True
    ), "write_intermediate_steps argument initialization failed"
    assert (
        args.Tm_chem_correction_param_detection_oligo__DMSO == 0
    ), "Tm_chem_correction_param_detection_oligo_DMSO argument initialization failed"
    assert (
        args.Tm_chem_correction_param_detection_oligo__DMSOfactor == 0.75
    ), "Tm_chem_correction_param_detection_oligo_DMSOfactor argument initialization failed"
    assert (
        args.Tm_chem_correction_param_detection_oligo__fmdfactor == 0.65
    ), "Tm_chem_correction_param_detection_oligo_fmdfactor argument initialization failed"
    assert (
        args.Tm_chem_correction_param_detection_oligo__fmd == 30
    ), "Tm_chem_correction_param_detection_oligo_fmd argument initialization failed"
    assert args.dir_output == "output", "dir_output argument initialization failed"
    assert (
        args.min_probes_per_gene == 0
    ), "min_probes_per_gene argument initialization failed"
    assert args.source == "custom", "source argument initialization failed"
    assert (
        args.source_params__files_source == "NCBI"
    ), "source_params_files_source argument initialization failed"
    assert (
        args.source_params__species == "Homo_sapiens"
    ), "source_params_species argument initialization failed"
    assert (
        args.source_params__annotation_release == 110
    ), "source_params_annotation_release argument initialization failed"
    assert (
        args.source_params__genome_assembly == "GRCh38"
    ), "source_params_genome_assembly argument initialization failed"
    assert (
        args.source_params__file_annotation
        == "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
    ), "source_params_file_annotation argument initialization failed"
    assert (
        args.source_params__file_sequence
        == "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
    ), "source_params_file_sequence argument initialization failed"


def test_extract_arguments_to_dict():
    expected_config = load_files_into_dict(config_dir)

    parser = ArgumentParser()
    initialize_parameters(parser, config_dir)
    args = parser.parse_args(
        [
            "--n-jobs",
            "1",
            "--dir-output",
            "output",
            "--min-probes-per-gene",
            "0",
            "--write-removed-genes",
            "True",
            "--write-intermediate-steps",
            "True",
            "--source",
            "custom",
            "--source-params--files-source",
            "NCBI",
            "--source-params--species",
            "Homo_sapiens",
            "--source-params--annotation-release",
            "110",
            "--source-params--genome-assembly",
            "GRCh38",
            "--source-params--file-annotation",
            "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf",
            "--source-params--file-sequence",
            "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna",
            "--region",
            "transcript",
            "--Tm-chem-correction-param-detection-oligo--DMSO",
            "0",
            "--Tm-chem-correction-param-detection-oligo--DMSOfactor",
            "0.75",
            "--Tm-chem-correction-param-detection-oligo--fmdfactor",
            "0.65",
            "--Tm-chem-correction-param-detection-oligo--fmd",
            "30",
        ]
    )

    extracted_config = extract_arguments_to_dict(parser, args)

    assert (
        extracted_config == expected_config
    ), f"Failed to convert the parsed arguments to dict. \n\nExpected:\n{expected_config}\n\nGot:\n{extracted_config}"


test_initialize_parameters()
test_extract_arguments_to_dict()
