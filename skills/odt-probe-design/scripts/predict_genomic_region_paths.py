#!/usr/bin/env python3
"""Print the FASTA paths the ODT genomic region generator will write for a `source: custom` config.

Usage: predict_genomic_region_paths.py genomic_region_generator.yaml

Output: one line per enabled region, `<region>\t<path>`, with paths relative to the working
directory the generator is run from. Needs PyYAML (available in the ODT runtime).

Only for `source: custom`: there the file names are fully determined by the config. For NCBI and
Ensembl, ODT resolves the annotation release and genome assembly on the FTP server at download time,
so the probe config uses @GRG:<region>@ placeholders instead (see resolve_genomic_region_paths.py).
Naming follows ODT v2.3.0.
"""

import os
import sys

import yaml

REGION_SUFFIX = {
    "gene": "gene",
    "intergenic": "intergenic",
    "exon": "exon",
    "intron": "intron",
    "cds": "cds",
    "utr": "utr__five_prime__three_prime",
}


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    with open(sys.argv[1]) as handle:
        config = yaml.safe_load(handle)

    if str(config["source"]).lower() != "custom":
        sys.exit(
            f"source '{config['source']}': file names are only known after the download, "
            "use @GRG:<region>@ placeholders in the probe config instead."
        )

    params = config["source_params"]
    file_info = (
        f"source__{params.get('files_source') or 'custom'}"
        f"__species__{params.get('species') or 'unknown'}"
        f"__annotation_release__{params.get('annotation_release') or 'unknown'}"
        f"__genome_assemly__{params.get('genome_assembly') or 'unknown'}"
    )
    dir_annotation = os.path.join(config["dir_output"], "annotation")
    for region, enabled in config["genomic_regions"].items():
        if not enabled:
            continue
        if region == "exon_exon_junction":
            suffix = f"exon_exon_junction__block_size__{config['exon_exon_junction_block_size']}"
        else:
            suffix = REGION_SUFFIX[region]
        print(f"{region}\t{os.path.join(dir_annotation, f'{file_info}__annotation_type__{suffix}.fna')}")


if __name__ == "__main__":
    main()
