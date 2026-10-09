#!/usr/bin/env python3
"""Replace @GRG:<region>@ placeholders in an ODT config with FASTA files from the genomic region generator.

Usage: resolve_genomic_region_paths.py <config_in.yaml> <config_out.yaml> <grg_dir_output>

Placeholders are used for NCBI/Ensembl downloads because the generated file names contain the
annotation release and genome assembly, which ODT only resolves on the FTP server at download time.
The genomic region generator writes files to <grg_dir_output>/annotation/ named
"...__annotation_type__<region>[__...].fna". Each placeholder must match exactly one file.
Only uses the standard library so it runs on the host or inside the ODT container.
"""

import glob
import os
import re
import sys

PATTERNS = {
    "gene": "*__annotation_type__gene.fna",
    "intergenic": "*__annotation_type__intergenic.fna",
    "exon": "*__annotation_type__exon.fna",
    "intron": "*__annotation_type__intron.fna",
    "cds": "*__annotation_type__cds.fna",
    "utr": "*__annotation_type__utr*.fna",
    "exon_exon_junction": "*__annotation_type__exon_exon_junction__block_size__*.fna",
}


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    file_in, file_out, dir_grg = sys.argv[1:]
    dir_annotation = os.path.join(dir_grg, "annotation")

    def resolve(match: re.Match) -> str:
        region = match.group(1)
        if region not in PATTERNS:
            sys.exit(f"Unknown genomic region placeholder '@GRG:{region}@'. Valid: {', '.join(PATTERNS)}")
        files = sorted(glob.glob(os.path.join(dir_annotation, PATTERNS[region])))
        if len(files) != 1:
            sys.exit(
                f"Expected exactly one '{region}' FASTA in {dir_annotation}, found {len(files)}: {files}. "
                f"Was '{region}' enabled in the genomic region generator config, and did the generator finish?"
            )
        return os.path.relpath(files[0])

    with open(file_in) as handle:
        text = handle.read()
    text = re.sub(r"@GRG:([a-z_]+)@", resolve, text)
    with open(file_out, "w") as handle:
        handle.write(text)
    print(f"Wrote {file_out}")


if __name__ == "__main__":
    main()
