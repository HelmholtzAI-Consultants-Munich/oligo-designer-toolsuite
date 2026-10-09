#!/usr/bin/env bash
# Run script generated for the Oligo Designer Toolsuite (ODT) v2.3.0.
# Pipeline: __PIPELINE_COMMAND__
# Usage:    bash run_odt.sh
# Long runs: nohup bash run_odt.sh > run_odt.log 2>&1 &   (or wrap in an sbatch script on a cluster)
#
# NOTE on @GRG:<region>@ placeholders in the probe config: when the genomic region generator downloads
# from NCBI or Ensembl, the names of the FASTA files it writes contain the annotation release and genome
# assembly, which are only known after the download (e.g. release "current" -> "116"). The probe config
# therefore uses placeholders such as @GRG:exon@. After step 1 this script replaces them with the real
# file names and writes <probe config>.resolved.yaml, which is the config that is actually run.
# The FASTA files only exist after the generator (step 1) has finished successfully.
set -euo pipefail

########## settings (filled in by Claude) ##########
RUNTIME="__RUNTIME__"            # docker | apptainer | conda | native
IMAGE="ghcr.io/helmholtzai-consultants-munich/odt:2.3.0"   # docker
SIF="__SIF_PATH__"               # apptainer, e.g. ./odt_2.3.0.sif
CONDA_ENV="__CONDA_ENV__"        # conda, e.g. odt
PROJECT_DIR="__PROJECT_DIR__"    # absolute path; mounted into the container at the same path
EXTRA_MOUNTS=(__EXTRA_MOUNTS__)  # directories outside PROJECT_DIR that inputs live in, e.g. ("/data/genomes")
# Docker/Apptainer: the container only sees PROJECT_DIR and EXTRA_MOUNTS (Apptainer also $HOME, /tmp and
# the current directory by default). If any file in the configs lies elsewhere, including the target of
# a symlink inside PROJECT_DIR, add its directory to EXTRA_MOUNTS, otherwise ODT reports it as not found.

RUN_GENOMIC_REGION_GENERATOR=__RUN_GRG__   # true | false
CONFIG_GRG="genomic_region_generator.yaml"
GRG_DIR_OUTPUT="__GRG_DIR_OUTPUT__"        # must equal dir_output in CONFIG_GRG

PIPELINE_COMMAND="__PIPELINE_COMMAND__"    # e.g. scrinshot_probe_designer
CONFIG_PROBES="__PROBE_CONFIG__"           # e.g. scrinshot_probe_designer.yaml (may contain @GRG:...@ placeholders)
####################################################

cd "$PROJECT_DIR"

odt() {
    case "$RUNTIME" in
        docker)
            local mounts=(-v "$PROJECT_DIR":"$PROJECT_DIR")
            for m in ${EXTRA_MOUNTS[@]+"${EXTRA_MOUNTS[@]}"}; do mounts+=(-v "$m":"$m"); done
            docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp "${mounts[@]}" -w "$PWD" "$IMAGE" "$@"
            ;;
        apptainer)
            local cmd
            cmd=$(command -v apptainer || command -v singularity)
            local binds="$PROJECT_DIR"
            for m in ${EXTRA_MOUNTS[@]+"${EXTRA_MOUNTS[@]}"}; do binds="$binds,$m"; done
            "$cmd" exec --bind "$binds" "$SIF" "$@"
            ;;
        conda)
            conda run -n "$CONDA_ENV" --no-capture-output "$@"
            ;;
        native)
            "$@"
            ;;
        *)
            echo "Unknown RUNTIME '$RUNTIME'" >&2; exit 1
            ;;
    esac
}

echo "[$(date '+%F %T')] Using runtime: $RUNTIME"

if [ "$RUN_GENOMIC_REGION_GENERATOR" = true ]; then
    echo "[$(date '+%F %T')] Step 1: genomic region generator"
    odt genomic_region_generator -c "$CONFIG_GRG"
fi

CONFIG_TO_RUN="$CONFIG_PROBES"
if grep -q "@GRG:" "$CONFIG_PROBES"; then
    CONFIG_TO_RUN="${CONFIG_PROBES%.yaml}.resolved.yaml"
    echo "[$(date '+%F %T')] Replacing @GRG:...@ placeholders with the generated FASTA files -> $CONFIG_TO_RUN"
    odt python3 resolve_genomic_region_paths.py "$CONFIG_PROBES" "$CONFIG_TO_RUN" "$GRG_DIR_OUTPUT"
fi

# The probe designer needs all FASTA files listed in its config; check they exist before the long run.
missing=()
while IFS= read -r file_fasta; do
    [ -f "$file_fasta" ] || missing+=("$file_fasta")
done < <(grep -oE "[^[:space:]\"']+\.(fna|fa|fasta)" "$CONFIG_TO_RUN" | sort -u)
if [ ${#missing[@]} -gt 0 ]; then
    echo "ERROR: FASTA files listed in $CONFIG_TO_RUN do not exist:" >&2
    printf '  %s\n' "${missing[@]}" >&2
    if [ "$RUN_GENOMIC_REGION_GENERATOR" = true ]; then
        echo "The genomic region generator has to finish successfully first. Files it created:" >&2
        ls "$GRG_DIR_OUTPUT/annotation/"*.fna >&2 2>/dev/null || echo "  (none in $GRG_DIR_OUTPUT/annotation/)" >&2
    fi
    exit 1
fi

echo "[$(date '+%F %T')] Step 2: $PIPELINE_COMMAND"
odt "$PIPELINE_COMMAND" -c "$CONFIG_TO_RUN"

echo "[$(date '+%F %T')] Done. Results are in the dir_output set in $CONFIG_TO_RUN"
