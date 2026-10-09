#!/usr/bin/env bash
# Report which runtimes are available for the Oligo Designer Toolsuite (ODT) and whether ODT is installed.
# Usage: check_odt.sh [conda_env_name] [path_to_sif]
# Read-only: downloads and installs nothing.

ODT_VERSION="2.3.0"
IMAGE="ghcr.io/helmholtzai-consultants-munich/odt:${ODT_VERSION}"
ENV_NAME="${1:-odt}"
SIF="${2:-}"

echo "== System"
echo "os: $(uname -s) $(uname -m)"
echo "cpus: $(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null)"

echo "== Docker"
if command -v docker >/dev/null 2>&1; then
    echo "docker: $(docker --version)"
    if docker info >/dev/null 2>&1; then
        echo "daemon: running"
        if docker image inspect "$IMAGE" >/dev/null 2>&1; then
            echo "odt_image: present ($IMAGE)"
        else
            echo "odt_image: missing ($IMAGE)"
        fi
    else
        echo "daemon: not running or no permission"
    fi
else
    echo "docker: not found"
fi

echo "== Apptainer/Singularity"
APPT=""
for c in apptainer singularity; do
    if command -v "$c" >/dev/null 2>&1; then APPT="$c"; echo "$c: $($c --version)"; break; fi
done
[ -z "$APPT" ] && echo "apptainer: not found"
if [ -z "$SIF" ]; then
    SIF=$(ls ./odt_${ODT_VERSION}.sif "$HOME"/odt_${ODT_VERSION}.sif 2>/dev/null | head -1)
fi
if [ -n "$SIF" ] && [ -f "$SIF" ]; then echo "odt_sif: present ($SIF)"; else echo "odt_sif: not found (looked for odt_${ODT_VERSION}.sif in . and \$HOME)"; fi

echo "== conda/mamba"
CONDA=""
for c in mamba micromamba conda; do
    if command -v "$c" >/dev/null 2>&1; then CONDA="$c"; echo "$c: $($c --version 2>&1 | head -1)"; break; fi
done
if [ -n "$CONDA" ]; then
    if $CONDA env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "env '$ENV_NAME': present"
        v=$($CONDA run -n "$ENV_NAME" pip show oligo-designer-toolsuite 2>/dev/null | awk '/^Version:/{print $2}')
        echo "odt_in_env: ${v:-not installed}"
        for t in blastn bedtools bowtie bowtie2 samtools bcftools; do
            $CONDA run -n "$ENV_NAME" which "$t" >/dev/null 2>&1 && echo "  $t: ok" || echo "  $t: MISSING"
        done
    else
        echo "env '$ENV_NAME': not found"
    fi
else
    echo "conda: not found"
fi

echo "== ODT on PATH (current shell)"
if command -v genomic_region_generator >/dev/null 2>&1; then
    v=$(pip show oligo-designer-toolsuite 2>/dev/null | awk '/^Version:/{print $2}')
    echo "odt_cli: $(command -v genomic_region_generator) (version ${v:-unknown})"
else
    echo "odt_cli: not on PATH"
fi

echo "== Wanted version: ${ODT_VERSION}"
