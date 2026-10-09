#!/usr/bin/env bash
# Copy a file from the ODT repository (config template, example table from data/functional_oligos/,
# or environment.yaml) to a destination.
# Usage: get_odt_file.sh <repo_path> <destination>
#   e.g. get_odt_file.sh data/configs/scrinshot_probe_designer.yaml ./scrinshot_probe_designer.yaml
#        get_odt_file.sh environment.yaml ./environment.yaml
#
# Source, in order of preference:
#   1. the ODT checkout this skill lives in (skills/odt-probe-design -> ../../), so the templates in
#      data/configs are the single source of truth,
#   2. GitHub at tag ODT_REF (default v2.3.0), when the skill was installed on its own (plugin, zip).
set -euo pipefail

ODT_REF="${ODT_REF:-v2.3.0}"
RAW_URL="https://raw.githubusercontent.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/${ODT_REF}"

if [ $# -ne 2 ]; then
    sed -n '2,11p' "$0" >&2
    exit 1
fi
repo_path="$1"
destination="$2"

skill_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
repo_root="$(cd "$skill_dir/../.." && pwd)"

if [ -f "$repo_root/$repo_path" ] && [ -f "$repo_root/pyproject.toml" ] \
    && grep -q 'name = "oligo-designer-toolsuite"' "$repo_root/pyproject.toml"; then
    cp "$repo_root/$repo_path" "$destination"
    echo "Copied $repo_path from local ODT checkout ($repo_root) to $destination"
elif command -v curl >/dev/null 2>&1 && curl -fsSL "$RAW_URL/$repo_path" -o "$destination"; then
    echo "Downloaded $repo_path from GitHub ($ODT_REF) to $destination"
elif command -v wget >/dev/null 2>&1 && wget -q "$RAW_URL/$repo_path" -O "$destination"; then
    echo "Downloaded $repo_path from GitHub ($ODT_REF) to $destination"
else
    rm -f "$destination"
    echo "Could not get $repo_path: no local ODT checkout and no access to $RAW_URL/$repo_path" >&2
    exit 1
fi
