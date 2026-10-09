# Installing ODT v2.3.0

Supported: Linux x64 and macOS arm64 (Intel macs are deprecated, Windows untested — use Docker/WSL2 there).
Python 3.10–3.12.

All pipelines are exposed as CLI commands taking one argument `-c/--config <yaml>`:

```
genomic_region_generator   scrinshot_probe_designer   merfish_probe_designer
seqfish_plus_probe_designer   hcr_probe_designer   cycle_hcr_probe_designer   oligo_seq_probe_designer
```

## 1. Docker (preferred)

The image is multi-arch (linux/amd64, linux/arm64) and contains ODT plus BLAST, bedtools,
bowtie/bowtie2, samtools, bcftools.

```bash
docker pull ghcr.io/helmholtzai-consultants-munich/odt:2.3.0
docker run --rm ghcr.io/helmholtzai-consultants-munich/odt:2.3.0 genomic_region_generator --help
```

Running with the project directory mounted at the same path (so relative and absolute paths in the
config keep working) and as the current user (so output files are not owned by root):

```bash
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PROJECT_DIR":"$PROJECT_DIR" -w "$PROJECT_DIR" \
  ghcr.io/helmholtzai-consultants-munich/odt:2.3.0 \
  scrinshot_probe_designer -c scrinshot_probe_designer.yaml
```

Notes:
- Docker Desktop on macOS limits RAM/CPU per VM (Settings → Resources). Make sure it is at least
  `n_jobs` CPUs and ~8–16 GB RAM for whole-transcriptome references.
- If `docker info` fails, the daemon is not running (start Docker Desktop / `sudo systemctl start docker`).

## 2. Apptainer / Singularity (HPC)

```bash
apptainer build odt_2.3.0.sif docker://ghcr.io/helmholtzai-consultants-munich/odt:2.3.0
apptainer exec odt_2.3.0.sif genomic_region_generator --help
```

(`singularity` works the same way.) Building needs ~3 GB in `$APPTAINER_CACHEDIR` / `$TMPDIR`; on
clusters set these to scratch if the home quota is small. Apptainer binds `$HOME`, `$PWD` and `/tmp`
automatically; other directories need `--bind /path`.

Running:

```bash
apptainer exec --bind "$PROJECT_DIR" odt_2.3.0.sif scrinshot_probe_designer -c scrinshot_probe_designer.yaml
```

For an interactive shell inside the container, micromamba has to be activated:

```bash
apptainer shell odt_2.3.0.sif
eval "$(micromamba shell hook --shell bash)"
micromamba activate
```

## 3. conda / mamba

Use Miniforge if the institution does not allow Anaconda.

```bash
bash <skill_dir>/scripts/get_odt_file.sh environment.yaml ./odt_environment.yaml
conda env create -n odt -f ./odt_environment.yaml
conda activate odt
pip install "oligo-designer-toolsuite==2.3.0"
genomic_region_generator --help
```

`environment.yaml` (from ODT v2.3.0) pins: python 3.12, blast ≥ 2.15, bedtools ≥ 2.30, bowtie ≥ 1.3.1,
bowtie2 ≥ 2.5, samtools ≥ 1.22, bcftools ≥ 1.22, channels bioconda + conda-forge.

Running non-interactively (e.g. from a script):

```bash
conda run -n odt --no-capture-output scrinshot_probe_designer -c scrinshot_probe_designer.yaml
```

## Verifying the version

- Docker/Apptainer: the tag/filename already pins it; additionally
  `python -c "import oligo_designer_toolsuite as o; print(o.__version__)"` inside the container.
- conda: `conda run -n odt pip show oligo-designer-toolsuite | grep Version`.
