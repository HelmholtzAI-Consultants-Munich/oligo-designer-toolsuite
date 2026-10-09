---
name: odt-probe-design
description: Design FISH/spatial transcriptomics probes (SCRINSHOT, MERFISH, SeqFISH+, HCR, CycleHCR, Oligo-Seq) with Oligo Designer Toolsuite (ODT) v2.3.0: install, inputs, parameters, config + run script.
---

# ODT probe design

This skill walks a user through designing probes with the
[Oligo Designer Toolsuite](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite)
(ODT), pinned to **v2.3.0**. The deliverables are:

1. one or two YAML config files (genomic region generator and/or probe designer), and
2. a shell script `run_odt.sh` that the user executes themselves.

Claude/Codex/ChatGPT does **not** run the (long, resource-heavy) design pipelines itself unless the user explicitly asks.

## Bundled resources

| Path | Use |
| --- | --- |
| `scripts/check_odt.sh` | Detects Docker / Apptainer / conda / native ODT and reports versions. Run it first. |
| `scripts/predict_genomic_region_paths.py` | For `source: custom`: prints the FASTA paths the genomic region generator will write, so the probe config can list them directly. Needs PyYAML. |
| `scripts/resolve_genomic_region_paths.py` | For NCBI/Ensembl: replaces `@GRG:<region>@` placeholders in the probe config with the FASTA files the generator wrote. Called from `run_odt.sh`. |
| `scripts/get_odt_file.sh` | Gets a file from the ODT repo: config templates from `data/configs/`, example initiator/readout tables from `data/functional_oligos/`, or `environment.yaml`. Uses the local checkout the skill lives in, otherwise GitHub at tag v2.3.0. **Always start from these templates**, never write a config from scratch. |
| `scripts/make_hcr_codebook.py` | Writes an HCR codebook (one initiator per gene) from the gene list and the initiator table. |
| `assets/run_odt_template.sh` | Template for the run script. |
| `references/installation.md` | Install commands for each runtime. |
| `references/genomic_regions.md` | Genomic region generator: sources, parameters, which regions each pipeline needs. |
| `references/pipelines.md` | Per pipeline: command, template, required text inputs, important numeric parameters, outputs. |
| `references/parameter_guide.md` | Background for discussing parameters (Tm, GC, BLAST specificity, set selection, runtime). |

Read the reference file for a step when you reach that step; do not load all of them up front.

## Execution environment

First determine where you are running:

- **Local mode** (Claude Code, Codex CLI, or any agent with a shell on the user's machine): follow
  all steps as written — run `check_odt.sh`, check the user's files, write the config and script
  into their project directory.
- **Sandbox mode** (claude.ai / Claude Desktop chat, ChatGPT, API code execution): you run in a
  cloud container, not on the user's computer. Signs: no access to the user's paths, no `docker`,
  a generic home directory such as `/home/user` or `/mnt/data`. Tell the user this once, then:
  - Step 2: skip `check_odt.sh`. Ask how ODT is installed (Docker / Apptainer / conda / not yet) and
    give the install commands from `references/installation.md` for them to run.
  - Steps 3–4: don't check the paths on disk. Ask the user to confirm that the files exist, and to
    paste the output of `grep -m3 '>' <fasta>` if the FASTA format is in doubt. Uploaded files (gene
    lists, tables, codebooks) can be checked in the sandbox.
  - Templates and example tables: `get_odt_file.sh` needs access to raw.githubusercontent.com. If the
    sandbox has no network, ask the user to upload the file from `data/configs/` or
    `data/functional_oligos/` of the ODT repo (tag v2.3.0). `make_hcr_codebook.py` runs in the sandbox
    on uploaded files; offer the codebook as a download.
  - Step 3 (regions to be generated): `predict_genomic_region_paths.py` (custom source) runs in the
    sandbox if PyYAML is available; otherwise build the paths from the pattern in `genomic_regions.md`.
  - Step 6: write the config(s), `run_odt.sh` and, if placeholders are used,
    `resolve_genomic_region_paths.py` to the sandbox output directory and offer them as downloads.
    Leave `PROJECT_DIR` as a value the user must fill in, and say so clearly. For Docker/Apptainer give the
    same `EXTRA_MOUNTS` warning as in Step 6.

## Workflow

Keep a running summary of the user's answers (runtime, pipeline, paths, parameters). Use
`AskUserQuestion` for multiple-choice decisions and plain conversation for free-text values.
Never guess a path or sequence: if a required text value is missing, ask.

### Step 1 — Which pipeline?

If not already clear from the request, ask which experiment the probes are for:
SCRINSHOT, MERFISH, SeqFISH+, HCR, CycleHCR, or Oligo-Seq. If the user is unsure, briefly
describe the options using `references/pipelines.md` and recommend one. Also ask for the
project directory where configs, inputs and outputs should live (default: current directory).

### Step 2 — Check / install ODT

Run `bash <skill_dir>/scripts/check_odt.sh` and interpret the output.

- ODT v2.3.0 already usable (image pulled, `.sif` present, or conda env with `oligo-designer-toolsuite==2.3.0`)
  → confirm with the user which runtime to use and continue.
- Not installed → offer the options **in this order of preference** and let the user pick:
  1. **Docker** (if `docker` is available and the daemon runs) — `ghcr.io/helmholtzai-consultants-munich/odt:2.3.0`
  2. **Apptainer/Singularity** (typical on HPC, where Docker is unavailable) — build a `.sif` from the Docker image
  3. **conda/mamba** — create env from ODT's `environment.yaml` (via `get_odt_file.sh environment.yaml ...`) and `pip install oligo-designer-toolsuite==2.3.0`
- If none of docker/apptainer/conda exists, tell the user what to install (Docker Desktop, or Miniforge) and stop until they have it.

Show the exact install commands from `references/installation.md` and ask before running them
(they download several GB). Once installed, verify with `genomic_region_generator --help`
through the chosen runtime.

### Step 3 — Genomic region files

Ask: **"Do you already have the genomic region FASTA files, or do they need to be generated?"**

- **Already prepared** → ask for the paths of
  - the FASTA file(s) to design probes from (target database), and
  - the FASTA file(s) used as reference for the specificity filters,
  and check they exist (`ls`, and `grep -m3 '>' file` to show the header format). Mention the
  headers must be ODT-style (`>GENE::...`); FASTAs from other tools will not work. If only a
  GTF + genome FASTA exist, that is the "needs to be generated" case with `source: custom`.
- **Needs to be generated** → follow `references/genomic_regions.md`:
  1. Ask for the annotation source: NCBI, Ensembl, or custom GTF + genome FASTA.
  2. Collect the text inputs for that source (species, release, taxon / accession, or file paths).
  3. Pick the genomic regions the chosen pipeline needs (table in `genomic_regions.md`); tell the user which you selected and why.
  4. Write `genomic_region_generator.yaml` from the matching template. Unused `source_params` must be `null`.
  5. Refer to the not-yet-existing FASTA files in the probe designer config (details and examples in
     "Referring to generated files" in `genomic_regions.md`):
     - **NCBI / Ensembl** → placeholders `@GRG:exon@`, `@GRG:exon_exon_junction@`, `@GRG:gene@`,
       `@GRG:cds@`, `@GRG:utr@`, `@GRG:intron@`, `@GRG:intergenic@`. `run_odt.sh` replaces them after the
       generator has run and writes `<config>.resolved.yaml`, which is the config that actually runs.
     - **custom** → the real paths, from
       `python3 <skill_dir>/scripts/predict_genomic_region_paths.py genomic_region_generator.yaml`
       (run in the project directory; use the ODT runtime's Python if the host has no PyYAML).
  6. **Tell the user** (and put it as a YAML comment above the FASTA lists in the probe config):
     - the files are created by the genomic region generator, which has to finish successfully before
       the probe designer can run (`run_odt.sh` runs it first and checks the files exist);
     - for NCBI/Ensembl, why placeholders are used: the file names contain the annotation release and
       genome assembly, which ODT only learns while downloading (e.g. `annotation_release: current`
       becomes the actual release number, NCBI `110` becomes `NCBI_Homo_sapiens_Annotation_Release_110`),
       so they cannot be written in advance. The resolved config with the real paths is
       `<config>.resolved.yaml` after the run.

### Step 4 — Minimal required inputs

From `references/pipelines.md`, collect **every text-valued parameter** (paths, IDs, sequences,
names) of the chosen pipeline that is not already known. Ask for them together in one message
as a short numbered list, with the template default shown where one exists (e.g. primer
sequences, channel names), so the user can answer "keep default". At minimum:

- `dir_output` (suggest `output_<pipeline>`),
- the gene list file (`file_regions` / `file_region_ids`): one gene ID per line, IDs must match
  the FASTA headers (NCBI → gene symbols like `AARS1`; Ensembl → `ENSG...` IDs). If the user gives
  gene names inline, offer to write the file for them. Empty/`null` means *all genes* — warn that this is very expensive.
- pipeline-specific files/sequences (readout/initiator tables, codebook, primers, linker, channels, VCF, ...).
  For HCR and CycleHCR, follow "Example tables and codebooks" below.

Also ask for `n_jobs` (CPU cores) — default to the number of cores available minus one, max ~16 unless they have a big machine.

#### Example tables and codebooks (HCR, CycleHCR)

ODT ships example tables in `data/functional_oligos/` (details in `references/pipelines.md`). Ask
whether the user has their own table; if not, offer the example, but **only use it after the user
confirms the sequences match the reagents they will order** (initiators must match their HCR
amplifier hairpins; readout probes must match the fluorescent readouts they use). Copy the
chosen file into the project directory with `get_odt_file.sh` so the config does not depend on the
ODT repo, and point the config at the copy.

- **HCR**: `file_initiator_table` → `data/functional_oligos/hcr_initiators.tsv` (5 initiators,
  `bit_1`–`bit_5`). A codebook is **required** in v2.3.0. If the user has none, ask which gene should
  get which initiator (or whether the order of the gene list is fine), then run
  `python3 <skill_dir>/scripts/make_hcr_codebook.py --genes <gene_list> --initiators <initiator_table> --out <project_dir>/hcr_codebook.tsv [--assign GENE=bit_k ...]`
  and set `file_codebook` to the output. More genes than initiators → the script stops; discuss
  adding initiators or reducing genes (`--allow-shared-bits` only if pooling genes into one channel
  is intended).
- **CycleHCR**: `file_readout_probe_table` → `data/functional_oligos/cycle_hcr_readout_probes_no_codebook.tsv`
  (180 readout probes) with `file_codebook: null`, so the pipeline generates the codebook. If the user
  brings their own codebook, the readout table needs a matching `bit` column (format as in
  `cycle_hcr_readout_probes_with_codebook.tsv`).
- The example codebooks (`hcr_codebook.tsv`, `cycle_hcr_codebook.tsv`) only contain the demo genes
  AARS1, DECR2, PRR35 — never use them for a real gene list; they are format examples only.

### Step 5 — Optional parameters and discussion

Show a compact table of the **user-facing numeric/boolean parameters** of the chosen pipeline
(the "USER PARAMETERS" block of the template, or for Oligo-Seq the `target_probe` section) with
their default values, grouped (probe length, GC, Tm, homopolymers, set selection, weights, ...).
Then ask: *"Do you want to change any of these, or discuss what values suit your experiment?"*

- Answer questions using `references/parameter_guide.md` and the comments in the template. Be
  concrete: explain the trade-off, suggest a value, and say which default the published protocol uses.
- Point out consistency constraints (e.g. `min <= opt <= max`, `set_size_min <= set_size_opt`,
  exon-exon-junction `block_size` vs. maximum probe length, HCR `L + gap + R` length).
- "Developer parameters" (BLAST, Tm tables, graph search) should stay at defaults unless the user
  asks or has a specific reason (e.g. hybridization buffer with different salt/formamide).
- Loop until the user says they are done.

### Step 6 — Write the config and the run script

1. Get the template for the pipeline into the project directory with
   `bash <skill_dir>/scripts/get_odt_file.sh data/configs/<template>.yaml <project_dir>/<pipeline>.yaml`
   (template names: see `references/pipelines.md` and `references/genomic_regions.md`) and edit
   only the values that change. Keep all keys and comments (the pipelines read keys directly, a
   missing key is a crash). For Oligo-Seq keep `schema_version: 2`, and disable filters whose inputs
   the user does not have (e.g. `variant_filter.enabled: false` without a VCF).
2. Paths: prefer paths relative to the project directory. With Docker/Apptainer the container only
   sees the project directory and the directories in `EXTRA_MOUNTS` (Apptainer additionally binds
   `$HOME`, `/tmp` and the current directory by default). Fill `EXTRA_MOUNTS` with the directories of
   input files you know are elsewhere, and **warn the user** (in particular with Docker) that they
   must add every directory outside these default paths to `EXTRA_MOUNTS`: absolute paths elsewhere,
   `../` paths, and targets of symlinks inside the project (a link `data/genome.fna -> /data/hg38/...`
   needs `/data/hg38`). Otherwise ODT inside the container reports the file as not found.
3. Create `run_odt.sh` from `assets/run_odt_template.sh`: fill in runtime, image/SIF/env, project
   directory, and enable the genomic-region step only if Step 3 required it. If placeholders are
   used, copy `scripts/resolve_genomic_region_paths.py` next to the script. The script checks that
   every FASTA file listed in the (resolved) probe config exists before starting the probe designer.
4. Sanity-check: `python3 -c "import yaml,sys; yaml.safe_load(open(sys.argv[1]))" <file>` for each
   config (if Python+PyYAML are present on the host), `bash -n run_odt.sh`, and diff against a fresh copy of the
   template (`get_odt_file.sh` to a temp file, then `diff <template> <config>`) — show the user the diff as the summary of their choices.
5. Tell the user how to run it (`bash run_odt.sh`, or wrap in `sbatch`/`nohup` for long jobs), the
   expected runtime drivers (number of genes, `n_jobs`, size of reference FASTA), and where to
   find the results (see "Outputs" in `references/pipelines.md`).

## Rules

- Pin everything to v2.3.0 (image tag `2.3.0`, `pip install oligo-designer-toolsuite==2.3.0`).
  Do not use `latest`.
- Never invent sequences (primers, readouts, initiators, codebooks). Use template defaults only when
  the user agrees, otherwise ask the user for theirs.
- Ask before downloading or installing anything; never run the design pipeline without the user's go-ahead.
- If the user hits an error later, read the log in `<dir_output>/` (`*.log`) before suggesting fixes.
