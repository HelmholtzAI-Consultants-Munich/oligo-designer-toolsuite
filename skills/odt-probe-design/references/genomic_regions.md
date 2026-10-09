# Genomic region generator

Command: `genomic_region_generator -c genomic_region_generator.yaml`

Downloads (NCBI/Ensembl) or reads (custom) a genome FASTA + GTF and writes one FASTA per requested
region type to `<dir_output>/annotation/`. File names look like

```
source__NCBI__species__Homo_sapiens__annotation_release__110__genome_assemly__GRCh38__annotation_type__exon.fna
...__annotation_type__exon_exon_junction__block_size__50.fna
...__annotation_type__utr__five_prime__three_prime.fna
```

FASTA headers are `>{gene_id}::source=...;regiontype=...;gene_id=...;transcript_id=...;exon_number=...::chr:start-end(strand)`.
The `gene_id` (NCBI: gene symbol, Ensembl: `ENSG...` ID, custom: whatever the GTF `gene_id` is) is
what the gene list file of the probe designers must contain.

Downloads are large (human genome ≈ 1 GB compressed); runtime is minutes to ~1 h, mostly I/O.

## Templates (in `data/configs/` of the ODT repo, get them with `scripts/get_odt_file.sh`)

| Source | Template |
| --- | --- |
| NCBI | `data/configs/genomic_region_generator_ncbi.yaml` |
| Ensembl | `data/configs/genomic_region_generator_ensemble.yaml` |
| custom GTF + FASTA | `data/configs/genomic_region_generator_custom.yaml` |

## Text inputs per source

Common: `dir_output` (e.g. `output_genomic_region_generator`).

**NCBI** (`source: ncbi`) — `source_params.mode` is `species` or `assembly`. **All parameters not used by the chosen mode must be `null`.**

| mode | parameter | example / notes |
| --- | --- | --- |
| species | `taxon` | `vertebrate_mammalian`; one of archaea, bacteria, fungi, invertebrate, metagenomes, plants, protozoa, unknown, vertebrate_mammalian, vertebrate_other, viral |
| species | `species` | NCBI FTP folder name, e.g. `Homo_sapiens`, `Mus_musculus` |
| species | `annotation_release` | e.g. `110` (only with `assembly_source: annotation_releases`/`auto`), otherwise `current` |
| species | `assembly_source` | `auto` (default), `annotation_releases`, `latest_assembly_versions`, `reference` |
| assembly | `refseq_assembly_accession` | e.g. `GCF_000001405.40` |
| assembly | `assembly_name` | e.g. `GRCh38.p14` |

For mode `assembly`, set `taxon`, `species`, `annotation_release`, `assembly_source` to `null`;
for mode `species`, set `refseq_assembly_accession`, `assembly_name` to `null`.

**Ensembl** (`source: ensembl`)

| parameter | example / notes |
| --- | --- |
| `species` | lower-case FTP name, e.g. `homo_sapiens`, `mus_musculus` (see ftp.ensembl.org/pub/current_gtf/) |
| `annotation_release` | `current` or a release number like `112` |

**Custom** (`source: custom`)

| parameter | required | notes |
| --- | --- | --- |
| `file_annotation` | yes | GTF with `gene_id`, `transcript_id`, `exon_number` attributes and exon (and CDS/UTR if needed) features |
| `file_sequence` | yes | genome FASTA whose sequence names match the GTF seqnames |
| `files_source` | no | e.g. `NCBI`, `Ensembl`, `GENCODE` — only used in file names/headers |
| `species` | no | |
| `annotation_release` | no | |
| `genome_assembly` | no | |

## Region types and what each pipeline needs

Set the needed ones to `true` under `genomic_regions:`. `exon_exon_junction_block_size` is the number
of bases taken on each side of a junction (so junction sequences are `2 × block_size` long). Setting
it larger than the maximum probe length makes no sense; the default 50 works for all pipelines.

| Pipeline | Target database (`files_fasta_target_probe_database`) | Reference database(s) for specificity | Regions to generate |
| --- | --- | --- | --- |
| SCRINSHOT | exon + exon_exon_junction | exon + exon_exon_junction | exon, exon_exon_junction |
| MERFISH | exon + exon_exon_junction | exon + exon_exon_junction (target, readout, primer) | exon, exon_exon_junction |
| SeqFISH+ | cds + utr | exon + exon_exon_junction | cds, utr, exon, exon_exon_junction |
| HCR | exon + exon_exon_junction | gene | exon, exon_exon_junction, gene |
| CycleHCR | exon + exon_exon_junction | gene | exon, exon_exon_junction, gene |
| Oligo-Seq | exon + exon_exon_junction | exon + exon_exon_junction | exon, exon_exon_junction |

These mirror the v2.3.0 example configs. Alternatives worth discussing:
- Targeting nascent/intronic RNA → add `intron` (or `gene`) to the target database.
- Whole-genome specificity (also DNA off-targets) → use `gene` + `intergenic` as reference, and set
  the BLAST `-strand` to `both` (see parameter guide).

## Referring to generated files in the probe config

When the generator still has to run, its output files do not exist yet when the probe config is
written. All files are written to `<dir_output>/annotation/` and named

```
source__<source>__species__<species>__annotation_release__<release>__genome_assemly__<assembly>__annotation_type__<type>.fna
```

with `<type>` one of `gene`, `intergenic`, `exon`, `intron`, `cds`, `utr__five_prime__three_prime`,
`exon_exon_junction__block_size__<exon_exon_junction_block_size>` (note ODT's spelling `assemly`).

### NCBI / Ensembl: placeholders

Release and assembly are resolved on the FTP server at download time (e.g. Ensembl `current` → `116`,
NCBI `110` → `NCBI_Homo_sapiens_Annotation_Release_110` with assembly `GRCh38.p14`), so the names cannot
be written in advance. Use placeholders and explain this to the user:

```yaml
# The FASTA files below are created by the genomic region generator (genomic_region_generator.yaml).
# @GRG:<region>@ are placeholders because the file names contain the annotation release and genome
# assembly, which are only known after the download. run_odt.sh replaces them after the generator has
# finished and runs the resolved copy <this config>.resolved.yaml.
files_fasta_target_probe_database:
  - "@GRG:exon@"
  - "@GRG:exon_exon_junction@"
```

`scripts/resolve_genomic_region_paths.py` replaces each `@GRG:<type>@` with the single matching file in
`<grg_dir_output>/annotation/` (types: gene, intergenic, exon, exon_exon_junction, cds, utr, intron)
and fails loudly if zero or several files match.

### custom: real paths

For `source: custom` the name is fully determined by the config: `<source>` = `files_source` (default
`custom`), and `species`, `annotation_release`, `genome_assembly` from the config (default `unknown`
each). Write the real paths, computed with

```bash
python3 <skill_dir>/scripts/predict_genomic_region_paths.py genomic_region_generator.yaml
```

(prints `<region>\t<path>` per enabled region; needs PyYAML, e.g. via the ODT runtime). Add a comment:

```yaml
# Created by the genomic region generator (genomic_region_generator.yaml): run it first.
files_fasta_target_probe_database:
  - output_genomic_region_generator/annotation/source__NCBI__species__Homo_sapiens__annotation_release__110__genome_assemly__GRCh38__annotation_type__exon.fna
```

Recompute the paths if the generator config changes (regions, block size, metadata fields).
