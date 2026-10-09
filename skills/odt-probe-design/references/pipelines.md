# ODT v2.3.0 probe design pipelines

All pipelines: `<command> -c <config.yaml>`. Paths in configs are resolved relative to the working
directory the command runs in. Templates are in `data/configs/` of the ODT repo; get them with `scripts/get_odt_file.sh data/configs/<template>.yaml <dest>`.

Common structure (all except Oligo-Seq): flat YAML with a "USER PARAMETERS" block and a
"DEVELOPER PARAMETERS" block. Oligo-Seq uses a nested, pydantic-validated schema (`schema_version: 2`).

"Text inputs" below are the parameters Claude/Codex/ChatGPT must collect in Step 4 (paths, IDs, sequences, names).
"Key numeric parameters" are the ones to show in Step 5. They are written as `parameter` (default);
`x_min/opt/max` (a/b/c) is short for `x_min`, `x_opt`, `x_max` with defaults a, b, c. Units: lengths in nt,
GC content and isoform consensus in %, Tm and temperatures in °C, ΔG in kcal/mol.

Common text inputs for every probe designer:

| Parameter | Notes |
| --- | --- |
| `dir_output` | Output directory. Suggest `output_<pipeline>`. |
| `file_regions` | Text file, one gene ID per line, IDs as in FASTA headers. `null` = all genes (very slow). |
| `files_fasta_target_probe_database` | List of FASTA files to design probes from (see genomic_regions.md). |
| `files_fasta_reference_database_target_probe` | List of FASTA files for the specificity filter. |

Common numeric parameters: `n_jobs` (4), `write_intermediate_steps` (true).

---

## SCRINSHOT — `scrinshot_probe_designer`

Padlock probes for SCRINSHOT (rolling-circle amplification smFISH). Each target probe is split into two
arms ligated on the target; a detection oligo (with uracils) binds the backbone.
Template: `scrinshot_probe_designer.yaml`.

**Text inputs**: common ones only.

**Key numeric parameters**
- Target probe: `target_probe_length_min/max` (40/45), `target_probe_isoform_consensus` (50),
  `target_probe_GC_content_min/opt/max` (40/50/60), `target_probe_Tm_min/opt/max` (65/70/75),
  `target_probe_homopolymeric_base_n` (A/T/C/G: 5/5/5/5).
- Padlock arms: `target_probe_padlock_arm_length_min` (10), `target_probe_padlock_arm_Tm_min/max` (50/60),
  `target_probe_padlock_arm_Tm_dif_max` (2), `target_probe_ligation_region_size` (5; 0 disables the
  ligation-site seed BLAST).
- Sets: `n_sets` (3), `set_size_min/opt` (3/5), `distance_between_target_probes` (0).
- Score weights: `target_probe_isoform_weight` (2), `target_probe_GC_weight` (1), `target_probe_Tm_weight` (1).
- Detection oligo: `detection_oligo_min_thymines` (2), `detection_oligo_length_min/max` (15/40),
  `detection_oligo_U_distance` (5), `detection_oligo_Tm_opt` (56).

**Outputs**: `padlock_probes.yml`, `padlock_probes` table (Excel, one sheet per gene),
`padlock_probes_order.yml` (padlock + detection oligo sequences).

---

## MERFISH — `merfish_probe_designer`

Encoding probes = 30 nt target region + readout sequences + PCR primers; combinatorial barcodes
(Hamming-distance code). Designs target probes, readout probes, codebook and primers.
Template: `merfish_probe_designer.yaml`.

**Text inputs** (in addition to common)
| Parameter | Default | Notes |
| --- | --- | --- |
| `files_fasta_reference_database_readout_probe` | exon + junction | Readouts must not hit the transcriptome. |
| `files_fasta_reference_database_primer` | exon + junction | |
| `reverse_primer_sequence` | `CCCTATAGTGAGTCGTATTA` | Rev. complement of the T7 promoter; change only if their protocol uses another one. |
| `channels_ids` | `[Alexa488, Cy3b, Alexa647]` | Fluorophore/channel names, list of strings. |

**Key numeric parameters**
- Target probe: `target_probe_length_min/max` (30/30), `target_probe_isoform_consensus` (50),
  `target_probe_GC_content_min/opt/max` (43/55/63), `target_probe_Tm_min/opt/max` (66/72/76),
  `target_probe_homopolymeric_base_n` (A/T/C/G: 6/6/6/6), `target_probe_T_secondary_structure` (76).
- Sets: `n_sets` (3), `set_size_min/opt` (50/50), `distance_between_target_probes` (0).
- Score weights: `target_probe_GC_weight` (1), `target_probe_Tm_weight` (1), `target_probe_isoform_weight` (2).
- Readout probe: `readout_probe_length` (20), `readout_probe_base_probabilities` (A/C/G/T: 0.25/0/0.5/0.25),
  `readout_probe_GC_content_min/max` (40/50), `readout_probe_homopolymeric_base_n` (G: 3),
  `readout_probe_set_size` (16), `readout_probe_homogeneous_properties_weights` (TmNN_oligo/GC_content_oligo: 1/1),
  `n_bits` (16), `min_hamming_dist` (4), `hamming_weight` (2).
  Constraints: `readout_probe_set_size` should equal `n_bits` (one readout per bit). The codebook is
  built greedily from all `n_bits`-long words with exactly `hamming_weight` ones and pairwise distance
  ≥ `min_hamming_dist`, so the number of genes is capped: the template default (16 bits, weight 2,
  MHD 4) only yields **8 codewords**; the classic MERFISH MHD4 code (16 bits, weight 4) allows up to
  ~140 (greedy may find fewer). Always check the gene count against this and discuss.
- Primer: `primer_length` (20), `primer_base_probabilities` (A/C/G/T: 0.25/0.25/0.25/0.25),
  `primer_GC_content_min/max` (50/65), `primer_number_GC_GCclamp` (1), `primer_number_three_prime_base_GCclamp` (2),
  `primer_homopolymeric_base_n` (A/T/C/G: 4/4/4/4), `primer_max_len_selfcomplement` (6),
  `primer_max_len_complement_reverse_primer` (5), `primer_Tm_min/max` (60/75),
  `primer_T_secondary_structure` (76), `primer_secondary_structures_threshold_deltaG` (0).

**Outputs**: `merfish_probes.yml`, `merfish_probes` table, `merfish_probes_order.yml`, plus readout/codebook info.

---

## SeqFISH+ — `seqfish_plus_probe_designer`

28 nt target probes with readout sites, pseudocolor barcoding over several rounds, primers.
Template: `seqfish_plus_probe_designer.yaml`. Targets CDS + UTR by default and prefers UTR probes
(`target_probe_UTR_weight`).

**Text inputs** (in addition to common): same as MERFISH (`files_fasta_reference_database_readout_probe`,
`files_fasta_reference_database_primer`, `reverse_primer_sequence`, `channels_ids`).

**Key numeric parameters**
- Target probe: `target_probe_length_min/max` (28/28), `target_probe_isoform_consensus` (100),
  `target_probe_GC_content_min/opt/max` (45/55/65), `target_probe_homopolymeric_base_n` (A/T/C/G: 5/5/5/5),
  `target_probe_T_secondary_structure` (76).
- Sets: `n_sets` (3), `set_size_min/opt` (24/24), `distance_between_target_probes` (2).
- Score weights: `target_probe_GC_weight` (1), `target_probe_UTR_weight` (10).
- Readout probe: `readout_probe_length` (15), `readout_probe_base_probabilities` (A/C/G/T: 0.25/0.25/0.25/0.25),
  `readout_probe_GC_content_min/max` (40/60), `readout_probe_homopolymeric_base_n` (G: 3),
  `n_barcode_rounds` (4), `n_pseudocolors` (20).
  Capacity: up to `n_channels × n_pseudocolors^(n_barcode_rounds − 1)` genes (last round is a parity
  check); readouts needed = channels × pseudocolors × rounds.
- Primer: same parameters and defaults as MERFISH (here `primer_secondary_structures_threshold_deltaG` (0)
  is in the developer block).

**Outputs**: `seqfish_plus_probes.yml`, table, `seqfish_plus_probes_order.yml`.

---

## HCR — `hcr_probe_designer`

HCR v3 split-initiator probe pairs: L and R halves (25 nt each, 2 nt gap) each carrying half of an
initiator. Template: `hcr_probe_designer.yaml`.

**Text inputs** (in addition to common)
| Parameter | Default | Notes |
| --- | --- | --- |
| `file_initiator_table` | `data/functional_oligos/hcr_initiators.tsv` | TSV/CSV with columns `bit`, `initiator_L_sequence`, `initiator_R_sequence` (one row per hairpin amplifier). The user's own table, or the example below after they confirm it matches their hairpins. |
| `file_codebook` | `data/functional_oligos/hcr_codebook.tsv` | **Required** (in v2.3.0, `null` aborts with "Generation of codebook is not yet implemented"). Columns `region_id`, `bit_1..bit_n`, 0/1 entries; region IDs = gene IDs, and every gene in the gene list must be in it. |
| `linker_sequence` | `AA` | Spacer between initiator half and target-binding half. |

Example files in the ODT repo (`data/functional_oligos/`, get with `scripts/get_odt_file.sh`):
- `hcr_initiators.tsv`: 5 initiator pairs `bit_1`–`bit_5`. Reusable if the user confirms they match
  their amplifier hairpins. First row:
  ```
  bit	initiator_L_sequence	initiator_R_sequence
  bit_1	gAggAgggCAgCAAACgggAAgAgTCTTCCTTTACg	gCATTCTTTCTTgAggAgggCAgCAAACgggAAgAg
  ```
- `hcr_codebook.tsv`: demo genes only (format example). For a real gene list create the codebook with
  `scripts/make_hcr_codebook.py --genes <gene_list> --initiators <initiator_table> --out hcr_codebook.tsv`
  (one initiator per gene, `--assign GENE=bit_k` to choose; at most as many genes as initiators).
  ODT keeps only rows with exactly one `1`; genes missing from the codebook abort the run.

**Key numeric parameters**
- General: `top_n_sets` (1; number of sets written to `hcr_probes.yml` and `hcr_probes_order.yml`).
- Target probe: `target_probe_isoform_consensus` (0), `target_probe_L/R_probe_sequence_length` (25/25),
  `target_probe_gap_sequence_length` (2), `target_probe_GC_content_min/max` (40/65), `target_probe_Tm_min/max` (50/80),
  `target_probe_homopolymeric_base_n` (A/T/C/G: 4/4/4/4), `target_probe_T_secondary_structure` (37),
  `target_probe_junction_region_size` (0).
- Sets: `n_sets` (30), `set_size_min/opt` (10/25), `distance_between_target_probes` (2).
- Score weights: `target_probe_isoform_weight` (1).
- Developer: `target_probe_secondary_structures_threshold_deltaG` (-7).

**Outputs**: `hcr_probes.yml`, table, `hcr_probes_order.yml`.

---

## CycleHCR — `cycle_hcr_probe_designer`

Long L/R probe pairs (45 + 2 + 45 nt) with readout probe binding sites for cyclic HCR imaging.
Template: `cycle_hcr_probe_designer.yaml`.

**Text inputs** (in addition to common)
| Parameter | Default | Notes |
| --- | --- | --- |
| `file_readout_probe_table` | `data/functional_oligos/cycle_hcr_readout_probes_with_codebook.tsv` | Columns `bit`, `channel`, `readout_probe_id`, `L/R`, `readout_probe_sequence`. Without a codebook, the `bit` column is omitted. |
| `file_codebook` | `data/functional_oligos/cycle_hcr_codebook.tsv` | Optional: `null` = codebook generated from the readout table. Columns `region_id`, `bit_1..bit_n`. |
| `linker_sequence` | `TT` | |
| `forward_primer_sequence` | `TAATACGACTCACTATAGCGTCATC` | T7 promoter-based. |
| `reverse_primer_sequence` | `CGACACCGAACGTGCGACAA` | |

Example files in the ODT repo (`data/functional_oligos/`, get with `scripts/get_odt_file.sh`):
- `cycle_hcr_readout_probes_no_codebook.tsv`: 180 readout probes (columns `channel`, `readout_probe_id`,
  `readout_probe_sequence`, `L/R`), use with `file_codebook: null` so ODT generates the codebook (bits
  are assigned from the sorted `readout_probe_id`/`channel`). Reusable if the user confirms the readouts.
- `cycle_hcr_readout_probes_with_codebook.tsv` + `cycle_hcr_codebook.tsv`: demo pair for the genes
  AARS1, DECR2, PRR35 (format example for users who bring their own codebook).

**Key numeric parameters**
- Target probe: `target_probe_isoform_consensus` (0), `target_probe_L/R_probe_sequence_length` (45/45),
  `target_probe_gap_sequence_length` (2), `target_probe_GC_content_min/max` (30/90),
  `target_probe_Tm_min/max` (90/200; effectively off), `target_probe_homopolymeric_base_n` (A/T/C/G: 6/6/6/6),
  `target_probe_T_secondary_structure` (90), `target_probe_junction_region_size` (13).
- Sets: `n_sets` (3), `set_size_min/opt` (10/25), `distance_between_target_probes` (2).
- Score weights: `target_probe_Tm_weight` (1), `target_probe_isoform_weight` (10).

**Outputs**: `cyclehcr_probes.yml`, table, `cyclehcr_probes_order.yml`.

---

## Oligo-Seq — `oligo_seq_probe_designer`

26–30 nt probes for sequencing-based readout. Nested config, validated with pydantic — a wrong key
or type aborts with "Invalid configuration file". Every filter has `enabled: true/false`.
Template: `oligo_seq_probe_designer.yaml`. **Keep `schema_version: 2`.**

**Text inputs**
| Key path | Notes |
| --- | --- |
| `general.dir_output` | |
| `target_probe.oligo_generation.file_region_ids` | gene list file (`null` = all genes) |
| `target_probe.oligo_generation.files_fasta_probe_database` | target FASTA list |
| `target_probe.specificity_filters.specificity_blastn_filter.files_fasta_reference_database` | reference FASTA list |
| `target_probe.specificity_filters.variant_filter.files_vcf_reference_database` | VCF list (e.g. dbSNP/gnomAD for the same assembly, chromosome names must match the FASTA). No VCF → `variant_filter.enabled: false`. |
| `target_probe.specificity_filters.variant_filter.action` | `flag` or `filter` |
| `target_probe.property_filters.prohibited_sequences_filter.prohibited_sequences` | list of forbidden motifs (default `["TCT", "CTC"]`), e.g. restriction sites |
| `target_probe.property_filters.targeted_exons_filter.targeted_exons` | exon numbers as strings, only if the filter is enabled |

**Key numeric parameters** (key paths below `target_probe.`, except `general.*`)
- General: `general.n_jobs` (4), `general.write_intermediate_steps` (true).
- Generation: `oligo_generation.probe_length_min/max` (26/30), `oligo_generation.probe_split_region` (4).
- Filter switches: `property_filters.<filter>.enabled` and `specificity_filters.<filter>.enabled` (all true,
  except `targeted_exons_filter.enabled` (false)).
- Property filters: `isoform_consensus_filter.isoform_consensus` (0),
  `homopolymeric_runs_filter.homopolymeric_base_n` (A/T/C/G: 6/6/6/6), `GC_content_filter.GC_content_min/max` (45/65),
  `prohibited_sequences_filter.kmer_abundance_threshold` (k = 3/4/5/6: 0.0469/0.0117/0.0029/0.00073),
  `self_complementarity_filter.max_len_selfcomplement` (10), `Tm_filter.Tm_min/max` (50/70),
  `secondary_structure_filter.T` (37), `secondary_structure_filter.thr_DG` (0)
  (all under `property_filters.`).
- Specificity filters: `read_length_bias_filter.read_length_bias` (20),
  `cross_hybridization_blastn_filter.hit_parameters.coverage` (50),
  `specificity_blastn_filter.hit_parameters.coverage` (50) (all under `specificity_filters.`).
- Sets: `probe_set_selection.independent_set_selection.n_sets` (3), `...independent_set_selection.set_size_min/opt` (3/5),
  `...independent_set_selection.distance_between_target_probes` (0).
- Score weights (under `probe_set_selection.`): `uniform_distance_score.weight` (1), `isoform_consensus_score.weight` (1),
  `targeted_exons_score.weight` (0), `GC_content_score.weight` (1), `GC_content_score.GC_content_min/opt/max` (45/55/65),
  `Tm_score.weight` (1), `Tm_score.Tm_min/opt/max` (50/60/70).
- Tm buffer: `global_parameters.Tm_parameters.Na` (1000 mM), `global_parameters.Tm_chem_correction_parameters.parameters.fmd` (20 %).

Note: the GC/Tm min/max in `probe_set_selection` should match the corresponding filters.

**Outputs**: `oligo_seq_probes.yml`, table, `oligo_seq_probes_order.yml`.

---

## Outputs common to all pipelines

In `dir_output/`: a log file, `regions_with_insufficient_oligos_for_<db>.txt` (genes dropped and at
which step — the first place to look if genes are missing), and, if `write_intermediate_steps: true`,
per-step intermediate databases (useful to see which filter removed most candidates, but large).
