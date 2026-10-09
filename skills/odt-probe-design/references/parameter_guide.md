# Parameter guide (for discussing choices with the user)

General stance: the template defaults reproduce the published protocols (SCRINSHOT, MERFISH,
SeqFISH+, HCR v3, ...). Recommend changing them only for a concrete reason: a different hybridization
buffer, too few probes per gene, a hard-to-target gene family, a different readout chemistry, or runtime.
When the user changes something, explain the trade-off in one or two sentences and suggest a value.

## How the pipelines work (to explain filters and failures)

1. **Oligo generation** — sliding window over the target FASTA for every length in `[length_min, length_max]`.
2. **Property filters** — sequence-intrinsic: length, hard/soft-masked bases, homopolymer runs,
   GC content, Tm, secondary structure (ΔG at temperature T), self-complementarity, prohibited motifs,
   isoform consensus.
3. **Specificity filters** — BLASTN against the reference FASTA(s): probes with off-target hits are
   removed; cross-hybridization between probes of different genes is removed; optional seed-region
   checks around the ligation/junction site; Oligo-Seq adds read-length bias and variant filters.
4. **Scoring + set selection** — each probe gets an efficiency score (weighted distance of GC/Tm from the
   optimum, isoform consensus, UTR, ...); a graph of non-overlapping probes is searched for cliques to pick
   `n_sets` sets of size `set_size_opt` (fallback down to `set_size_min`).
5. **Assembly** — add readouts / initiators / primers / backbone and write order files.

Genes with fewer than `set_size_min` surviving probes are dropped and listed in
`regions_with_insufficient_oligos_for_*.txt` together with the step that removed them.

## Common knobs and their trade-offs

**Probe length** — longer → higher Tm and binding stability, more sequence available per probe but
fewer non-overlapping positions per transcript and higher cost. Keep within the chemistry's constraints
(padlock arms, split probes, synthesis length limits of the oligo pool vendor).

**GC content** (min/opt/max, %) — 40–60 % is the usual safe band. Widen (e.g. 35–65) for AT- or GC-rich
genes that lose too many probes; narrowing makes Tm more uniform across probes.

**Tm** (min/opt/max, °C) — computed with Biopython `Tm_NN` using `*_Tm_parameters` (nearest-neighbour
table, oligo concentrations `dnac1/dnac2` in nM, salts in mM) plus optional chemical correction
(formamide `fmd` in %, DMSO). Only change the buffer parameters if the user's hybridization buffer differs
from the protocol (e.g. different formamide % or salt) — then Tm min/opt/max should be re-centred
because absolute Tm values shift. Rule of thumb: ~0.6–0.7 °C lower Tm per 1 % formamide.

**Homopolymer runs** (`homopolymeric_base_n`) — a probe is rejected if it contains a run of ≥ n identical
bases. Lower = stricter (synthesis quality, non-specific binding); raise by 1 for AT-rich targets.

**Secondary structure** (`T_secondary_structure`, `secondary_structures_threshold_deltaG`) — rejects probes
whose predicted folding ΔG at temperature T is below the threshold (kcal/mol). Use the hybridization
temperature; a more negative threshold (e.g. −7 for HCR) is more permissive.

**Isoform consensus** (%) — fraction of a gene's annotated transcripts that contain the probe sequence.
100 = probes hit all isoforms (good for gene-level counts, fewer candidates); 0 = no constraint (more probes,
some isoform-specific). Increase its score weight to prefer constitutive exons without hard filtering.

**Specificity (BLAST)** — `search_parameters` are passed to `blastn` (`-perc_identity`, `-word_size`,
`-strand`, ...); `hit_parameters` define what counts as an off-target hit: `coverage` (% of probe length
aligned) or `min_alignment_length` (nt). Stricter (lower coverage / shorter min alignment / lower
perc_identity) → fewer off-targets but fewer probes; looser → more probes but more risk. For RNA targets
with transcriptome references keep `-strand minus`; for genome references (gene/intergenic) consider
`both`. Reference choice matters most: exon+junction catches mRNA off-targets, `gene` also catches
pre-mRNA/intronic ones.

**Set selection**
- `set_size_opt` / `set_size_min` — probes per gene. More probes → brighter, more robust signal; short
  or highly homologous genes may not reach it. Lower `set_size_min` to keep such genes (with weaker signal),
  or accept their loss.
- `distance_between_target_probes` — negative allows overlap, 0 adjacent, positive = min gap in nt
  (needed if neighbouring probes sterically hinder each other, e.g. HCR pairs).
- `n_sets` — alternative sets per gene (backups or for choosing). More sets → longer runtime.
- Weights — relative importance of GC/Tm/isoform/UTR in the score. Only relative values matter.
- `n_attempts_graph`, `n_attempts_clique_enum`, `diversification_fraction`, `jaccard_opt/step` —
  search effort and diversity between sets; increase only if sets are poor or too similar.

## Typical troubleshooting conversations

- **"Many genes end up in regions_with_insufficient_oligos"** → check which step removed them; relax that
  step (GC/Tm range, homopolymer, isoform consensus, BLAST coverage) or lower `set_size_min`. Paralog-rich
  genes often fail at specificity — this is usually correct behaviour.
- **"It takes too long / runs out of memory"** → fewer genes per run, smaller reference (exon+junction
  instead of gene), `write_intermediate_steps: false`, lower `n_sets`; `n_jobs` up for speed, down for memory
  (cache holds 2·n_jobs+1 regions).
- **"Different buffer"** → adjust `*_Tm_parameters` / `*_Tm_chem_correction_parameters` and re-centre Tm min/opt/max.
- **"Different species / assembly"** → regenerate genomic regions; all FASTAs (target + references) must
  come from the same annotation; gene IDs in the gene list must match the new headers.
- **"Gene not found"** → `grep -m1 '^>GENE::' target.fna`; NCBI uses symbols, Ensembl uses `ENSG` IDs.
