# Data Resources

## Folder `genomic_regions`

The file "sequences_ncbi_exon_exon_junctions.fna" was derived from "data/genomic_regions/exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna" to contain only sequences for *ABAT* via:

```
grep ">ABAT" data/genomic_regions/exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna -A 1 > tests/data/genomic_regions/sequences_ABAT_ncbi_exon_exon_junctions.fna
```

The file "sequences_AARS1_ncbi_exon_exon_junctions_short.fna" was generated with the `genomic_region_generator` pipeline and the config file `tests/data/configs/genomic_region_generator_ncbi_exon_exon_junctions_short.yaml`. The output of the pipeline was filtered for AARS1 with the following command:

```
grep ">AARS1" exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.p14.fna -A 1 > sequences_AARS1_ncbi_exon_exon_junctions_short.fna
```

The file "exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.p14.fna" was subsequently deleted.

## Folder ```annotations```

The file "GCF_000001405.40_GRCh38.p14_genomic.gtf" was downloaded from NCBI via the `NcbiGenomicRegionGenerator`.
Then, chromsomes 16 and KI270728.1 were extracted via:

```
awk '$1 ~ /^#/ {print $0;next} {if ($1 == "16" || $1 == "KI270728.1") print}' GCF_000001405.40_GRCh38.p14_genomic.gtf > custom_GCF_000001405.40_GRCh38.p14_genomic_chr16_KI270728-1.gtf
```
