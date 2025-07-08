# Data Resources

## Folder `genomic_regions`

The file "sequences_ncbi_exon_exon_junctions.fna" was derived from "data/genomic_regions/exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna" to contain only sequences for *ABAT* via:

```
grep ">ABAT" data/genomic_regions/exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna -A 1 > tests/data/genomic_regions/sequences_ABAT_ncbi_exon_exon_junctions.fna
```
