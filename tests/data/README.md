# Data Resources

## Folder `genomic_regions`

The file "sequences_AARS1_ncbi_exon_exon_junctions_short.fna" was generated with the `genomic_region_generator` pipeline and the config file `tests/data/configs/genomic_region_generator_ncbi_exon_exon_junctions_short.yaml`. The output of the pipeline was filtered for AARS1 with the following command:

```
grep ">AARS1" exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.p14.fna -A 1 > sequences_AARS1_ncbi_exon_exon_junctions_short.fna
```

The file "exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.p14.fna" was subsequently deleted.
