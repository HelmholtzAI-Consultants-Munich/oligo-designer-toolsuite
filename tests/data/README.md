# Data Resources

## Folder ```annotations```

The file "GCF_000001405.40_GRCh38.p14_genomic.gtf" was downloaded from NCBI via the `NcbiGenomicRegionGenerator`.
Then, chromsomes 16 and KI270728.1 were extracted via:

```
awk '$1 ~ /^#/ {print $0;next} {if ($1 == "16" || $1 == "KI270728.1") print}' GCF_000001405.40_GRCh38.p14_genomic.gtf > custom_GCF_000001405.40_GRCh38.p14_genomic_chr16_KI270728-1.gtf
```
