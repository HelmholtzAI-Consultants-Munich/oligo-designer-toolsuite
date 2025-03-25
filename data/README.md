# Data Resources

## Folder ```annotations```

This folder contains example annotations from NCBI release 110 for *H.Sapiens*. For all annotation the RefSeq-Accn was replaced with chromsome names or GenBank-Accn (e.g. from NC_000016.10 to 16). It includes FASTA, GTF and VCF example files, filtered for chromosome 16 with:

```
awk '$1 ~ /^#/ {print $0;next} {if ($1 == "16") print}' GCF_000001405.40_GRCh38.p14_genomic.gtf > custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf
```

using ```awk``` or

```
seqkit grep -i -r -p '^16' GCF_000001405.40_GRCh38.p14_genomic.fna -o custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna
```

using ```seqkit``` tool or 

```
bcftools view -r 16 --output-file custom_GCF_000001405.40.chr16.vcf --output-type v GCF_000001405.40.gz
```

using ```bcftools```. We further subselect 50.000 variants from the chromosome 16 vcf file, to reduce the file size .



## Folder ```configs```

This folder contains the configuration files for different ready-to-use genomic region generator or oligo design pipelines with default parameters that can be adjusted by the users.

## Folder ```genes```

This folder contains example lists of gene names that can be used to test run the pipelines.

## Folder ```genomic_region```

This folder contains different pre-generated genomic region files that were generated from the Fasta and GTF file in the `annotations` folder using the `_genomic_region_generator.py` pipeline.
