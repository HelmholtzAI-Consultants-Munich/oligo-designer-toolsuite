# Data Resources

## Folder ```annotations```

This folder contains example annotations from NCBI release 110 for *H.Sapiens*. It includes FASTA and GTF example files, filtered for chromosome 16 with:

```
awk '$1 ~ /^#/ {print $0;next} {if ($1 == "16") print}' annotation.gff > annotation_chr16.gtf
```

using ```awk``` or

```
seqkit grep -i -r -p '^16' genome.fna -o genome_chr16.fna
```

using ```seqkit``` tool.

## Folder ```configs```

This folder contains the configuration files for different ready-to-use genomic region generator or oligo design pipelines with default parameters that can be adjusted by the users.

## Folder ```genes```

This folder contains example lists of gene names that can be used to test run the pipelines.

## Folder ```genomic_region```

This folder contains different pre-generated genomic region files that were generated from the Fasta and GTF file in the `annotations` folder using the `_genomic_region_generator.py` pipeline.
