# Data Resources

## Folder ```annotations```

This folder contains example annotations from NCBI release 110 for *H.Sapiens* and Ensembl release 108 for *H.Sapiens*. It includes BED, fasta, GTF and GFF example files, filtered for chromosome 16 with:

```
awk '$1 ~ /^#/ {print $0;next} {if ($1 == "16") print}' annotation.gff > annotation_chr16.gtf
```

using ```awk``` or

``` 
seqkit grep -i -r -p '^16' genome.fna -o genome_chr16.fna
```

using ```seqkit``` tool.

## Folder ```configs```

This folder contains the configuration files for different ready-to-use oligo design pipelines. For each pipeline, there is a user configuration file, which allows the user to adjust pipeline specific parameters as well as a developer configuration file, which allows advanced parameter configuration that deviate from the standard parameters recommended for the pipeline. 

## Folder ```genes```

This folder contains example lists of gene names, that are used for testing and the tutorials, for NCBI and Ensembl.