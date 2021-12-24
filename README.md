# Probe Design

In a targeted spatial transcriptomics experiment, you need to know a priori what genes you would like to see. These genes are targeted by probes that target mRNA sequences by hybridization. However, we cannot design probes for all genes as some genes are just too similar in sequence to keep apart. A probe for one of these genes would also bind to the other, similar gene and therefore cannot give a specific readout.

This package provides a probe designer that designs all valid probes for a given gene. Probes are designed with a specific length and are pre-filtered based on their GC content and melting temperature. To ensure that the potential probes don't map to any other region on the genome, the probes are blasted with ```blastN``` against the transcriptome of the given species (e.g. human transcriptome). All probes that map to other genomic regions (outside their gene region) with more than 80% similarity are discarded, given that the probe and matching genomic regions have at least 50% overlap and the matching region covers -4 to +5nt of the probe center. All remaining probes are reported as suitable probes for a given gene. Genes can be filtered by the number of possible probes, i.e. all genes with less than *n* possible probes can be filtered out. 

## Installation

The code has been implemented using Python 3.8. To install the necessary packages for this framework with conda run:

```
conda env create -f environment.yaml
```

To install the package via pip:

```
pip install .        (Installation as python package: run inside directory)
``` 
or if you want to develop the package:
```
pip install -e .        (Installation as python package: run inside directory)
``` 

To install the package via setup.py:
```
python setup.py install
```

## Usage

### Input parameters

All input parameters are defined in a config yaml file, e.g. ```probe_design.yaml```.
The user can choose if he want to download the gene annotation from *NCBI* or *ensemble* (the user can specify the enseble version, e.g. 104). The annotation will be downloaded automatically via ftp download server and the gene *gtf* file and genome *fasta* file will be saved to the output directory. Note, annotation can be downloaded for the following species: *Homo Sapiens*. Alternatively, the user can specify a genome *fasta* and a gene *gtf* file, that should be used for the probe design. 

Probes can be designed for all genes in the gene *gtf* file, or alternatively for only a selected number of genes. If the user wants to use all genes defined in the gene *gtf* file, the parameter ```gene_list``` should be empty, i.e. ```gene_list: []```. If the user wants to specifiy a list of genes, for which probes should be designed, he should provide the gene identifiers in ```gene_list```, e.g. for ensemble gene annotation ```gene_list: ['ENSG00000000003', 'ENSG00000000938', 'ENSG00000001631', 'ENSG00000003393']```. 

To design probes, the user has to specify the length of the potential probes, e.g. 45 nt, provide the GC content and melting temperature filtering criteria, e.g. GC content should be between 40 and 60%, as well as the minimum number of probes per gene in oder to filter out genes for which no probes can be designed (if all genes should be returned set this parameter to 0). For the nucleotide blast search (```blastN```), the user has to specify three blast parameters: the word size for the ```blastN``` seed (exact match of probe to genomic region, recommended: 10), the minimum coverage between the probe and genomic region (ranging from 0 to 1 - full coverage), as well as the maximum similarity between the probe and genomic region (ranging from 0 to 1 - no missmatch).

Optional parameters that the user can define include the number of genes that should be processed in one batch (for multiprocessing) as well as different parameters for the computation of the melting temperature. 


### Pipeline Description

The pipeline has four major steps: 

1) Downloading annotations: within this step the gene and genome annotations are downloaded from *NCBI* or *ensemble*. The main functions used for this step are in ```src/load_annotations.py```. Annotations are downloaded via ftp server (```ftp_download``` function). One functions downloads the gene *gtf* (```download_gene_gtf``` function) and another funtion downloads the genome *fasta* file (```download_genome_fasta``` function). In case, *NCBI* annotation is chosen for download, the annotation has to be post-processed. The chromosome names have to be mapped from *RefSeq* accession number to either *sequence-name* (for chromosomes) or to *GenBank* accession number (for scaffolds) in order to be used by ```bedtools```. The required mapping is downloaded automatically from *NCBI* (```download_chr_mapping``` function). 

2) Get gene list and probe list per gene: Either get the full list of genes from gene *gtf* annotation (```get_gene_list``` function in ```src/utils.py```) or load the list of genes provided by the user. The main functions to retrieve the list of all potential probes per gene are implemented in ```src/get_probes.py```. The main function ```get_probes``` is called and generates batches of *n* (user-defined) genes that are processed in paralelle with multiprocessing, using all available cores on the server. Each job will run three functions consecutively: ```_get_exome_fasta```, ```_GC_Tm_filter``` and ```_get_probes_fasta```. The first function will obtain the sequence of all annotated exons of each gene and write them into a *fasta* file. The second function will create a sliding window over each exon sequence (with a window step of 1nt) to create a list of potential probes per exon, excluding probes that contain masked nucleotides (*N*). All created probes are filtered based on GC content and melting temperature. The probes passing both filters are aggregated over each gene and all probes that have exactly the same sequence (only considering probes of the same gene) are merged into one fasta entry, saving the information of exon and transcript-of-origin as well as start and end position of the respective probes. The last functions outputs all probe sequences with corresponsing headers, containing information about *gene ID*, *transcript ID*, *exon ID*, *probe ID*, *chromosome*, *start*, *end*, *strand*, *GC content* and *melting temperature* in their header. 

3) Run Blast search: for each *fasta* file with potential probe sequences run blast search against transcriptome. The main functions to run the blast search are implemented in ```src/filter_probes.py```. In a first step the transcriptome is created using the function ```_get_transcriptome_gtf``` and the blast database of the transcriptome is created using the ```Bio.Blast.Applications.NcbimakeblastdbCommandline``` api. The main function ```run_blast_search``` is called, which executes the blast search on each batch (one *fasta* file per batch was created in step 2) in paralelle with multiprocessing, using all available cores on the server. Therefore, we use the ```Bio.Blast.Applications.NcbiblastnCommandline``` api, which runs the bioconda installation of blastn with the defined parameters. Here, we hardcoded 4 threads per blast search (this can be changed in the code in the ```probe_design.py``` python script).

4) Filter probes with Blast search results: once the blast search is finished for all batches, the results are used to filter out probes that map to too similar sequences within the genome (but not within it's own gene). This part is divided into three steps (all functions are implemented in ```src/filter_probes.py```): load blast results (with ```_read_blast_output```), find probes without matches (with ```_get_probes_wo_match```) and save suitable probes per gene into outputfile (```_write_output```). We filter the probes based on the percent identity that should not be greater than the user-defined threshold, the alignment length (coverage) that should not be greater the user-defined threshold and coverage at the region around the probe center (the -4 to +5 nt region should not be covered). In a last step, we create one output entry per unique probe (defined by start and end position of the probe) and list all transcripts that contain this probe (i.e. the probe can be suitable for multiple transcripts if those transcripts contain the same exon).


### Running the Pipeline

We tested the pipeline on a 96 core CPU server. Therefore, 

- we connect to ```vicb-submit-01``` or ```vicb-submit-02``` and open a ```screen``` (this allows to run srun in background without crashing when loosing the VPN connection)
- we submit a job to the server, e.g. ```srun -p cpu_p -c 96 --mem 200 --nice=10000 -w cpusrv20 --pty /bin/bash``` 
- we activate the conda environment, e.g. ```conda activate environment```
- we run the pipeline, i.e. ```python probe_design.py -c probe_design.yaml```

If the user choses to download the gene and genome annotation from *NCBI* or *ensemble*, the downloaded files will be saved to the user-defined output directory.

Running the pipeline will create two subdirectories in the specified output directory:

- ```probes```: batches of exon *gtf* files, exon *fasta* files, probe *fasta* files and blast output files as well as the transcriptome *gtf* and *fasta* files retrieved from the provided or downloaded gene and genome annotation
- ```results```: one file with gene that were filtered out due to insufficient number of probes as well as one file per gene, containing all possible probes with *probe sequence*, *chromosome*, *start*, *end*, *strand*, *gene ID*, *transcript ID*, *exon ID*, *probe ID* (artificial ID generated for all probes of one gene), *GC content* and *melting temperature*

All steps and config parameters will be documented in a log file, that is saved in the directory where the pipeline is executed from. The logging file will have the format: ```log_probe_design_{year}-{month}-{day}-{hour}-{minute}.txt```.
