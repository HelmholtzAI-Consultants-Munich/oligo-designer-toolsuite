![stability-wip](https://img.shields.io/badge/stability-work_in_progress-lightgrey.svg)

# Oligo Designer Toolsuite

Oligonucleotides (abbrev. oligos) are short, synthetic strands of DNA or RNA that have many application areas, ranging from research to disease diagnosis or therapeutics. Oligos can be used as primers during DNA amplification, as probes for *in situ* hybridization or as guide RNAs for CRISPR-based gene editing. Based on the intended application and experimental design, researchers can customize the length, sequence composition, and thermodynamic properties of the designed oligos.

Various tools exist that provide custom design of oligo sequences depending on the area of application. Interestingly, all those pipelines have many common basic processing steps, ranging from the generation of custom-length oligo sequences, the filtering of oligo sequences based on thermodynamic properties as well as the selection of an optimal set of oligos. Despite the fact that most tools apply the same basic processing steps, each newly developed tool usually uses its own implementation and different versions of package dependencies for those basic processing steps. As a consequence, the comparability of tools that differ only in certain steps is hampered, but also the development of new tools and the update of existing tools is slowed down, because developers do not have a common resource for basic functionalities to fall back on. We tackle this issue by providing such a common resource in our *Oligo Designer Toolsuite*. This Toolsuite is a collection of modules that provide all basic functionalities for custom oligo design pipelines within a flexible Python framework. All modules have a standardized I/O format and can be combined individually depending on the required processing steps.

![](docs/figures/oligo_design.png)


## Documentation

For the complete documentation of the package, working principles and toutorials refer to the [Read the Docs documentation](https://oligo-designer-toolsuite.readthedocs.io/en/latest/).



## Installation

The installation of the package can be partially done via pip, however some additional components need to be installed separately.

### Pip installation

This package was build with Python 3.10

| Package  | Version |
| ------------- | ------------- |
| argparse  | 1.4.0  |
| Bio  | 1.3.8  |
| datetime | 4.4 |
| gtfparse  | 1.2.1 |
| iteration_utilities  | 0.11.0 |
| networkx  | 2.8.1 |
| pandas  | 1.4.2 |
| pybedtools  | 0.9.0 |
| pyfaidx  | 0.6.4 |
| pyyaml  | 6.0 |
| joblib | 1.2.0 |
| bcbio-gff  | 0.6.9 |
| six  | 1.16.0 |


All packages listed above are automatically installed if the installation is done via ```pip```.


**Install Options:**

PyPI install:

```
pip install oligo-designer-toolsuite
```

Installation of the package via pip from source:

```
git clone https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite.git

pip install .        (Installation as python package: run inside directory)

pip install -e .        (Development Installation as python package: run inside directory)
```

Note: if you are using conda, first install pip with: ```conda install pip```

### Additional packages

In addition to the packages listed above, you need to install **Blast Software**, **BedTools**, **Bowtie** and **Bowtie2**.

- **Blast** can be instelled via [NCBI webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  or via ```Bioconda``` installation of Blast with:

		conda install -c bioconda blast

- **BedTools** can be installed via [BedTools GitHub](https://bedtools.readthedocs.io/en/latest/content/installation.html) or via Bioconda installation of BedTools with:

		conda install -c bioconda bedtools

- **Bowtie** and **Bowtie2** can be installed with :

		conda install -c bioconda bowtie to install Bowtie package
		conda install -c bioconda bowtie2 to install the Bowtie 2 package
