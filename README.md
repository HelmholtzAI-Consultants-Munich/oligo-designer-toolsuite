<div align="center">

# *Oligo Designer Toolsuite* - Lightweight Development of Custom Oligo Design Pipelines

[![test](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test.yml/badge.svg)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test.yml)
[![PyPI](https://img.shields.io/pypi/v/oligo-designer-toolsuite.svg)](https://pypi.org/project/oligo-designer-toolsuite)
[![codecov](https://codecov.io/gh/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/branch/pipelines/graph/badge.svg)](https://codecov.io/gh/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite)
[![stars](https://img.shields.io/github/stars/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite?logo=GitHub&color=yellow)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/stargazers)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/397343029.svg)](https://zenodo.org/badge/latestdoi/397343029)

[Docs] | [Tutorials]

[Docs]: https://oligo-designer-toolsuite.readthedocs.io/
[Tutorials]: https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/tree/dev/tutorials

</div>

Oligonucleotides (abbrev. oligos) are short, synthetic strands of DNA or RNA that are designed with respect to a specific target region and have many application areas,
ranging from research to disease diagnosis or therapeutics. Oligos can be used as primers during DNA amplification, as probes for in situ hybridization or as guide RNAs for CRISPR-based gene editing.
Based on the intended application and experimental design, researchers have to customize the length, sequence composition, and thermodynamic properties of the designed oligos.

<div align="center">

<img src="https://raw.githubusercontent.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/dev/docs/source/_figures/oligo_design.png" width="800">

</div>


Various tools exist that provide custom design of oligo sequences depending on the area of application. Interestingly, all those pipelines have many common basic processing steps,
ranging from the generation of custom-length oligo sequences, the filtering of oligo sequences based on thermodynamic properties as well as the selection of an optimal set of oligos.
Despite the fact that most tools apply the same basic processing steps, each newly developed tool usually uses its own implementation and different versions of package dependencies for those basic processing steps.
As a consequence, the comparability of tools that differ only in certain steps is hampered, but also the maintenance of existing tools and the development of new tools is slowed down,
because developers do not have a common resource for basic functionalities to use. We tackle this issue by providing such a common resource in our open-source *Oligo Designer Toolsuite*.

***Oligo Designer Toolsuite*** **is a collection of modules that provide all basic functionalities for custom oligo design pipelines within a flexible Python framework.**
Furthermore, we introduce a common underlying data structure, which allows the user to easily combine different modules, depending on the required processing steps.
We also provide ready-to-use oligo design pipelines for specific experimental setups, e.g. SCRINSHOT or SeqFISH+ probe design for Spatial Transcriptomics.



## Installation

**Requirements:**

This packages was tested for ```Python 3.9 - 3.10``` on ubuntu and macos. For stable installatio, we recommend to first setup a conda environment, e.g.:

```
conda create -n odt python=3.10
conda activate odt
```

It depends on the following additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** that need to be installed independently. To install those tools via conda, please activate the Bioconda and conda-forge channels in your conda environment with and update conda and all packages in your environment:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all
```

Follow this instruction to install the required additional tools:

- **Blast** (2.15 or higher) can be instelled via [NCBI webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) or via [Bioconda](http://bioconda.github.io/recipes/blast/README.html) installation of Blast with:

		conda install "blast>=2.15.0"

- **BedTools** (2.30 or higher) can be installed via [BedTools GitHub](https://bedtools.readthedocs.io/en/latest/content/installation.html) or via [Bioconda](http://bioconda.github.io/recipes/bedtools/README.html) installation of BedTools with:

		conda install "bedtools>=2.30"

- **Bowtie** (1.3 or higher) can be installed via [Bowtie webpage](https://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie) or via [Bioconda](http://bioconda.github.io/recipes/bowtie/README.html) installation of Bowtie with:

		conda install "bowtie>=1.3.1"

- **Bowtie2** (2.5 or higher) can be installed via [Bowtie2 webpage](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2) or via [Bioconda](http://bioconda.github.io/recipes/bowtie2/README.html) installation of Bowtie2 with:

		conda install "bowtie2>=2.5"

All other required packages are automatically installed if installation is done via :code:`pip`.

**Install Options:**

The installation of the package is done via pip. Note: if you are using conda, first install pip with: ```conda install pip```.

PyPI install:

```
pip install oligo-designer-toolsuite
```


Installation from source:

```
git clone https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite.git
cd oligo-designer-toolsuite
git switch pipelines
```

- Installation as python package (run inside directory):

		pip install .


- Development Installation as python package (run inside directory):

		pip install -e .[dev]


## Implemented Oligo Design Pipelines

### Scrinshot Probe Design

A padlock probe contains a constant backbone sequence of 53 nucleotides (nt) and the 5’- and 3’- arms, which are complementary to the corresponding mRNA sequence. 
The gene-specific arms of padlock probes are around 20nt long each, thus the total length of the gene-specific sequence of each padlock is around 40nt.


#### Usage

*Command-Line Call:*

To create scrinshot probes you can run the pipeline with 

```
scrinshot_probe_designer -c data/configs/scrinshot_probe_designer.yaml
````

where:

- ```-c```: config file, which contains parameter settings, specific to scrinshot probe design, *scrinshot_probe_designer.yaml* contains default parameter settings

All steps and config parameters will be documented in a log file, that is saved in the directory where the pipeline is executed from. 
The logging file will have the format: ```log_scrinshot_probe_designer_{year}-{month}-{day}-{hour}-{minute}.txt```.

### Oligo-Seq Probe Design

An oligo-seq probe is an oligo hybridization probe, which is optimized for probe-based targeted sequencing to measure RNA expression.

#### Usage

*Command-Line Call:*

To create oligo-seq probes you can run the pipeline with 

```
oligo_seq_probe_designer -c data/configs/oligo_seq_probe_designer.yaml
````

where:

- ```-c```: config file, which contains parameter settings, specific to oligo-seq probe design, *oligo_seq_probe_designer.yaml* contains default parameter settings

All steps and config parameters will be documented in a log file, that is saved in the directory where the pipeline is executed from. 
The logging file will have the format: ```log_oligo_seq_probe_designer_{year}-{month}-{day}-{hour}-{minute}.txt```.


## Contributing

Contributions are more than welcome! Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you can't contribute. To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.

## How to cite

If the Oligo Designer Toolsuite is useful for your research, consider citing the package:

```
@software{campi_2023_7823048,
    author       = { Isra Mekki,
		     Francesco Campi,  
                     Louis Kümmerle,
		     Chelsea Bright,
		     Malte Lücken
                     Fabian Theis,
                     Marie Piraud,
                     Lisa Barros de Andrade e Sousa
                     },
    title        = {{Oligo Designer Toolsuite}},
    year         = 2023,
    publisher    = {Zenodo},
    version      = {v0.1.3},
    doi          = {10.5281/zenodo.7823048},
    url          = {https://doi.org/10.5281/zenodo.7823048}
}
```

If you are using one of the spatial transcriptomics pipelines provided along the Oligo Designer Toolsuite, consider citing in addition the paper:

```
@article {Kuemmerle2022.08.16.504115,
    author 	 = { Louis B. Kuemmerle,
		     Malte D. Luecken,
		     Alexandra B. Firsova
		     Lisa Barros de Andrade e Sousa
		     Lena Strasser
		     Lukas Heumos
		     Ilhem Isra Mekki
		     Krishnaa T. Mahbubani
		     Alexandros Sountoulidis
		     Tamas Balassa
		     Ferenc Kovacs
		     Peter Horvath
		     Marie Piraud
		     Ali Ertürk
		     Christos Samakovlis
		     Fabian J. Theis},
    title 	 = {{Probe set selection for targeted spatial transcriptomics}},
    year 	 = {2022},
    publisher 	 = {Cold Spring Harbor Laboratory},
    journal 	 = {bioRxiv},
    doi 	 = {10.1101/2022.08.16.504115},
    URL 	 = {https://www.biorxiv.org/content/early/2022/08/17/2022.08.16.504115}
}
```

## License

```oligo-designer-toolsuite``` is released under the MIT license. See [LICENSE](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/dev/LICENSE) for additional details about it.
