<div align="center">

# *Oligo Designer Toolsuite* - Lightweight Development of Custom Oligo Design Pipelines

[![Documentation Status](https://readthedocs.org/projects/oligo-designer-toolsuite/badge/?version=latest)](https://oligo-designer-toolsuite.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/oligo-designer-toolsuite.svg)](https://pypi.org/project/oligo-designer-toolsuite)
[![DOI](https://zenodo.org/badge/397343029.svg)](https://zenodo.org/badge/latestdoi/397343029)
[![stars](https://img.shields.io/github/stars/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite?logo=GitHub&color=yellow)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/stargazers)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![TestUbuntuX64](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test_ubuntu_x64.yml/badge.svg)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test_ubuntu_x64.yml)
[![TestMacOsArm64](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test_macos_arm64.yml/badge.svg)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test_macos_arm64.yml)
[![codecov](https://codecov.io/gh/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/branch/main/graph/badge.svg)](https://codecov.io/gh/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite)

<!-- LINK INTRODUCTION START -->

</div>

Oligonucleotides (abbrev. oligos) are short, synthetic strands of DNA or RNA that are designed with respect to a specific target region and have many application areas,
ranging from research to disease diagnosis or therapeutics. Oligos can be used as primers during DNA amplification, as probes for in situ hybridization or as guide RNAs for CRISPR-based gene editing. Based on the intended application and experimental design, researchers have to customize the length, sequence composition, and thermodynamic properties of the designed oligos.

<div align="center">

<img src="https://raw.githubusercontent.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/main/docs/source/_static/oligo_design.png" width="800">

</div>


Various tools exist that provide custom design of oligo sequences depending on the area of application. Even though most tools apply the same basic processing steps, ranging from the generation of custom-length oligo sequences, the filtering of oligo sequences based on thermodynamic properties as well as the selection of an optimal set of oligos, each newly developed tool uses its own implementation and different package dependencies. Consequently, not only the development of new tools is slowed down, but also the maintenance and modification of existing tools is hampered, because developers do not have a common resource for those functionalities to use. We tackle this issue with our open-source *Oligo Designer Toolsuite*.

🚀 ***Oligo Designer Toolsuite*** **is a collection of modules that provides all basic functionalities for custom oligo design pipelines as well as advanced experiment-specific functionalities like machine learning models for oligo specificity prediction within a flexible Python framework.** 

To allow the flexible usage of different modules, depending on the required processing steps, we developed a common underlying data structure that ensures the cross-compatibility of all modules within the framework. This data structure is runtime and memory optimized to enable the processing of large sequence dataset in a reasonable time frame. With our Oligo Designer Toolsuite we aim to set new standards in the development of oligo design pipelines, helping to accelerate the development of new tools and facilitate the upgrade of existing tools with the latest developments in the field. We also provide ready-to-use oligo design pipelines for specific experimental setups, e.g. MERFISH or SeqFISH+ probe design for Spatial Transcriptomics.

<!-- LINK INTRODUCTION END -->

## Implemented Oligo Design Pipelines

The following pipelines are pre-implemented and ready-to-use:

[🧬 Genomic Region Generator](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/genomic_region_generator.html)

[🧪 Scrinshot Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/scrinshot_probe_designer.html)

[🧪 SeqFISH+ Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/seqfishplus_probe_designer.html)

[🧪 MERFISH Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/merfish_probe_designer.html)

[🧫 Oligo-Seq Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/oligoseq_probe_designer.html)

If you would like to modify an existing pipeline for adjusted experimental settings or to design a new pipeline for a different experimental setup, feel free to reach out to us [Lisa Barros de Andrade e Sousa](mailto:lisa.barros@helmholtz-munich.de) or [Isra Mekki](mailto:isra.mekki@helmholtz-munich.de).

## Installation

<!-- LINK INSTALLATION START -->

### Requirements

This packages was tested for ***Python 3.9 - 3.12*** on ***Linux (x64)*** and ***MacOS (osx64 and arm64)***.

For stable installation, we recommend to first setup a conda environment.

*Note: if your institution does not support anaconda, you can use [miniforge](https://github.com/conda-forge/miniforge) instead to run the conda installations.*

**🖥️ Linux Requirements**

First create a conda environment:

```
conda create -n odt python=3.11
conda activate odt
```

To install the additional required tools via conda, please activate the *bioconda* and *conda-forge* channels in your conda environment and update conda and all packages in your environment:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all
```

The additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** need to be installed independently:

```
conda install "blast>=2.15.0" 
conda install "bedtools>=2.30"
conda install "bowtie>=1.3.1"
conda install "bowtie2>=2.5"
```

All other required packages are automatically installed if installation is done via ```pip``` (see below).


**🖥️ MacOS M Chip (arm64) Requirements**

For the Apple M chips, there is currently no Blast installation available via conda. Hence, we need a workaround for the Blast installation. Therefore, we have two options:

***Option 1: Blast installation via Homebrew***

*Pro: Allows to install a native M Chip (arm64) version of Blast and supports installation of* ```torch > 2.2.2```.   
*Con: Requires installation Homebrew, which needs sudo rights.*

First create a conda environment:

```
conda create -n odt-arm64 python=3.11
conda activate odt-arm64
```

To install the additional required tools via conda, please activate the *bioconda* and *conda-forge* channels in your conda environment and update conda and all packages in your environment:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all
```

The additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** need to be installed independently. **BedTools**, **Bowtie** and **Bowtie2** can be installed via conda:

```
conda install "bedtools>=2.30"
conda install "bowtie>=1.3.1"
conda install "bowtie2>=2.5"
```

To install the M Chip (arm64) version of Blast you need Homebrew, which can be installed as described [here](https://brew.sh/). Blast can then be installed via Homebrew:

```
brew install blast
```

All other required packages are automatically installed if installation is done via ```pip``` (see below).

***Option 2: Blast installation via conda installation through Intel Chip (osx64) environment emulation***

*Pro: Works with conda and requires no extra dependencies.*  
*Con: Runs via Rosetta osx64 emulation, which does not support installation of* ```torch > 2.2.2```

First we need to create an conda environment that emulates the osx64 processor:

```
CONDA_SUBDIR=osx-64 conda create -n odt-osx64 python=3.11
conda activate odt-osx64
conda config --env --set subdir osx-64
```

To install the additional required tools via conda, please activate the *bioconda* and *conda-forge* channels in your conda environment and update conda and all packages in your environment:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all
```

The following additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** need to be installed independently:

```
conda install "blast>=2.15.0" 
conda install "bedtools>=2.30"
conda install "bowtie>=1.3.1"
conda install "bowtie2>=2.5"
```

Since ```torch > 2.2.2``` installation is not provided anymore for Max Intel Chips (osx64 processor), we need to make sure to have ```numpy < 2.0``` to avoid conflicts which ```torch <= 2.2.2```:

```
pip install "numpy<2.0"
```

All other required packages are automatically installed if installation is done via ```pip``` (see below).


**🖥️ MacOS Intel Chip (osx64) Requirements**

First create a conda environment:

```
conda create -n odt-x64 python=3.11
conda activate odt-x64
```

To install the additional required tools via conda, please activate the *bioconda* and *conda-forge* channels in your conda environment and update conda and all packages in your environment:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all
```

The following additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** need to be installed independently:

```
conda install "blast>=2.15.0" 
conda install "bedtools>=2.30"
conda install "bowtie>=1.3.1"
conda install "bowtie2>=2.5"
```

Since ```torch > 2.2.2``` installation is not provided anymore for Max Intel Chips (osx64 processor), we need to make sure to have ```numpy < 2.0``` to avoid conflicts which ```torch <= 2.2.2```:

```
pip install "numpy<2.0"
```

All other required packages are automatically installed if installation is done via ```pip``` (see below).

**🖥️ Windows Requirements**

TBD

The following additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** need to be installed independently:

- **Blast** (2.15 or higher) can be installed via [NCBI webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

- **BedTools** (2.30 or higher) can be installed via [BedTools GitHub](https://bedtools.readthedocs.io/en/latest/content/installation.html)

- **Bowtie** (1.3 or higher) can be installed via [Bowtie webpage](https://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie)

- **Bowtie2** (2.5 or higher) can be installed via [Bowtie2 webpage](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2) 


### Install Options

The installation of the package is done via pip. Note: if you are using conda, first install pip with: ```conda install pip```.

PyPI install:

```
pip install oligo-designer-toolsuite
```


Installation from source:

```
git clone https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite.git
cd oligo-designer-toolsuite
```

- Installation as python package (run inside directory):

		pip install .


- Development installation as python package (run inside directory):

		pip install -e .[dev]


<!-- LINK INSTALLATION END -->

## Contributing

<!-- LINK CONTRIBUTION START -->

Contributions are more than welcome! Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you can't contribute. To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.

For any further inquiries please send an email to [Lisa Barros de Andrade e Sousa](mailto:lisa.barros@helmholtz-munich.de) or [Isra Mekki](mailto:isra.mekki@helmholtz-munich.de).

<!-- LINK CONTRIBUTION END -->

## How to cite

<!-- LINK CITE START -->

If the Oligo Designer Toolsuite is useful for your research, consider citing the package:

```
@software{campi_2023_7823048,
	author   = {Isra Mekki,
		     Francesco Campi,  
		     Louis Kümmerle,
		     Chelsea Bright,
		     Malte Lücken
		     Fabian Theis,
		     Marie Piraud,
		     Lisa Barros de Andrade e Sousa},
    title        = {{Oligo Designer Toolsuite}},
    year         = 2023,
    publisher    = {Zenodo},
    version      = {v0.1.3},
    doi          = {10.5281/zenodo.7823048},
    url          = {https://doi.org/10.5281/zenodo.7823048}
}
```

<!-- LINK CITE END -->

If you are using the SCRINSHOT, MERFISH or SeqFISH+ pipeline provided along the Oligo Designer Toolsuite, consider citing in addition the paper:

```
@article {kuemmerle2024probe,
    author 	 = { Louis B. Kuemmerle,
		     Malte D. Luecken,
		     Alexandra B. Firsova
		     Lisa Barros de Andrade e Sousa
		     Lena Strasser
                     Ilhem Isra Mekki
                     Francesco Campi
		     Lukas Heumos
		     Maiia Shulman
                     Valentina Beliaeva
                     Soroor Hediyeh-Zadeh
                     Anna C. Schaar
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
    year 	 = {2024},
    publisher 	 = {Nature Publishing Group US New York},
    journal 	 = {Nature methods},
    doi 	 = {10.1038/s41592-024-02496-z},
    URL 	 = {https://doi.org/10.1038/s41592-024-02496-z}
}
```

## License

<!-- LINK LICENSE START -->

```oligo-designer-toolsuite``` is released under the MIT license. See [LICENSE](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/dev/LICENSE) for additional details about it.

<!-- LINK LICENSE END -->
