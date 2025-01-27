<div align="center">

# *Oligo Designer Toolsuite* - Lightweight Development of Custom Oligo Design Pipelines

[![test](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test.yml/badge.svg)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test.yml)
[![PyPI](https://img.shields.io/pypi/v/oligo-designer-toolsuite.svg)](https://pypi.org/project/oligo-designer-toolsuite)
[![codecov](https://codecov.io/gh/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/branch/main/graph/badge.svg)](https://codecov.io/gh/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite)
[![stars](https://img.shields.io/github/stars/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite?logo=GitHub&color=yellow)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/stargazers)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/397343029.svg)](https://zenodo.org/badge/latestdoi/397343029)

[📃 Documentation]

[📃 Documentation]: https://oligo-designer-toolsuite.readthedocs.io/


</div>

Oligonucleotides (abbrev. oligos) are short, synthetic strands of DNA or RNA that are designed with respect to a specific target region and have many application areas,
ranging from research to disease diagnosis or therapeutics. Oligos can be used as primers during DNA amplification, as probes for in situ hybridization or as guide RNAs for CRISPR-based gene editing. Based on the intended application and experimental design, researchers have to customize the length, sequence composition, and thermodynamic properties of the designed oligos.

<div align="center">

<img src="https://raw.githubusercontent.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/main/docs/source/_figures/oligo_design.png" width="800">

</div>


Various tools exist that provide custom design of oligo sequences depending on the area of application. Even though most tools apply the same basic processing steps, ranging from the generation of custom-length oligo sequences, the filtering of oligo sequences based on thermodynamic properties as well as the selection of an optimal set of oligos, each newly developed tool uses its own implementation and different package dependencies. Consequently, not only the development of new tools is slowed down, but also the maintenance and modification of existing tools is hampered, because developers do not have a common resource for those functionalities to use. We tackle this issue with our open-source *Oligo Designer Toolsuite*.

🚀 ***Oligo Designer Toolsuite*** **is a collection of modules that provides all basic functionalities for custom oligo design pipelines as well as advanced experiment-specific functionalities like machine learning models for oligo specificity prediction within a flexible Python framework.** 

To allow the flexible usage of different modules, depending on the required processing steps, we developed a common underlying data structure that ensures the cross-compatibility of all modules within the framework. This data structure is runtime and memory optimized to enable the processing of large sequence dataset in a reasonable time frame. With our Oligo Designer Toolsuite we aim to set new standards in the development of oligo design pipelines, helping to accelerate the development of new tools and facilitate the upgrade of existing tools with the latest developments in the field. We also provide ready-to-use oligo design pipelines for specific experimental setups, e.g. SCRINSHOT or SeqFISH+ probe design for Spatial Transcriptomics.

## Installation

**Requirements:**

This packages was tested for ```Python 3.9 - 3.10``` on ubuntu and macos. For stable installation, we recommend to first setup a conda environment, e.g.:

```
conda create -n odt python=3.10
conda activate odt
```

*Note: if your institution does not support anaconda, you can use [miniforge](https://github.com/conda-forge/miniforge) instead to run the conda installations.*

If you have an Apple M chip, you need to create an environment simulating an x86 processor to be able to install **Blast**. This can be done as follows:

```
CONDA_SUBDIR=osx-64 conda create -n odt python=3.10
conda activate odt
conda config --env --set subdir osx-64
```

It depends on the following additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** that need to be installed independently. To install those tools via conda, please activate the Bioconda and conda-forge channels in your conda environment with and update conda and all packages in your environment:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all
```

Follow this instruction to install the required additional tools:

- **Blast** (2.15 or higher) can be installed via [NCBI webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) or via [Bioconda](http://bioconda.github.io/recipes/blast/README.html) installation of Blast with:

		conda install "blast>=2.15.0"

- **BedTools** (2.30 or higher) can be installed via [BedTools GitHub](https://bedtools.readthedocs.io/en/latest/content/installation.html) or via [Bioconda](http://bioconda.github.io/recipes/bedtools/README.html) installation of BedTools with:

		conda install "bedtools>=2.30"

- **Bowtie** (1.3 or higher) can be installed via [Bowtie webpage](https://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie) or via [Bioconda](http://bioconda.github.io/recipes/bowtie/README.html) installation of Bowtie with:

		conda install "bowtie>=1.3.1"

- **Bowtie2** (2.5 or higher) can be installed via [Bowtie2 webpage](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2) or via [Bioconda](http://bioconda.github.io/recipes/bowtie2/README.html) installation of Bowtie2 with:

		conda install "bowtie2>=2.5"

All other required packages are automatically installed if installation is done via ```pip```.

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
```

- Installation as python package (run inside directory):

		pip install .


- Development installation as python package (run inside directory):

		pip install -e .[dev]


## Implemented Oligo Design Pipelines

The following pipelines are pre-implemented and ready-to-use:

[🧬 Genomic Region Generator](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/genomic_region_generator.html)

[🧪 Scrinshot Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/scrinshot_probe_designer.html)

[🧪 SeqFISH+ Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/seqfishplus_probe_designer.html)

[🧪 MERFISH Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/merfish_probe_designer.html)

[🧫 Oligo-Seq Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/oligoseq_probe_designer.html)

If you would like to modify an existing pipeline for your specific experimental settings or would like to design a new pipeline for a different experimental setup, feel free to reach out to us [Lisa Barros de Andrade e Sousa](mailto:lisa.barros@helmholtz-munich.de) or [Isra Mekki](mailto:isra.mekki@helmholtz-munich.de).

## Contributing

Contributions are more than welcome! Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you can't contribute. To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.

For any further inquiries please send an email to [Lisa Barros de Andrade e Sousa](mailto:lisa.barros@helmholtz-munich.de) or [Isra Mekki](mailto:isra.mekki@helmholtz-munich.de).


## How to cite

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

If you are using one of the spatial transcriptomics pipelines provided along the Oligo Designer Toolsuite, consider citing in addition the paper:

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

```oligo-designer-toolsuite``` is released under the MIT license. See [LICENSE](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/dev/LICENSE) for additional details about it.
