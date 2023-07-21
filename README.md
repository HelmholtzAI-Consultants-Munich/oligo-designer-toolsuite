<div align="center">

# *Oligo Designer Toolsuite* - Lightweight Development of Custom Oligo Design Pipelines

[![test](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test.yml/badge.svg)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test.yml)
[![PyPI](https://img.shields.io/pypi/v/oligo-designer-toolsuite.svg)](https://pypi.org/project/oligo-designer-toolsuite)
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

This packages was tested for ```Python 3.9 - 3.10``` on ubuntu. It depends on the following additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** that need to be installed independently. To install those tools via conda, please activate the Bioconda and conda-forge channels in your conda environment with and update conda and all packages in your environment:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda update conda
conda update --all
```

Follow this instruction to install the required additional tools:

- **Blast** (2.12 or higher) can be instelled via [NCBI webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) or via [Bioconda](http://bioconda.github.io/recipes/blast/README.html) installation of Blast with:

		conda install "blast>=2.12"

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
```

- Installation as python package (run inside directory):

		pip install .   


- Development Installation as python package (run inside directory):

		pip install -e .[dev]



## Contributing

Contributions are more than welcome! Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you can't contribute. To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.

## How to cite

If the Oligo Designer Toolsuite is useful for your research, consider citing the package:

```
@software{lisa_sousa_2023_7823048,
    author       = { Isra Mekki, 
		     Francesco Campi, 
		     Louis Kümmerle,
		     Hanane Mohaouchane, 
		     Maksym Tretiakov, 
		     Anna Starovoit,
		     Cheng-Wei Liao,
			 Marie Piraud,
		     Lisa Barros de Andrade e Sousa},
    title        = {{Oligo Designer Toolsuite}},
    month        = april,
    year         = 2023,
    publisher    = {Zenodo},
    version      = {v0.1.3},
    doi          = {10.5281/zenodo.7823048},
    url          = {https://doi.org/10.5281/zenodo.7823048}
}
```

## License

```oligo-designer-toolsuite``` is released under the MIT license. See [LICENSE](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/dev/LICENSE) for additional details about it.
