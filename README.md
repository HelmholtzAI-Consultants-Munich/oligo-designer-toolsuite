<div align="center">

# *Oligo Designer Toolsuite* - Lightweight Development of Custom Oligo Design Pipelines

[![Docs](https://img.shields.io/badge/docs-latest-blue?style=flat&logo=readthedocs)](https://oligo-designer-toolsuite.readthedocs.io/en/latest/)
[![PyPI](https://img.shields.io/pypi/v/oligo-designer-toolsuite.svg)](https://pypi.org/project/oligo-designer-toolsuite)
[![PyPI Downloads](https://static.pepy.tech/badge/oligo-designer-toolsuite)](https://pepy.tech/projects/oligo-designer-toolsuite)
[![DOI](https://zenodo.org/badge/397343029.svg)](https://zenodo.org/badge/latestdoi/397343029)
[![stars](https://img.shields.io/github/stars/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite?logo=GitHub&color=yellow)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/stargazers)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![codecov](https://codecov.io/gh/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/branch/main/graph/badge.svg)](https://codecov.io/gh/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite)
[![TestUbuntuX64](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test_ubuntu_x64.yml/badge.svg)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test_ubuntu_x64.yml)
[![TestMacOsArm64](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test_macos_arm64.yml/badge.svg)](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/actions/workflows/test_macos_arm64.yml)

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

[🧫 CycleHCR Probe Designer](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/cyclehcr_probe_designer.html)

If you would like to modify an existing pipeline for adjusted experimental settings or to design a new pipeline for a different experimental setup, feel free to reach out to us [Lisa Barros de Andrade e Sousa](mailto:lisa.barros@helmholtz-munich.de) or [Jonas Hagenberg](mailto:jonas.hagenberg@helmholtz-munich.de).

## Installation

<!-- LINK INSTALLATION START -->

### Requirements

This package was tested for ***Python 3.10 - 3.12*** on ***Linux (x64)*** and ***macOS (arm64)*** but not tested on Windows.

*Note: Intel-based macOS (osx64) installation is ⚠️ DEPRECATED ⚠️*

For stable installation, we recommend to first setup a conda environment.

*Note: if your institution does not support anaconda, you can use [miniforge](https://github.com/conda-forge/miniforge) instead to run the conda installations.*

First create a conda environment:

```
conda create -n odt python=3.12
conda activate odt
```

To install the additional required tools via conda, please activate the *bioconda* and *conda-forge* channels in your conda environment and update conda and all packages in your environment:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all
```

The additional tools need to be installed independently:

```
conda install "blast>=2.15.0"
conda install "bedtools>=2.30"
conda install "bowtie>=1.3.1"
conda install "bowtie2>=2.5"
conda install "bcftools>=1.22"
conda install "samtools>=1.22"
```

All other required packages are automatically installed if installation is done via ```pip``` (see below).

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

For any further inquiries please send an email to [Lisa Barros de Andrade e Sousa](mailto:lisa.barros@helmholtz-munich.de) or [Jonas Hagenberg](mailto:jonas.hagenberg@helmholtz-munich.de).

<!-- LINK CONTRIBUTION END -->

## How to cite

<!-- LINK CITE START -->

If the Oligo Designer Toolsuite is useful for your research, consider citing the package:

```
@software{
	author		= 	{	Barros de Andrade e Sousa L.,
  						Mekki I.,
						Campi F.,
						Kümmerle L.,
						Bright C.,
						Lücken M.,
						Theis F.,
						Piraud M.
					},
	title		= 	{Oligo Designer Toolsuite},
	year		= 	{2025},
	publisher	= 	{GitHub},
	journal 	= 	{GitHub repository},
	url 		= 	{https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite}
}
```

<!-- LINK CITE END -->

If you are using the SCRINSHOT, MERFISH or SeqFISH+ pipeline provided along the Oligo Designer Toolsuite, consider citing in addition the paper:

```
@article{
	author 	 	= 	{
						Louis B. Kuemmerle,
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
		     			Fabian J. Theis
					},
    title 	 	= 	{Probe set selection for targeted spatial transcriptomics},
    year 	 	= 	{2024},
    publisher 	= 	{Nature Publishing Group US New York},
    journal 	= 	{Nature methods},
    doi 	 	= 	{10.1038/s41592-024-02496-z},
    URL 	 	= 	{https://doi.org/10.1038/s41592-024-02496-z}
}
```

## License

<!-- LINK LICENSE START -->

```oligo-designer-toolsuite``` is released under the MIT license. See [LICENSE](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/dev/LICENSE) for additional details about it.

<!-- LINK LICENSE END -->
