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

To allow the flexible usage of different modules, depending on the required processing steps, we developed a common underlying data structure that ensures the cross-compatibility of all modules within the framework. This data structure is runtime and memory optimized to enable the processing of large sequence dataset in a reasonable time frame. With our Oligo Designer Toolsuite we aim to set new standards in the development of oligo design pipelines, helping to accelerate the development of new tools and facilitate the upgrade of existing tools with the latest developments in the field. *We also provide ready-to-use oligo design pipelines for multiple experimental applications*.

## ODT-Cloud

Use Oligo Designer Toolsuite directly in your browser through **ODT-Cloud**:

🌐 https://odt.helmholtz-munich.de/

ODT-Cloud provides a graphical user interface for running oligo design pipelines without local installation or dependency management. For advanced customization, genome-wide applications, and integration into automated analysis pipelines, we recommend using the command-line interface.


<!-- LINK INTRODUCTION END -->

## Installation

<!-- LINK INSTALLATION START -->

### Requirements

This package is continuously tested with **Python 3.10 – 3.12** on **Linux (x86_64)** and **macOS (Apple Silicon / arm64)**.
Windows is currently not supported.

> **Note:** Intel-based macOS (x86_64) is deprecated and no longer actively tested.
>
> *How to verify that you are using Apple Silicon?*
> ```bash
> python -c "import platform; print(platform.machine())"
> ```
> Expected output: `arm64`

For a stable installation, we recommend using a dedicated Conda environment via [miniforge](https://github.com/conda-forge/miniforge). Create a new Conda environment using the recommended channels and activate the environment:
```
conda create -n odt "python=3.12" "mamba" "pip"
conda activate odt
```

Install the required third-party tools:
```
mamba install -y \
  -c conda-forge \
  -c bioconda \
  --override-channels \
  "blast>=2.15.0" \
  "bedtools>=2.30" \
  "bowtie>=1.3.1" \
  "bowtie2>=2.5" \
  "samtools>=1.22" \
  "bcftools>=1.22"
```

All other required packages are automatically installed if installation is done via ```pip``` (see below).

### Install Options

Install from PyPI:
```
pip install oligo-designer-toolsuite
```

Install from Source:
```
git clone https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite.git
cd oligo-designer-toolsuite
pip install .
```

<!-- LINK INSTALLATION END -->

## Quick Start

<!-- LINK QUICKSTART START -->

ODT pipelines are configured through YAML files. Example configuration files are available in: `data/configs/`.

The probe design pipelines require genomic input sequences. These can be generated with the **Genomic Region Generator**. Start from one of the example configuration files. Adjust the configuration to your species, genome annotation and output directory, then run:

```bash
genomic_region_generator --config data/configs/genomic_region_generator_ncbi.yaml
```

This generates the genomic region files that are used as input for downstream probe design pipelines.

After generating the required genomic input files, choose the configuration file for your pipeline, e.g.: `data/configs/oligo_seq_probe_designer.yaml`. Adjust the input paths and design parameters, then run:

```bash
oligo_seq_probe_designer --config data/configs/oligo_seq_probe_designer.yaml
```

The pipeline writes the designed probes and summary files to the configured output directory.

Pipelines can also be called directly from Python. For Oligo-Seq:

```python
import yaml
from oligo_designer_toolsuite.config.pipelines.oligo_seq_probe_designer import OligoSeqProbeDesignerConfig
from oligo_designer_toolsuite.pipelines import oligo_seq_probe_designer

with open("data/configs/oligo_seq_probe_designer.yaml", "r") as handle:
    config_raw = yaml.safe_load(handle)

config = OligoSeqProbeDesignerConfig.model_validate(config_raw)
oligo_seq_probe_designer(config)
```

For a complete description of all configuration parameters, see the documentation of the respective pipeline.

<!-- LINK QUICKSTART END -->

## Implemented Oligo Design Pipelines

The following pipelines are pre-implemented, documented, and ready-to-use:

**Genomic Region Generator**

Extract genomic target regions from genome annotations and reference genomes. The pipeline can automatically download genome assemblies and annotations from Ensembl or NCBI, or work with user-provided FASTA and GTF files.

📖 [Documentation](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/genomic_region_generator.html)

**SCRINSHOT Probe Designer**

Design padlock probes for SCRINSHOT (Single-Cell RNA In-Situ Hybridization and Sequencing On Tissue) experiments. The pipeline generates gene-specific padlock probes for highly multiplexed and spatially resolved RNA detection at single-cell resolution.

📖 [Documentation](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/scrinshot_probe_designer.html)

**MERFISH Probe Designer**

Design encoding probes for MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization) experiments. Generated probes contain gene-specific targeting sequences and barcode regions enabling highly multiplexed spatial transcriptomics measurements.

📖 [Documentation](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/merfish_probe_designer.html)

**SeqFISH+ Probe Designer**

Design encoding probes for SeqFISH+ (Sequential Fluorescence In Situ Hybridization) experiments. The pipeline generates barcoded probe sets for large-scale spatially resolved transcriptomics with single-cell resolution.

📖 [Documentation](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/seqfishplus_probe_designer.html)

**CycleHCR Probe Designer**

Design barcoded probe sets for CycleHCR experiments. CycleHCR combines Hybridization Chain Reaction (HCR) with iterative barcoding strategies to enable highly multiplexed RNA imaging while minimizing molecular crowding effects.

📖 [Documentation](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/cyclehcr_probe_designer.html)

**Oligo-Seq Probe Designer**

Design probe libraries for Oligo-seq, a targeted RNA profiling technology that combines in situ hybridization with sequencing-based readout. The pipeline generates probes targeting exonic regions and exon-intron junctions for sensitive and reproducible transcript detection from low-input samples.

📖 [Documentation](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/oligoseq_probe_designer.html)

For tutorials, configuration examples, API documentation, and detailed pipeline descriptions, visit: https://oligo-designer-toolsuite.readthedocs.io


## Contributing

<!-- LINK CONTRIBUTION START -->

Contributions are more than welcome. Everything from code to notebooks, examples, and documentation is equally valuable.

To contribute:

1. Fork the repository
2. Create a feature branch
3. Implement your changes
4. Submit a pull request

We will do our best to review contributions and help integrate them into the project.

For questions or collaboration inquiries, please contact:

* Lisa Barros de Andrade e Sousa ([lisa.barros@helmholtz-munich.de](mailto:lisa.barros@helmholtz-munich.de))
* Jonas Hagenberg ([jonas.hagenberg@helmholtz-munich.de](mailto:jonas.hagenberg@helmholtz-munich.de))


<!-- LINK CONTRIBUTION END -->

## How to cite

<!-- LINK CITE START -->

If Oligo Designer Toolsuite is useful for your research, please consider citing the software:

```bibtex
@software{
    author      = {
                    Barros de Andrade e Sousa L.,
                    Mekki I.,
                    Campi F.,
                    Kümmerle L.,
                    Bright C.,
                    Lücken M.,
                    Theis F.,
                    Piraud M.
                  },
    title       = {Oligo Designer Toolsuite},
    year        = {2025},
    publisher   = {GitHub},
    journal     = {GitHub repository},
    url         = {https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite}
}
```

<!-- LINK CITE END -->

## License

<!-- LINK LICENSE START -->

```oligo-designer-toolsuite``` is released under the MIT license. See [LICENSE](https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/dev/LICENSE) for additional details about it.

<!-- LINK LICENSE END -->
