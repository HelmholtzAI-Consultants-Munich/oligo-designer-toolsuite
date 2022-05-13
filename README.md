# Oligo Designer Toolsuite

Oligonucleotides (abbrev. oligos) are short, synthetic strands of DNA or RNA that have many application areas, ranging from research to disease diagnosis or therapeutics. Oligos can be used as primers during DNA amplification, as probes for *in situ* hybridization or as guide RNAs for CRISPR-based gene editing. Based on the intended application and experimental design, researchers can customize the length, sequence composition, and thermodynamic properties of the designed oligos.

Various tools exist that provide custom design of oligo sequences depending on the area of application. Interestingly, all those pipelines have many common basic processing steps, ranging from the generation of custom-length oligo sequences, the filtering of oligo sequences based on thermodynamic properties as well as the selection of an optimal set of oligos. Despite the fact that most tools apply the same basic processing steps, each newly developed tool usually uses its own implementation and different versions of package dependencies for those basic processing steps. As a consequence, the comparability of tools that differ only in certain steps is hampered, but also the development of new tools and the update of existing tools is slowed down, because developers do not have a common resource for basic functionalities to fall back on. We tackle this issue by providing such a common resource in our *Oligo Designer Toolsuite*. This Toolsuite is a collection of modules that provide all basic functionalities for custom oligo design pipelines within a flexible Python framework. All modules have a standardized I/O format and can be combined individually depending on the required processing steps. 

## Installation

**Requirements:**

- >= Python 3.8 
- ```'datetime```, ```argparse```, ```pyyaml```, ```iteration_utilities```, ```pandas```, 
- ```Bio```, ```gtfparse```, ```pyfaidx```,  ```pybedtools```, ```networkx```

All required packages are automatically installed if installation is done via ```pip```.

**Install Options:**

PyPI install:

```
pip install oligo-designer-toolbox
```

Installation of the package via pip from source:

Clone the git repo and install the downloaded package with:

```
git clone https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite.git
```

```
pip install .        (Installation as python package: run inside directory)
```

or if you want to develop the package:

```
pip install -e .        (Installation as python package: run inside directory)
```

Note: if you are using conda, first install pip with: ```conda install pip```


# Implemented Oligo Design Pipelines

## Padlock Probe Design

A padlock probe contains a constant backbone sequence of 53 nucleotides (nt) and the 5’- and 3’- arms, which are complementary to the corresponding mRNA sequence. The gene-specific arms of padlock probes are around 20nt long each, thus the total length of the gene-specific sequence of each padlock is 40nt.


### Usage

To create padlock probes you can run the pipeline with 

```
padlock_probe_designer -c ./config/padlock_probe_designer.yaml -o output/ [-d False]
````

where:

- ```-c```: config file, which contains parameter settings, specific to padlock probe design, default: ```./config/padlock_probe_designer.yaml``` with default parameter settings
- ```-o```: output folder, where results of pipeline are stored
  - ```annotations```folder: downloaded gene and genome annotation as well as constructed transcriptome
  - ```probes```folder: list of probes per gene, which fulfill user-defined criteria, given in config file
  - ```probesets```folder: sets of non-overlapping probes per gene, ranked by best set criteria
  - ```padlock_probes```folder: final padlock probe sequences per gene, ready to order
- ```-d```: optional, 'download only' option, where only gene and genome annotation files are downloaded but no probes generated

All steps and config parameters will be documented in a log file, that is saved in the directory where the pipeline is executed from. The logging file will have the format: ```log_padlock_probe_designer_{year}-{month}-{day}-{hour}-{minute}.txt```.
