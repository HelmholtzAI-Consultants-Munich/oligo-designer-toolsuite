Installation
============

Requirements:
-------------------

This packages was tested for ``Python 3.9 - 3.10`` on ubuntu. It depends on the following additional tools **Blast**, **BedTools**, **Bowtie** and **Bowtie2** that need to be installed independently. To install those tools via conda, please activate the Bioconda and conda-forge channels in your conda environment with and update conda and all packages in your environment:

::

	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda update conda
	conda update --all


Follow this instruction to install the required additional tools:

- **Blast** (2.12 or higher) can be instelled via `NCBI webpage <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__ or via `Bioconda <http://bioconda.github.io/recipes/blast/README.html>`__ installation of Blast with:

	::

		conda install "blast>=2.12"

- **BedTools** (2.30 or higher) can be installed via `BedTools GitHub <https://bedtools.readthedocs.io/en/latest/content/installation.html>`__ or via `Bioconda <http://bioconda.github.io/recipes/bedtools/README.html>`__ installation of BedTools with:

	::

		conda install "bedtools>=2.30"

- **Bowtie** (1.3 or higher) can be installed via `Bowtie webpage <https://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie>`__ or via `Bioconda <http://bioconda.github.io/recipes/bowtie/README.html>`__ installation of Bowtie with:

	::

		conda install "bowtie>=1.3.1"

- **Bowtie2** (2.5 or higher) can be installed via `Bowtie2 webpage <https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2>`__ or via `Bioconda <http://bioconda.github.io/recipes/bowtie2/README.html>`__ installation of Bowtie2 with:

	::

		conda install "bowtie2>=2.5"

All other required packages are automatically installed if installation is done via :code:`pip`.

Install Options:
-------------------

The installation of the package is done via pip. Note: if you are using conda, first install pip with: ``conda install pip``.

PyPI install:

::

	pip install oligo-designer-toolsuite


Installation from source:

::

	git clone https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite.git


- Installation as python package (run inside directory):

	::

		pip install .


- Development Installation as python package (run inside directory):

	::

		pip install -e . [dev]
