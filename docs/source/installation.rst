Installation
============

**Requirements:**

This package was build with Python 3.10

=================== =======
Package             Version
=================== =======
argparse            1.4.0
Bio                 1.3.8
datetime            4.4
gtfparse            1.2.1
iteration_utilities 0.11.0
networkx            2.8.1
pandas              1.4.2
pybedtools          0.9.0
pyfaidx             0.6.4
pyyaml              6.0
joblib		    1.2.0
bcbio-gff 	    0.6.9
six                 1.16.0
=================== =======

All required packages are automatically installed if installation is
done via ``pip``.

**Install Options:**

PyPI install:

::

   pip install oligo-designer-toolbox

Installation of the package via pip from source:

::

   git clone https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite.git

   pip install .        (Installation as python package: run inside directory)

   pip install -e .        (Development Installation as python package: run inside directory)

Note: if you are using conda, first install pip with:
``conda install pip``

In addition to the packages listed above, you need to install **Blast
Software**, **BedTools**, **Bowtie** and **Bowtie2**.

**Blast** can be instelled via `NCBI
webpage <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__
or via ``Bioconda`` installation of Blast with
``conda install -c bioconda blast``.

**BedTools** can be installed via `BedTools GitHub <https://bedtools.readthedocs.io/en/latest/content/installation.html>`__
or via Bioconda installation of BedTools with ``conda install -c bioconda bedtools``.

**Bowtie** and **Bowtie2** can be installed with ``conda install -c bioconda bowtie to install Bowtie package`` and ``conda install -c bioconda bowtie2 to install the Bowtie 2 package``.
