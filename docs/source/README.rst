Installation
============

**Requirements:**

This package was build with Python 3.8

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

In addition to the packages listed above, you need to install *Blast
Software*. This can be done via `NCBI
webpage <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__
or via ``Bioconda`` installation of Blast with
``conda install -c bioconda blast``.

Implemented Oligo Design Pipelines
==================================

Padlock Probe Design
--------------------

A padlock probe contains a constant backbone sequence of 53 nucleotides
(nt) and the 5’- and 3’- arms, which are complementary to the
corresponding mRNA sequence. The gene-specific arms of padlock probes
are around 20nt long each, thus the total length of the gene-specific
sequence of each padlock is 40nt.

Usage
~~~~~

**Command-Line Call:**

To create padlock probes you can run the pipeline with

::

   padlock_probe_designer -c ./config/padlock_probe_designer.yaml -o output/ [-d False]

where:

-  ``-c``: config file, which contains parameter settings, specific to
   padlock probe design, *./config/padlock_probe_designer.yaml* contains
   default parameter settings
-  ``-o``: output folder, where results of pipeline are stored

   -  ``annotations``\ folder: downloaded gene and genome annotation as
      well as constructed transcriptome
   -  ``probes``\ folder: list of probes per gene, which fulfill
      user-defined criteria, given in config file
   -  ``probesets``\ folder: sets of non-overlapping probes per gene,
      ranked by best set criteria
   -  ``padlock_probes``\ folder: final padlock probe sequences per
      gene, ready to order

-  ``-d``: optional, ‘download only’ option, where only gene and genome
   annotation files are downloaded but no probes generated, default:
   False

All steps and config parameters will be documented in a log file, that
is saved in the directory where the pipeline is executed from. The
logging file will have the format:
``log_padlock_probe_designer_{year}-{month}-{day}-{hour}-{minute}.txt``.

**Python Import:**

Import padlock probe design pipeline as python package:

::

   import oligo_designer_toolsuite.pipelines.padlock_probe_designer as packlock_probe_designer

   config = './config/padlock_probe_designer.yaml'
   dir_output = './padlock_probes'

   annotations = packlock_probe_designer.download_annotations(config, dir_output)
   packlock_probe_designer.filter_probes(config, annotations, dir_output)
   del annotations # free memory

   packlock_probe_designer.generate_probe_sets(config, dir_output)
   packlock_probe_designer.design_padlock_probes(config, dir_output)
