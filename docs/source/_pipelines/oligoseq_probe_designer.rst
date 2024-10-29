Oligo-Seq Probe Designer
==========================

An oligo-seq probe is an oligo hybridization probe, which is optimized for probe-based targeted sequencing to measure RNA expression.

Usage
--------

Command-Line Call
^^^^^^^^^^^^^^^^^^^^

To create oligo-seq probes you can run the pipeline with 

::

    oligo_seq_probe_designer -c data/configs/oligo_seq_probe_designer.yaml

where:

``-c``: config file, which contains parameter settings, specific to oligo-seq probe design, *oligo_seq_probe_designer.yaml* contains default parameter settings

All steps and config parameters will be documented in a log file, that is saved in the directory where the pipeline is executed from. 
The logging file will have the format: ``log_oligo_seq_probe_designer_{year}-{month}-{day}-{hour}-{minute}.txt``.

Python API
^^^^^^^^^^^^^^^^^^^

TBD

Pipeline Description
-----------------------

TBD