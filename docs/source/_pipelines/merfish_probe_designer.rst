SeqFISH+ Probe Designer
==========================


A SeqFISH+ probe is a flourescent probe that contains a 28-nt gene-specific sequence complementary to the mRNA, four 15-nt barcode sequences, which are read out by fluorescent secondary readout probes, single T-nucleotide spacers between readout and gene-specific regions, and two 20-nt PCR primer binding sites. 
The specific readout sequences contained by an encoding probe are determined by the binary barcode assigned to that RNA.

Usage
-------

Command-Line Call
^^^^^^^^^^^^^^^^^^^

To create SeqFISH+ probes you can run the pipeline with 

::

    seqfish_plus_probe_designer -c data/configs/seqfish_plus_probe_designer.yaml


where:

``-c``: config file, which contains parameter settings, specific to SeqFISH+ probe design, *seqfish_plus_probe_designer.yaml* contains default parameter settings

All steps and config parameters will be documented in a log file, that is saved in the directory where the pipeline is executed from. 
The logging file will have the format: ``log_seqfish_plus_probe_designer_{year}-{month}-{day}-{hour}-{minute}.txt``.

Python API
^^^^^^^^^^^^^^^^^^^



Pipeline Description
-----------------------

The pipeline has four major steps:

1) probe generation,

2) probe filtering by sequence property and binding specificity, 

3) probe set selection for each gene, and

4) final probe sequence generation.

.. image:: ../_figures/pipeline_merfish.jpg


For the probe generation step, the user has to provide a FASTA file with sequences of genomic regions which is used as reference for the generation of probe sequences. 
The probe sequences are generated using the ``OligoSequenceGenerator``. 
Therefore, the user has to define the probe length (can be given as a range), and optionally provide a list of gene identifiers (matching the gene identifiers of the annotation file) for which probes should be generated. 
If no gene list is given, probes are generated for all genes in the reference. 
The probe sequences are generated in a sliding window fashion from the DNA sequence of the non-coding strand, assuming that the sequence of the coding strand represents the target sequence of the probe. 
The generated probes are stored in a FASTA file, where the header of each sequence stores the information about its reference region and genomic coordinates. 
In a next step, this Fasta file is used to create an ``OligoDatabase``, which is the underlying data structure of the framework that allows combining different filter and selection functionalities in a custom fashion. 
When the probe sequences are loaded into the database, all probes of one gene having the exact same sequence are merged into one entry, saving the transcript, exon and genomic coordinate information of the respective probes. 
Creating the database which contains all possible probes for a given set of genes, concludes the first step of each probe design pipeline. 

In the second step, the number of probes per gene is reduced by applying different sequence properties by applying the ``PropertyFilter`` and binding specificity filters by applying the ``SpecificityFilter``. 
For the MERFISH protocol, the following filters are applied: removal of sequences that contain unidentified nucleotides (``HardMaskedSequenceFilter``), that have a GC content (``GCContentFilter``) or melting temperature (``MeltingTemperatureNNFilter``) outside a user-specified range, that contain homopolymeric runs of any nucleotide longer than a user-specified threshold (``HomopolymericRunsFilter``), that contain secondary structures like hairpins below a user-defined free energy threshold (``SecondaryStructureFilter``).
After removing probes with undesired sequence properties from the database, the probe database is checked for probes that potentially cross-hybridize, i.e. probes from different genes that have the exact same or similar sequence. 
Those probes are removed from the database to ensure uniqueness of probes for each gene. 
Cross-hybridizing probes are identified with the ``CrossHybridizationFilter`` that uses a BlastN alignment search to identify similar sequences and removes those hits with the ``RemoveByBiggerRegionPolicy`` that sequentially removes the probes from the genes that have the bigger probe sets. 
Next, the probes are checked for off-target binding with any other region of a provided background reference as user-input FASAT file. 
Off-target regions are sequences of the background reference (e.g. transcriptome or genome) which match the probe region with a certain degree of homology but are not located within the gene region of the probe. 
Those off-target regions are identified with the ``BlastNFilter`` that removes probes where a BlastN alignment search found off-target sequence matches with a certain coverage and similarity, for which the user has to define thresholds. 
Filtering the probe database for off-target binding concludes the second step of the probe design pipeline. 

In the third step of the pipeline, the best sets of non-overlapping probes are identified for each gene. 
The ``OligosetGeneratorIndependentSet`` class is used to generate ranked, non-overlapping probe sets where each probe and probe set is scored according to a protocol dependent scoring function, i.e. by the distance to the optimal GC content and melting temperature, weighted by the number of targeted transcripts of the probes in the set. 
The identification of the best scored non-overlapping set of probes for each gene concludes the third step of the pipeline. 
Following this step all genes with insufficient number of probes (user-defined) are removed from the database and stored in a separate file for user-inspection.

In the last step of the pipeline, the ready-to-order probe sequences containing all additional required sequences are designed for the best non-overlapping sets of each gene. 
For the SCRINSHOT protocol, the padlock backbone is added to each probe and for each probe a detection oligo is created, by cropping the probe with even nucleotide removal from both ends, exchanging Thymines to Uracils, and placing the fluorescent dye at the side with the closest Uracil as described in Sountoulidis et al. [1]. 
The output is stored in two seperate files: 

- ``padlock_probes_order.yml``: contains for each probe the sequences of the padlock probe and the detection oligo.
- ``padlock_probes.yml``: contains a detailed description for each probe, including the sequences of each part of the probe and probe specific attributes.

All default parameters can be found in the ``scrinshot_probe_designer.yaml`` config file provided along the repository.


If you are using the SCRINSHOT Probe Design Pipeline, consider citing in addition Kuemmerle et al. [2]

.. [1] Sountoulidis, A., Liontos, A., Nguyen, H. P., Firsova, A. B., Fysikopoulos, A., Qian, X., ... & Samakovlis, C. (2020). SCRINSHOT enables spatial mapping of cell states in tissue sections with single-cell resolution. PLoS biology, 18(11), e3000675. https://doi.org/10.1371/journal.pbio.3000675
.. [2] Kuemmerle, L. B., Luecken, M. D., Firsova, A. B., Barros de Andrade e Sousa, L., Straßer, L., Heumos, L., ... & Theis, F. J. (2022). Probe set selection for targeted spatial transcriptomics. bioRxiv, 2022-08. https://doi.org/10.1101/2022.08.16.504115 

