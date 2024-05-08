How to Run Ready-To-Use Pipelines? Commandline Calls
=====================================================

In this section we will show you how to run ready-to-use pipelines to design oligonucleotides for specific experiments.

Probe Design Pipelines for *in-situ* Hybridization
----------------------------------------------------

SCRINSHOT Probe Design
^^^^^^^^^^^^^^^^^^^^^^^^

Usage:
""""""

::

    scrinshot_probe_designer [options]* -o <output_dir>

Main Arguments:
"""""""""""""""

--o  name of the directory where the output files will be written

Optional Arguments:
"""""""""""""""""""

**Parameters for Configuration**

--config    path to config file. All parameters listed below can be given as a simple config.yaml file. You can find examples for different config files in the repository under data/configs.
            When no config file is specified, the pipeline automatically generates a default config file in the background and overwrites the default parameters with the user defined parameters listed below.

**Parameters for annotation loader**

--source    define annotation source (currently supported: ncbi, ensembl and custom)


*NCBI annnotation parameters:*

--taxon                 taxon of the species, valid taxa are: archaea, bacteria, fungi, invertebrate, mitochondrion, plant, plasmid, plastid, protozoa, vertebrate_mammalian, vertebrate_other, viral
--species               species name in NCBI download format, e.g. 'Homo_sapiens' for human; see `ncbi page <https://ftp.ncbi.nlm.nih.gov/genomes/refseq/>`_ for available species name
--annotation_release    release number (e.g. 109 or 109.20211119 for ncbi) of annotation or 'current' to use most recent annotation release. Check out release numbers for NCBI at ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/


*Ensembl annotation parameters:*

--species               species name in ensembl download format, e.g. 'homo_sapiens' for human; see `ensembl page <http://ftp.ensembl.org/pub/release-108/gtf/>`_ for available species names
--annotation_release    release number of annotation, e.g. 'release-108' or 'current' to use most recent annotation release. Check out release numbers for ensemble at ftp.ensembl.org/pub/


*Custom annotation parameters:*

--file_annotation       GTF file with gene annotation
--file_sequence         FASTA file with genome sequence
--files_source          original source of the genomic files -> optional
--species               species of provided annotation, leave empty if unknown -> optional
--annotation_release    release number of provided annotation, leave empty if unknown -> optional
--genome_assembly       genome assembly of provided annotation, leave empty if unknown -> optional

**Parameters for Probe Sequences Database:**

--probe_length_min      minimum length of probes
--probe_length_max      maximum length of probes
--min_probes_per_gene   minimum number of probes that a gene must have before it gets deleted


**Parameters for Property Filters:**

*Parameters for Probe Sequence:*

--GC_content_min        minimum GC content of probes
--GC_content_max        maximum GC content of probes
--Tm_min                minimum melting temperature of probes
--Tm_max                maximum melting temperature of probes

*Parameters for Padlock Arms:*

--min_arm_length        minimum length of each arm
--max_arm_Tm_dif        maximum melting temperature difference of both arms
--arm_Tm_min            minimum melting temperature of each arm (difference shouldn't be higher than 5! But range is not super important, the lower the better)
--arm_Tm_max            maximum melting temperature of each arm

*Parameters for Melting Temperature:*

--Tm_parameters_probe             melting temperature parameters for probe design
--Tm_chem_correction_param_pobe   parameters for chemical correction of melting temperature for probe design


**Parameters for Specificity Filters:**

*BlastN Similarity Filter:*

--blast_word_size           word size for the blastn seed (exact match to target)
--blast_percent_identity    maximum similarity between oligos and target sequences, ranging from 0 to 100% (no missmatch)
--blast_coverage            minimum coverage between oligos and target sequence, ranging from 0 to 100% (full coverage)

*Bowtie Ligation Region filter:*

--ligation_region_size      size of the seed region around the ligation site for bowtie seed region filter


**Parameters for Oligo Efficiency Score:**

--Tm_min                minimum melting temperature of probes
--Tm_max                maximum melting temperature of probes
--Tm_opt                optimal melting temperature of probes
--Tm_weight             weight of the Tm of the probe in the efficiency score
--GC_content_min        minimum GC content of probes
--GC_content_max        maximum GC content of probes
--GC_content_opt        optimal GC content of probes
--GC_weight             weight of the GC content of the probe in the efficiency score

**Parameters for Oligosets Generation:**

--probeset_size_opt     ideal number of oligos per probeset
--probeset_size_min     minimum number of oligos per probeset
--n_sets                maximum number of sets per gene


**Parameters for Padlock Final Sequence Design:**

--detect_oligo_length_min                       minimum length of detection oligo
--detect_oligo_length_max                       maximum length of detection oligo
--detect_oligo_Tm_opt                           optimal melting temperature of detection oligo
--Tm_parameters_detection_oligo                 melting temperature parameters for detection oligo design
--Tm_chem_correction_param_detection_oligo      parameters for chemical correction of melting temperature for detection oligo design



MERFISH Probe Design
^^^^^^^^^^^^^^^^^^^^^^^^

TBD


SeqFISH+ Probe Design
^^^^^^^^^^^^^^^^^^^^^^^^

TBD


Probe Design Pipelines for CRISPR experiments
----------------------------------------------

Coming Soon
