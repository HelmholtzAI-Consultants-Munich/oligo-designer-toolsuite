Oligo-Seq Probe Designer
==========================

An oligo-seq probe is an oligo hybridization probe, which is optimized for probe-based targeted sequencing to measure RNA expression.

If you are using the Oligo-Seq Probe Design Pipeline, consider citing the Oligo Designer Toolsuite package [1].


Command-Line Call
------------------

To create oligo-seq probes you can run the pipeline with 

::

    oligo_seq_probe_designer -c data/configs/oligo_seq_probe_designer.yaml

where:

``-c``: config file, which contains parameter settings, specific to oligo-seq probe design, `oligo_seq_probe_designer.yaml <https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/main/data/configs/oligo_seq_probe_designer.yaml>`__ contains default parameter settings

All steps and config parameters will be documented in a log file, that is saved in the defined output directory. 
The logging file will have the format: ``log_oligo_seq_probe_designer_{year}-{month}-{day}-{hour}-{minute}.txt``.


Python API
------------------

The oligo-seq probe design pipeline can also be integrated directly into Python code.
Below is an example demonstrating how this can be done.
For a complete explanation of all function parameters, refer to the API documentation.

.. code-block:: python
    
    ##### Initialize Probe Designer Pipeline #####
    # We create an instance of the OligoSeqProbeDesigner class. This pipeline handles all steps 
    # required to design probes for oligo-seq experiments, including target probes and final output. 
    # - write_intermediate_steps: whether to save intermediate results (True/False)
    # - dir_output: directory path where output files will be stored
    # - n_jobs: number of CPU cores/threads to use for parallel tasks
    pipeline = OligoSeqProbeDesigner(
        write_intermediate_steps=True,
        dir_output="output_oligo_seq_probe_designer",
        n_jobs=2,
    )

    # Optional: If you need to customize certain developer parameters (for debugging, advanced usage, etc.),
    # call set_developer_parameters(...) with any overrides. By default, the pipeline uses internal defaults.
    pipeline.set_developer_parameters(...)

    ##### Design Target Probes #####
    # We first generate probes that hybridize specifically to target genes sequences.
    # The pipeline will generate multiple candidate sets (n_sets) and return them as part of the probe database.
    oligo_database = pipeline.design_target_probes(
        files_fasta_target_probe_database=...,                  # List of FASTA files with target gene sequences
        files_fasta_reference_database_targe_probe=...,         # List of FASTA files for specificity reference 
        gene_ids=...,                                           # List of gene symbols or identifiers
        target_probe_length_min=26,
        target_probe_length_max=30,
        target_probe_split_region=4,
        target_probe_targeted_exons=[1,2,3],
        target_probe_isoform_consensus=0,
        target_probe_GC_content_min=45,
        target_probe_GC_content_opt=55,
        target_probe_GC_content_max=65,
        target_probe_Tm_min=50,
        target_probe_Tm_opt=60,
        target_probe_Tm_max=70,
        target_probe_secondary_structures_T=37,
        target_probe_secondary_structures_threshold_deltaG=0,
        target_probe_homopolymeric_base_n={"A": 6, "T": 6, "C": 6, "G": 6},
        target_probe_max_len_selfcomplement=10,
        target_probe_hybridization_probability_threshold=0.001,
        target_probe_GC_weight=1,
        target_probe_Tm_weight=1,
        set_size_min=3,
        set_size_opt=5,
        distance_between_target_probes=0,
        n_sets=100,
    )

    ##### Generate Final Output #####
    # Once the oligo database is complete, generate_output() compiles 
    # final reports and stores the top N probe sets. 
    pipeline.generate_output(oligo_database=oligo_database, top_n_sets=3)


Pipeline Description
-----------------------

The pipeline has four major steps:

1) probe generation (dark blue),

2) probe filtering by sequence property and binding specificity (light blue), 

3) probe set selection for each gene (green), and

4) final probe sequence generation (yellow).


.. image:: ../_static/pipeline_oligoseq.jpg
    :width: 500px
    :align: center


For the probe generation step, the user has to provide a FASTA file with genomic sequences which is used as reference for the generation of probe sequences. 
The probe sequences are generated using the ``OligoSequenceGenerator``. 
Therefore, the user has to define the probe length (can be given as a range), and optionally provide a list of gene identifiers (matching the gene identifiers of the annotation file) for which probes should be generated. 
If no gene list is given, probes are generated for all genes in the reference. 
The probe sequences are generated in a sliding window fashion from the DNA sequence of the non-coding strand, assuming that the sequence of the coding strand represents the target sequence of the probe. 
The generated probes are stored in a FASTA file, where the header of each sequence stores the information about its reference region and genomic coordinates. 
In a next step, this FASTA file is used to create an ``OligoDatabase``, which contains all possible probes for a given set of genes. 
When the probe sequences are loaded into the database, all probes of one gene having the exact same sequence are merged into one entry, saving the transcript, exon and genomic coordinate information of the respective probes. 

In the second step, the number of probes per gene is reduced by applying different sequence property (``PropertyFilter``) and binding specificity filters (``SpecificityFilter``). 
For the MERFISH protocol, the following filters are applied: removal of sequences that contain unidentified nucleotides (``HardMaskedSequenceFilter``), that contain low-complexity region like repeat regions (``SoftMaskedSequenceFilter``), that have a GC content (``GCContentFilter``) or melting temperature (``MeltingTemperatureNNFilter``) outside a user-specified range, that contain homopolymeric runs of any nucleotide longer than a user-specified threshold (``HomopolymericRunsFilter``), that contain secondary structures like hairpins below a user-defined free energy threshold (``SecondaryStructureFilter``).
After removing probes with undesired sequence properties from the database, the probe database is checked for probes that potentially cross-hybridize, i.e. probes from different genes that have the exact same or similar sequence. 
Those probes are removed from the database to ensure uniqueness of probes for each gene. 
Cross-hybridizing probes are identified with the ``CrossHybridizationFilter`` that uses a BlastN alignment search to identify similar sequences and removes those hits with the ``RemoveByBiggerRegionPolicy`` that sequentially removes the probes from the genes that have the bigger probe sets. 
Next, the probes are checked for off-target binding with any other region of a provided background reference. 
Off-target regions are sequences of the background reference (e.g. transcriptome or genome) which match the probe region with a certain degree of homology but are not located within the gene region of the probe. 
Those off-target regions are identified with the ``BlastNFilter`` or ``BowtieFilter`` (users choice) and further refined using the ``HybridizationProbabilityFilter`` which calculates the probability of the probe hybridizing to the identified potential off-target sequences.  
Probes with a hybridization probability greater than the user-defined trheshold are removed from the database. Refining the alignment hits with the ``HybridizationProbabilityFilter`` helps to retain more probes in the database.

In the third step of the pipeline, the best sets of non-overlapping probes are identified for each gene. 
The ``OligosetGeneratorIndependentSet`` class is used to generate ranked, non-overlapping probe sets where each probe and probe set is scored according to a protocol dependent scoring function, i.e. by weighted GC content and melting temperature score, of the probes in the set. 
Following this step all genes with insufficient number of probes (user-defined) are removed from the database and stored in a separate file for user-inspection.

In the last step of the pipeline, the ready-to-order probe sequences are reported for the best non-overlapping sets of each gene. 

The output is stored in two separate files: 

- ``oligo_seq_probes.tsv``: contains a table with all generated probes.
- ``oligo_seq_probesets.yml``: contains a detailed description for each probe set, including the sequences of each part of the probe and probe specific attributes.

All default parameters can be found in the `oligo_seq_probe_designer.yaml <https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/main/data/configs/oligo_seq_probe_designer.yaml>`__ config file provided along the repository.


.. [1] Mekki, I., Campi, F., Kuemmerle, L. B., ... & Barros de Andrade e Sousa, L. (2023). Oligo Designer Toolsuite. Zenodo, https://doi.org/10.5281/zenodo.7823048 