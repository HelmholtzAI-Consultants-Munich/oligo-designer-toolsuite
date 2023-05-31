How to Run Ready-To-Use Pipelines? Python API
==============================================

In this section we will show you how to integrate ready-to-use pipelines to design oligonucleotides for specific experiments into your Python code.

SCRINSHOT Probe Design
-----------------------

Import Packages:

..  code-block:: python

    from oligo_designer_toolsuite.pipelines import ScrinshotProbeDesigner


Probeset Design
^^^^^^^^^^^^^^^^^

To design probesets for a given set of genes, we first create an instance of the ``ScrinshotProbeDesigner`` class. 
We need to define an output directory and can set the parameters ``write_removed_genes`` (if true, save gene with insufficient probes in a file) and
``write_intermediate_steps`` (if true, save the probe database after each processing step, such that the pipline can resumed from a certain step onwards).

..  code-block:: python

    probe_designer = ScrinshotProbeDesigner(dir_output="./output", write_removed_genes=True, write_intermediate_steps = True)

After instatiating the ``ScrinshotProbeDesigner`` class, we need to load the annotation we are using. As an example we will use the NCBI gene annotation. 
Hence, we define *ncbi* as source and define the NCBI-specific parameters *taxon*, *species* and *annotation_release*.
Apart from *NCBI* annotation, we can also choose an *Ensembl* annotation. If ``source="ncbi"`` or ``source="ensembl"`` is choosen, 
the annotation files are automatically downloaded from their servers.
In addition, we can provide a custom annotation when specifying ``sourec="custom"``. 

Parameters for annotation loader
"""""""""""""""""""""""""""""""""""""""""""""
- source: define annotation source (currently supported: ncbi, ensembl and custom)

*NCBI annnotation parameters:*

- ``taxon``: taxon of the species, valid taxa are: archaea, bacteria, fungi, invertebrate, mitochondrion, plant, plasmid, plastid, protozoa, vertebrate_mammalian, vertebrate_other, viral
- ``species``: species name in NCBI download format, e.g. 'Homo_sapiens' for human; see [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/) for available species name
- ``annotation_release``: release number (e.g. 109 or 109.20211119 for ncbi) of annotation or 'current' to use most recent annotation release. Check out release numbers for NCBI at ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/

*Ensembl annotation parameters:*

- ``species``: species name in ensembl download format, e.g. 'homo_sapiens' for human; see http://ftp.ensembl.org/pub/release-108/gtf/ for available species names
- ``annotation_release``: release number of annotation, e.g. 'release-108' or 'current' to use most recent annotation release. Check out release numbers for ensemble at ftp.ensembl.org/pub/

*Custom annotation parameters:*

- ``file_annotation``: GTF file with gene annotation
- ``file_sequence``: FASTA file with genome sequence
- ``files_source``: original source of the genomic files (optional)
- ``species``: species of provided annotation, leave empty if unknown (optional)
- ``annotation_release``: release number of provided annotation, leave empty if unknown (optional)
- ``genome_assembly``: genome assembly of provided annotation, leave empty if unknown (optional)


..  code-block:: python

    # example for ncbi annotation loader
    source = "ncbi"
    params = {
        "taxon": "vertebrate_mammalian",
        "species": "Homo_sapiens",
        "annotation_release": "110",
    }

    # example for ensembl annotation loader
    # source = "ensembl"
    # params = {
    #     "species": "homo_sapiens",
    #     "annotation_release": "109",
    # }

    # example for custom annotation loader
    # source = "custom"
    # params = {
    #     "file_annotation": "./output/annotation/GCF_000001405.40_GRCh38.p14_genomic.gtf",
    #     "file_sequence": "./output/annotation/GCF_000001405.40_GRCh38.p14_genomic.fna",
    #     "files_source": "NCBI",
    #     "species": "Homo_sapiens",
    #     "annotation_release": "110",
    #     "genome_assembly": "GRCh38.p14",
    # }

    probe_designer.load_annotations(source=source, source_params=params)


After downloading the annotations, we have to create the oligo database. 
Running the function below, will automatically create a transcriptome from the given annotation (therefore, the provided GTF file must contain transcript and exon information) 
and use this transcriptome to create all possible probes for each gene, that is provided in the *gene* list. 

Parameters for Probe Sequences Database
"""""""""""""""""""""""""""""""""""""""""""""

- ``probe_length_min``: minimum length of probes
- ``probe_length_max``: maximum length of probes
- ``min_probes_per_gene``: minimum number of probes that a gene must have before it gets deleted

..  code-block:: python

    probe_length_min = 38
    probe_length_max = 45
    min_probes_per_gene = 3
    genes = ...

    probe_database, file_database = probe_designer.create_probe_database(genes=genes, probe_length_min=probe_length_min, probe_length_max=probe_length_max, min_probes_per_gene=min_probes_per_gene, n_jobs=4)

*Note: Instead of creating a new probe database, we can also load an existing databases.*  

Loading a database can be useful when starting the pipeline from a certain step, e.g. load a database which was already filtered by probe properties and continue immediately with the specificity filter step. 
We can load an existing database by calling ``load_probe_database()``. See example code in the cells below (commented).

In order to create experiment-specific probes, we have to apply several filter to each probe, e.g. melting temperature or GC content filters. 

Parameters for Property Filters
"""""""""""""""""""""""""""""""""""""""""""""

*Parameters for Probe Sequence:*

- ``GC_content_min``: minimum GC content of probes
- ``GC_content_max``: maximum GC content of probes
- ``Tm_min``: minimum melting temperature of probes
- ``Tm_max``: maximum melting temperature of probes

*Parameters for Padlock Arms:*

- ``min_arm_length``: minimum length of each arm
- ``max_arm_Tm_dif``: maximum melting temperature difference of both arms
- ``arm_Tm_min``: minimum melting temperature of each arm (difference shouldn't be higher than 5! But range is not super important, the lower the better)
- ``arm_Tm_max``: maximum melting temperature of each arm

*Parameters for Melting Temperature:*

- ``Tm_parameters_probe``: melting temperature parameters for probe design
- ``Tm_chem_correction_param_pobe``: parameters for chemical correction of melting temperature for probe design

*Note: The melting temperature is used in 2 different stages (probe and detection oligo design), where a few parameters are shared and the others differ. 
Parameters for melting temperature - for more information on parameters, see:* `here <https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN>`_

..  code-block:: python

    ####### Load existing database #######
    # file_database = "./output/oligo_database/probe_database_initial.txt"
    # min_probes_per_gene = 3
    # probe_database = probe_designer.load_probe_database(file_database=file_database, min_probes_per_gene=min_probes_per_gene)

    ####### Apply Property Filter #######
    GC_content_min=40
    GC_content_max=60
    Tm_min=52
    Tm_max=67
    min_arm_length=10
    max_arm_Tm_dif=2
    arm_Tm_min=38
    arm_Tm_max=49

    probe_database, file_database = probe_designer.filter_probes_by_property(probe_database, GC_content_min=GC_content_min, GC_content_max=GC_content_max,
                                                                            Tm_min=Tm_min, Tm_max=Tm_max, min_arm_length=min_arm_length, max_arm_Tm_dif=max_arm_Tm_dif, arm_Tm_min=arm_Tm_min, arm_Tm_max=arm_Tm_max, n_jobs=4)



Parameters for Specificity Filters
"""""""""""""""""""""""""""""""""""""""""""""

*BlastN Similarity Filter:*

- ``blast_word_size``: word size for the blastn seed (exact match to target)
- ``blast_percent_identity``: maximum similarity between oligos and target sequences, ranging from 0 to 100% (no missmatch)
- ``blast_coverage``: minimum coverage between oligos and target sequence, ranging from 0 to 100% (full coverage)

*Bowtie Ligation Region filter:*

- ``ligation_region_size``: size of the seed region around the ligation site for bowtie seed region filter

*Note: Depending on the number of genes, this step might be time and memory consuming. For high number of genes, you might want to run this step on a bigger machine!*


..  code-block:: python

    ####### Load existing database #######
    # load annotation files for Reference Database
    # source = "custom"
    # custom_params = {
    #     "file_annotation": "./output/annotation/GCF_000001405.40_GRCh38.p14_genomic.gtf",
    #     "file_sequence": "./output/annotation/GCF_000001405.40_GRCh38.p14_genomic.fna",
    #     "files_source": "NCBI",
    #     "species": "Homo_sapiens",
    #     "annotation_release": "110",
    #     "genome_assembly": "GRCh38.p14",    
    # }
    # probe_designer.load_annotations(source=source, source_params=custom_params)

    # load existing database
    # file_database = "./output/oligo_database/probe_database_property_filter.txt"
    # min_probes_per_gene = 3
    # probe_database = probe_designer.load_probe_database(file_database=file_database, min_probes_per_gene=min_probes_per_gene)

    ####### Apply Specificity Filter #######
    ligation_region_size=5
    blast_word_size=10
    blast_percent_identity=80
    blast_coverage=50

    probe_database, file_database = probe_designer.filter_probes_by_specificity(probe_database, ligation_region_size=ligation_region_size, 
                                                                                blast_word_size=blast_word_size, blast_percent_identity=blast_percent_identity, blast_coverage=blast_coverage, n_jobs=2)
        
After applying different sets of filters to the probe database, we will create probesets for each gene, 
which are sets of probes that do not overlap and have a high efficiency score (calculated from melting temperature and GC content).

Parameters for Oligo Efficiency Score
"""""""""""""""""""""""""""""""""""""""""""""

- ``Tm_min``: minimum melting temperature of probes
- ``Tm_max``: maximum melting temperature of probes
- ``Tm_opt``: optimal melting temperature of probes
- ``Tm_weight``: weight of the Tm of the probe in the efficiency score
- ``GC_content_min``: minimum GC content of probes
- ``GC_content_max``: maximum GC content of probes
- ``GC_content_opt``: optimal GC content of probes
- ``GC_weight``: weight of the GC content of the probe in the efficiency score

Parameters for Oligosets Generation
"""""""""""""""""""""""""""""""""""""""""""""

- ``probeset_size_opt``: ideal number of oligos per probeset
- ``probeset_size_min``: minimum number of oligos per probeset
- ``n_sets``: maximum number of sets per gene


..  code-block:: python

    ####### Load existing database #######
    # file_database = "./output/oligo_database/oligo_database_specificity_filters.txt"
    # min_probes_per_gene = 3
    # probe_database = probe_designer.load_probe_database(file_database=file_database, min_probes_per_gene=min_probes_per_gene)

    ####### Apply Probe Set Selection #######
    probeset_size_opt=5
    probeset_size_min=2
    n_sets=100
    Tm_min=52
    Tm_max=67
    Tm_opt=60
    Tm_weight=1
    GC_content_min=40
    GC_content_max=60
    GC_content_opt=50
    GC_weight=1

    probe_database, file_database, dir_oligosets = probe_designer.create_probe_sets(probe_database, 
                                                                                    probeset_size_opt=probeset_size_opt, 
                                                                                    probeset_size_min=probeset_size_min, 
                                                                                    n_sets=n_sets, 
                                                                                    Tm_min=Tm_min, 
                                                                                    Tm_max=Tm_max, 
                                                                                    Tm_opt=Tm_opt, 
                                                                                    Tm_weight=Tm_weight, 
                                                                                    GC_content_min=GC_content_min, 
                                                                                    GC_content_max=GC_content_max, 
                                                                                    GC_content_opt=GC_content_opt, 
                                                                                    GC_weight=GC_weight, 
                                                                                    n_jobs=2)

In the probe database, the gene names are the keys of the database. All genes that do not have sufficient probes were removed from the database.  
Once we hve all genes with sufficient probes, we create the final "read to order" probe sequences. 
Calling the fuction below will produce two files, *padlock_probes* and *padlock_probes_order*.
The latter file contains the ready to order probe sequences for each gene.

Parameters for Padlock Final Sequence Design
"""""""""""""""""""""""""""""""""""""""""""""

- ``detect_oligo_length_min``: minimum length of detection oligo
- ``detect_oligo_length_max``: maximum length of detection oligo
- ``detect_oligo_Tm_opt``: optimal melting temperature of detection oligo
- ``Tm_parameters_detection_oligo``: melting temperature parameters for detection oligo design
- ``Tm_chem_correction_param_detection_oligo``: parameters for chemical correction of melting temperature for detection oligo design

*Note: The melting temperature is used in 2 different stages (probe and detection oligo design), where a few parameters are shared and the others differ. 
Parameters for melting temperature - for more information on parameters, see:* `here <https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN>`_

..  code-block:: python

    ##### Design final sequences #####
    detect_oligo_length_min = 18
    detect_oligo_length_max = 25
    detect_oligo_Tm_opt = 32

    probe_designer.create_final_sequences(probe_database, detect_oligo_length_min, detect_oligo_length_max, detect_oligo_Tm_opt)

