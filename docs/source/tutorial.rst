Toutorial
==================================
This tutorial implements a pipeline to designe Padlock probes using the ``oligo_designer_toolbox`` package.


Padlock Probe Design
--------------------

A padlock probe contains a constant backbone sequence of 53 nucleotides
(nt) and the 5’- and 3’- arms, which are complementary to the
corresponding mRNA sequence. The gene-specific arms of padlock probes
are around 20nt long each, thus the total length of the gene-specific
sequence of each padlock is 40nt.

Define the parameters
`````````````````````

The first thing to do is to define the parameters we want to use to generate the Padlock probes.
A possible way to define all the parameters, that is flexible and reusable, is to use a configuration file.
For this toutorial we will use the YAMl file,
have a look at the end of the tutorial to the configuration file is structured.

Once the configuration file has been set up we have to read its content:

::

    config_file = "../config/padlock_probe_designer.yaml"
    with open(config_file, 'r') as y:
    config = yaml.safe_load(y)
    dir_output = os.path.join(os.path.dirname(os.getcwd()), config["dir_output"]) # create the complete path for the output directory

Oligo sequences generation
``````````````````````````

Now we can start to actually build the pipeline, we will start by generating all the possible probes with length between the maximum and minimum value given belonging to the genes defined in the config file. The probes will be saved in a nested dicionary with the following structure:

[gene][probe_id][probe_feature].

For this we need a ``CustomDB``, or a class that inherits form it (e.g. ``NCBIDB`` and ``EnsemblDB``) and call the method ``create_oligos_DB``. These classes differ on how the the fasta and the gtf filesused to compute the sequences are obtained. The first one oses local files while the others dowload them form the NCBI or Ensambl ftp server.

::

	from oligo_designer_toolsuite.IO import CustomDB, NcbiDB, EnsemblDB

	# define the database class
	if config["source"] == "ncbi":
	    # dowload the fasta files formthe NCBI server
	    database = NcbiDB(
		probe_length_min=config["probe_length_min"],
		probe_length_max=config["probe_length_max"],
		species=config["species"],
		annotation_release=config["annotation_release"],
		n_jobs=config["n_jobs"],
		dir_output=dir_output,
		min_probes_per_gene=config["min_probes_per_gene"],
		)
	elif config["source"] == "ensembl":
	    # dowload the fasta files formthe NCBI server
	    database = EnsemblDB(
		probe_length_min=config["probe_length_min"],
		probe_length_max=config["probe_length_max"],
		species=config["species"],
		annotation_release=config["annotation_release"],
		n_jobs=config["n_jobs"],
		dir_output=dir_output,
		min_probes_per_gene=config["min_probes_per_gene"],
		)
	elif config["source"] == "custom":
	    # use already dowloaded files
	    database = CustomDB(
		probe_length_min=config["probe_length_min"],
		probe_length_max=config["probe_length_max"],
		species=config["species"],
		genome_assembly=config["genome_assembly"],
		annotation_release=config["annotation_release"],
		annotation_source=config["annotation_source"],
		file_annotation=config["file_annotation"],
		file_sequence=config["file_sequence"],
		n_jobs=config["n_jobs"],
		dir_output=dir_output,
		min_probes_per_gene=config["min_probes_per_gene"],
		)

	# read the genes file
	with open(config["file_genes"]) as handle:
	    lines = handle.readlines()
	    genes = [line.rstrip() for line in lines]

	#generate the sequences
	database.create_oligos_DB(genes=genes)


Dictionary structure
''''''''''''''''''''

Here is an example of how the nested dictionary is structured. We print a very simple dictionary with 1 probe.

::

	{'AGRN':
		{'AGRN_1':
			{'probe_sequence': Seq('AGTCCCGTCCCCGGCGCGGCCCGCGCGCTCCTCCGCCG'),
			 'transcript_id': ['NM_001305275.2:XM_005244749.4:XM_011541429.3:NM_198576.4:XM_047419836.1'],
			 'exon_id': ['NM_001305275.2_exon1:XM_005244749.4_exon1:XM_011541429.3_exon1:NM_198576.4_exon1:XM_047419836.1_exon1'],
			 'chromosome': '1',
			 'start': [1020119],
			 'end': [1020157],
			 'strand': '+',
			 'length': 38
			}
		}
	}

Read and write
''''''''''''''

These classes deal with everything that is related with the management of the dataset. In particular, beyond creatig the dataset, they can also read and write the oligo sequences in a **tsv** or **gtf** fromat. The methods ``read_oligos_DB`` and ``write_oligos_DB`` have exactly this purpose.

Therefore, it is possible to save the current state of the dictionary during the pipeline and to retrive form a previous stage if an error uccurred.

::

	if config["write_intermediate_steps"]:
	    database.write_oligos_DB(format=config["file_format"], dir_oligos_DB="oligos_creation")


Property filter
```````````````

Once all the possible sequences are created, we apply a first filtering process based on the sequences properties (e.g. melting temperature or GC content). This is useful to reduce the amount of sequences we have to deal with in the next stages and discard all the sequences that are not suited for the experiment scope.

Each property filter is a calss that inherits from the Abstact Base Class ``PreFilterBase`` They have a method called ``apply`` that takes the ``oligos_DB`` and returns it filtered. To make this process smooth and modular the class ``PropertyFilter`` allows to apply several filters one after the other. It takes in input a list of filter classes and applies them sequentailly to the ``oligos_DB`` returning the final filterd version of the database. Additionally, all the necessary sequence features computed by the filters are stored in the ``oligos_DB`` for possible later use.

To create new property filters follow the Abstact Base Class requirements in ``PreFilterBase``.


::

	from oligo_designer_toolsuite.oligo_property_filter import (
	    PropertyFilter,
	    MaskedSequences,
	    GCContent,
	    MeltingTemperature,
	    PadlockArms
	)

	# the melting temperature params need to be preprocessed
	Tm_params = config["Tm_parameters"]["shared"].copy()
	Tm_params.update(config["Tm_parameters"]["property_filter"])
	Tm_params["nn_table"] = getattr(mt, Tm_params["nn_table"])
	Tm_params["tmm_table"] = getattr(mt, Tm_params["tmm_table"])
	Tm_params["imm_table"] = getattr(mt, Tm_params["imm_table"])
	Tm_params["de_table"] = getattr(mt, Tm_params["de_table"])

	Tm_correction_param = config["Tm_correction_parameters"]["shared"].copy()
	Tm_correction_param.update(config["Tm_correction_parameters"]["property_filter"])

	# initialize the filters clasees
	masked_sequences = MaskedSequences()
	gc_content = GCContent(GC_content_min=config["GC_content_min"], GC_content_max=config["GC_content_max"])
	melting_temperature = MeltingTemperature(
	    Tm_min=config["Tm_min"],
	    Tm_max=config["Tm_max"],
	    Tm_parameters=Tm_params,
	    Tm_correction_parameters=Tm_correction_param
	)
	padlock_arms = PadlockArms(
	    min_arm_length=config["min_arm_length"],
	    max_arm_Tm_dif=config["max_arm_Tm_dif"],
	    arm_Tm_min=config["arm_Tm_min"],
	    arm_Tm_max=config["arm_Tm_max"],
	    Tm_parameters=Tm_params,
	    Tm_correction_parameters=Tm_correction_param,
	)
	# create the list of filters
	filters = [masked_sequences, gc_content, melting_temperature, padlock_arms]

	# initialize the property filter class
	property_filter = PropertyFilter(filters=filters, write_genes_with_insufficient_probes=config["write_removed_genes"])
	# filter the database
	database = property_filter.apply(database=database, n_jobs=config["n_jobs"])
	# write the intermediate result in a file
	if config["write_intermediate_steps"]:
	    database.write_oligos_DB(format=config["file_format"], dir_oligos_DB="property_filter")


Specificity filters
```````````````````

Generally, in experiments using DNA sequences one of the main problems that can occur are off-target binding of the oligo sequences designed. To avoid this one can decide to delete all the sequences that also match regions of the DNA outside the gene they belong to.

The classes in the subpackage ``oligo_speificity_filters`` detect these probes using aligne methods such as Blast and Bowtie and remove them from the database. The currently implemeted classes are: ``ExactMatches``, ``Blastn``, ``Bowtie``, ``Bowtie2``, ``BowtieSeedRegion``, look at the documentation to understand more in detail their features. As before a second class ``SpecificityFilter`` takes a list of all the filters we want to apply, and applies them sequentially to the ``oligos_DB``.

For our pipeline, we will use ``ExactMatches``, ``Blastn``, ``BowtieSeedRegion``, where for the latter the seed region will be generated with ``LigationRegionCreation`` (look at the documentation to understand what seed region means).

However, alignement methods need a reference fasta file to detect the off-target regions. The ``Custom_DB`` provides also the possibility to generate this reference region with the method ``create_reference_DB``.

**Remark:** in future versions the database classes will be split, one will deal only with the oligo sequences, while the other only the the reference.


::

	from oligo_designer_toolsuite.oligo_specificity_filter import (
	    SpecificityFilter,
	    ExactMatches,
	    LigationRegionCreation,
	    BowtieSeedRegion,
	    Blastn,
	)

	dir_specificity = os.path.join(dir_output, "specificity_temporary") # folder where the temporary files will be written

	# generate the reference
	database.create_reference_DB() # use standard parameters

	# intialize the filter classes
	exact_mathces = ExactMatches(dir_specificity=dir_specificity)
	seed_ligation = LigationRegionCreation(ligation_region_size=config["ligation_region_size"])
	seed_region = BowtieSeedRegion(dir_specificity=dir_specificity, seed_region_creation=seed_ligation)
	blastn = Blastn(
	    dir_specificity=dir_specificity,
	    word_size=config["word_size"],
	    percent_identity=config["percent_identity"],
	    coverage=config["coverage"],
	)
	filters = [exact_mathces, seed_region, blastn]

	# initialize the specificity filter class
	specificity_filter = SpecificityFilter(filters=filters, write_genes_with_insufficient_probes=config["write_removed_genes"])
	# filte r the database
	database = specificity_filter.apply(database=database, n_jobs=config["n_jobs"])
	# write the intermediate result
	if config["write_intermediate_steps"]:
	    database.write_oligos_DB(format=config["file_format"], dir_oligos_DB="specificity_filter")



Probeset generation
```````````````````

In the next step of the pipeline the probes will be choosen according to their theoretical efficiency in the experiment scope (e.g. how well they bind to the target in the DNA). Each probe will receive a score computed by a class that inherits from ``ProbeScoringBase``. Later, the sequences will be organized in sets and a class inheriting from ``SetScoringBase`` will give a general efficiency score to the set. At the end the best sets will be selected and sored.

It is required that the each array of sequences  contains probes that do not overlap. In fact, if two probes were overlapping, they would compete to bind to the same section of DNA and their efficiency would drop significantly. Therefore, it is extremely important to consider only sets of non-overlapping sequences.

The class ``ProbesetGenerator`` takes the scoring strategies and tries to find, among all the feasible non-overlapping sets of probes, the sets with the best efficiency scores. These sets will be save in a pandas DataFrame with the following structure:

+-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
| probeset_id | probe_0  | probe_1  | probe_2  |  ...  | probe_n  | set_score_1 | set_score_2 |  ...  |
+-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
| 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      | 1.2332      |  ...  |
+-------------+----------+----------+-----+----+-------+----------+-------------+-------------+-------+


::

	from oligo_designer_toolsuite.oligo_efficiency import(
	    PadlockProbeScoring,
	    PadlockSetScoring,
	)
	from oligo_designer_toolsuite.oligo_selection import ProbesetGenerator, padlock_heuristic_selection

	# initialize the scoring classes
	probes_scoring = PadlockProbeScoring(
	    Tm_min=config["Tm_min"],
	    Tm_opt=config["Tm_opt"],
	    Tm_max=config["Tm_max"],
	    GC_content_min=config["GC_content_min"],
	    GC_content_opt=config["GC_content_opt"],
	    GC_content_max=config["GC_content_max"],
	    Tm_weight=config["Tm_weight"],
	    GC_weight=config["GC_weight"],
	)
	set_scoring = PadlockSetScoring()

	# initialize the probeset generator class
	probeset_generator = ProbesetGenerator(
	    probeset_size=config["probeset_size"],
	    min_probeset_size=config["min_probeset_size"],
	    probes_scoring=probes_scoring,
	    set_scoring=set_scoring,
	    heurustic_selection=padlock_heuristic_selection,
	    write_genes_with_insufficient_probes=config["write_removed_genes"]
	)

	# generate the probeset
	database = probeset_generator.apply(database=database, n_sets=config["n_sets"], n_jobs=config["n_jobs"])
	# write the intermediate result
	if config["write_intermediate_steps"]:
	    database.write_probesets(dir_probesets="probesets")


Last step
`````````

	Once the best probesets are generated each experiment design might require (or not) an addtional step. In the case of the Padlock probe designer the last step consists in designing the final padlock probe sequences.

::

	from oligo_designer_toolsuite.experiment_specific import PadlockSequenceDesigner

	# preprocessing of themelting temperature parameters
	Tm_params = config["Tm_parameters"]["shared"].copy()
	Tm_params.update(config["Tm_parameters"]["detection_oligo"])
	Tm_params["nn_table"] = getattr(mt, Tm_params["nn_table"])
	Tm_params["tmm_table"] = getattr(mt, Tm_params["tmm_table"])
	Tm_params["imm_table"] = getattr(mt, Tm_params["imm_table"])
	Tm_params["de_table"] = getattr(mt, Tm_params["de_table"])

	Tm_correction_param = config["Tm_correction_parameters"]["shared"].copy()
	Tm_correction_param.update(config["Tm_correction_parameters"]["detection_oligo"])

	# initilize the padlock sequence designer class
	padlock_sequence_designer = PadlockSequenceDesigner(
	    detect_oligo_length_min=config["detect_oligo_length_min"],
	    detect_oligo_length_max=config["detect_oligo_length_max"],
	    detect_oligo_Tm_opt=config["detect_oligo_Tm_opt"],
	    Tm_parameters=Tm_params,
	    Tm_correction_parameters=Tm_correction_param
	)
	# generate the padlock sequence
	padlock_sequence_designer.design_padlocks(database=database)


Configuration File
------------------

::

	### General parameters
	dir_output: output # name of the directory where the output files will be written
	n_jobs: 12 # number of cores used to run the pipeline
	write_removed_genes: True # write in a file the removed genes
	write_intermediate_steps: True # writes the oligo sequences after each step of the pipeline
	file_format: tsv # fromat of the written files, can be "tas" or "gtf"


	### Parameters for genome and gene annotation
	# Which annotation should be used? Ensemble or NCBI?
	source: ncbi # available sources: 'ncbi', 'ensemble' or 'custom', if you use custom please provide file_gene_gtf and file_genome_fasta
	species: human # available species: human or mouse
	annotation_release: current # release number (e.g. 109 or 109.20211119 for ncbi) of annotation or 'current' to use most recent annotation release. Check out release numbers for NCBI at ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/ or for ensemble at ftp.ensembl.org/pub/
	min_probes_per_gene: 0 # minimum number of probes that a gene must have before it gets deleted
	# If 'custom' was choosen,  provide a gene annotation (gtf) and genome file (fasta), e.g. file_gene_gtf: ./data/Homo_sapiens.GRCh38.104.gtf and the relative information
	file_annotation: # ./data_ncbi/annotations/GCF_000001405.39_GRCh38.p13_genomic.gtf
	file_sequence: # ./data_ncbi/annotations/GCF_000001405.39_GRCh38.p13_genomic.fna
	genome_assembly: GRCh38.p14
	annotation_source: NCBI



	### Parameters for oligo sequences generation
	probe_length_min: 38 #min length of probes
	probe_length_max: 45 #max length of probes
	file_genes: ../data/genes_ncbi_10.txt # file with a list the genes used to generate the oligos sequences, set as None if al the genes are used


	### Parameters for the property filers
	GC_content_min: 40 #minimum GC content of probes
	GC_content_max: 60 #maximum GC content of probes
	Tm_min: 52 #55 #minimum melting temperature of probes
	Tm_max: 67 #63 #maximum melting temperature of probes
	# probe arms
	min_arm_length: 10 #min length of each arm
	max_arm_Tm_dif: 2 #maximum melting temperature difference of both arms
	arm_Tm_min: 38 #41 #minimum melting temperature of each arm (difference shouldn't be higher than 5! But range is not super important, the lower the better)
	arm_Tm_max: 49 #46 #maximum melting temperature of each arm


	### Parameters for the specificity filters
	# Blastn parameters
	word_size: 10 #word size for the blastn seed (exact match to target)
	coverage: 50 #minimum coverage between probes and target sequence, ranging from 0 to 100% (full coverage)
	percent_identity: 80 #maximum similarity between probes and target sequences, ranging from 0 to 100% (no missmatch)
	ligation_region_size: 10 #coverage between probes and target sequence should not span region around ligation site (e.g. ligation_region = 5 would correspond to -4 to +5 nt around ligation site), if ligation_region = 0, omit this requirement




	### Parameters for the oligo efficiency
	# Here also Tm_min, Tm_max, GC_content_min, GC_content_max are used, but have been defined before
	Tm_opt: 60 #60 #optimal melting temperature of probes
	GC_content_opt: 50 #optimal GC content of probes
	Tm_weight: 1 # weight of the Tm of the probe in the efficiency score
	GC_weight: 1 # weight of the GC content of the probe in the efficiency score




	### Parameters for the probesets generation
	probeset_size: 5 #ideal number of probes per probeset
	min_probeset_size: 2 #minimum number of probes per probeset
	n_sets: 100 # maximum number of sets per gene


	### Parameters for the padlock detection oligo design
	detect_oligo_length_min: 18 # min length of detection oligo
	detect_oligo_length_max: 25 # max length of detection oligo
	detect_oligo_Tm_opt: 32 # optimal melting temperature of detection oligo



	### Shared parameters
	# The melting temperature is used in 2 different stages (property filters and padlock detection oligo design), where a few parameters are shared and the others differ.
	# parameters for melting temperature -> for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
	Tm_parameters:
	    shared:
		check: True
		strict: True
		c_seq: null
		shift: 0
		nn_table: DNA_NN3
		tmm_table: DNA_TMM1
		imm_table: DNA_IMM1
		de_table: DNA_DE1
		dnac1: 50 #[nM]
		dnac2: 0
		selfcomp: False
		dNTPs: 0
		saltcorr: 7
	    property_filter:
		Na: 1.25 #[mM]
		K: 75 #[mM]
		Tris: 20 #[mM]
		Mg: 10 #[mM]
	    detection_oligo:
		Na: 39 #[mM]
		K: 0 #[mM]
		Tris: 0 #[mM]
		Mg: 0 #[mM]

	Tm_correction_parameters:
	    shared:
		DMSO: 0
		DMSOfactor: 0.75
		fmdfactor: 0.65
		fmdmethod: 1
		GC: null
	    property_filter:
		fmd: 20
	    detection_oligo:
		fmd: 30
