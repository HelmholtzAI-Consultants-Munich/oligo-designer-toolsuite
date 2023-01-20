User Guide: How to build your own guide design pipeline?
========================================================

In this section we will show you how to use the Oligo Designer Toolsuite to develop your own custom oligo design pipeline.

The Oligo Designer Toolsuite is a collection of modules that provide basic functionalities for custom oligo design pipelines
within a flexible Python framework. All modules have a standardized API and can be combined individually
depending on the required processing steps. Each step takes as input a database containg oligo sequences and
returns a database with the same structure.

For more detailed infomration have a look at one of the tutorials.

Oligos DB
---------

The class that contains the database of the oligos and all the related functionalities (e.g. read, witrite, ...) is ''database.CustomOligoDB''.
It creates the olios sequences for a given set of genes and stores them in a dictionary (``oligos_DB``) with a JSON format. [UPDATE]

The dictionary contains all the oligos and their features, the structure is the following:

::

    {"gene":
    	{"ID":
    		{"Sequence": "ACAGCGCGCCG"
    		 "transcript_id": ["XM_006710333.4", ...]
    		 "exon_id": ["XM_006710333.4_exon2_XM_006710333.4_exon1", ...]
    		 "chromosome": "1"
    		 "start": [1267959, ...]
    		 "end": [1267997, ...]
    		 "strand": "-"
    		 "length": 38
    		 ...
    		 "additional features": value
    		 "GC_content": 52.0 # Example of additional features
     		 }
    	}
    }


Reference DB
------------

TBD


Working principle
-----------------

On a high level, the steps of the pipeline correspond to the modules of the package, in particular:

- **Database**: generation of the oligo sequences and reference sequence database

- **Property Filters**: oligo selection based on specific oligo properties (e.g. 52.0 < Melting Temperature < 57.0 )

- **Specificity Filters**: filering of oligos with high off-target hits using alignement methods such as Bowtie and Blast

- **Efficiency Filters**: filtering of oligos with low efficiency, based on experiment specific scoring methods

- **Selection**: generation of the best sets of oligos based on an application-spercific scoring function


Each step of the pipeline is made of:

- **Application-Specific Modules**: experiment specific functionalities (e.g. filters for GC content)

- **General Module**: modules that combine the application-specific modules in a modular way. They take the database containing the oligo sequences and the experiment specific functionalities and apply the latter to the database.


At the beginning all the possible oligo sequences are created and stored
in the database class (``CustomDB``, ``NCBIDB``, ``EnsembleDB``) in a dictionary (``oligos_DB``). [UPDATE]

At each step a step-specific module will take as input the database class, perform the necessary computations and update the ``oligos_DB``. [UPDATE]

Output
------

The pipeline returns the best set of sequences for each requested gene. They are stored in a ``pandas.DataFrame`` with the following structure:

+-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
| oligoset_id | oligo_0  | oligo_1  | oligo_2  |  ...  | oligo_n  | set_score_1 | set_score_2 |  ...  |
+-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
| 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      | 1.2332      |  ...  |
+-------------+----------+----------+-----+----+-------+----------+-------------+-------------+-------+
