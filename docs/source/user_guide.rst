User Guide
==========

This Tool suite is a collection of modules that provide
basic functionalities for custom oligo design pipelines
within a flexible Python framework. All modules have a
standardized API and can be combined individually
depending on the required processing steps. Each step takes in input a database containg oligo sequences and
returns a database with the same structure.

For more detailed infomration look at the toutorial in the next section.

Oligos DB
---------

The class that contains the database of the oligos and all the related functionalities (e.g. read, witrite, ...) is ''IO.CustomDB''.
It creates the olios sequences for a given set of genes and stores them in a disctionay (``oligos_DB``) with a JSON format.

The dictionary contains all the probes and their features, the structure is the following:

::

    {"Gene name":
    	{"Sequence ID":
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


Working principle
-----------------

On a high level, the steps of the pipeline correspond to the modules od the package, in particular:

- **IO**: generation of the oligo sequences

- **Pre filters**: selection of the probes  based on feature based filters (e.g. 52.0 < Melting Temperature < 57.0 )

- **Specificity filters**: filering of probes with high off-target locations using alignement methods such as Bowtie and Blast

- **Selection**: generation of the best sets of oligos based on an application-spercific scoring function


Each step of the pipeline is mede of:

- **application-specific modules**: experiment specific functionalities (e.g. filters fro GC content)

- **general module**: modules that combine the application-specific modules in a modular way. They take the database containing the oligo sequences and the experiment specific functionalities and apply the latter to the database.


At the beginning all the possible oligo sequences are created and stored
in the database class (``CustomDB``, ``NCBIDB``, ``EnsembleDB``) in a dictionary (``oligos_DB``).

At each step a step-specific module will take as input the database class, perform the necessary computations and update the ``oligos_DB``.

Output
------

The pipeline returns the best set of sequences for each requested gene. They are stored in a ``pandas.DataFrame`` with the following structure:

+-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
| probeset_id | probe_0  | probe_1  | probe_2  |  ...  | probe_n  | set_score_1 | set_score_2 |  ...  |
+-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
| 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      | 1.2332      |  ...  |
+-------------+----------+----------+-----+----+-------+----------+-------------+-------------+-------+
