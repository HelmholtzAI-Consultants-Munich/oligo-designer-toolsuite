User Guide
==========

This Tool suite is a collection of modules that provide
basic functionalities for custom oligo design pipelines
within a flexible Python framework. All modules have a
standardized API and can be combined individually
depending on the required processing steps. Each step takes in input a database containg oligo sequences and
returns a database with the same structure.

Oligos DB
---------

The class that contains the database of teh oligos and all the related functionalities (e.g. read, witrite, ...) is ''IO.CustomDB''.
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

On a high level, the steps of the pipeline correspond to teh modules od the package, in particular:

- **IO**: generation of the oligo sequences

- **Pre filters**: selection of the probes  based on feature based filters (e.g. 52.0 < Melting Temperature < 57.0 )

- **Specificity filters**: filering of probes with high off-target locations using alignement methods such as Bowtie and Blast

- **Selection**: generation of the best sets of oligos based on an application-spercific scoring function


Each step of the pipeline is mode of:

- **general module**:

- **application-specific modules**:

The pipeline is split into different steps, each of which has a specific
task. For each step, we have

-  a class containing the basic functionalities in common between all
   the possible diverse experimental designs
-  Other classes that share the same structure and contain the
   functionalities specific to an experimental design are given as input
   to the basic class. The basic structure is explicitly stated in ABC
   classes, hence is possible to develop new functionalities which would
   perfectly integrate in the pipeline.

At the beginning all the possible oligo sequences are created and stored
in a dictionary, then at each step are filtered in order to keep only
the best ones for the purpose.

Passing the whole class allows at each step the write the
intermediate result in a consistent format. Moreover, in case of
necessity, it would be possible to create a new DB class by simply
reading a previous intermediate step and running the pipeline from that
point onwards.

Possible use cases of this functionality are:

-  save time in running two different pipelines which share some of the
   fist steps
-  recover from a broken pipeline

Output
------

The pipeline returns the best set of sequences for each requested gene.
