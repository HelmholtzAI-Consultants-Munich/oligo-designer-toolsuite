How to Build a Custom Oligo Design Pipeline
===========================================

The Oligo Designer Toolsuite is a collection of modules providing core functionalities 
for custom oligo design pipelines within a flexible Python framework. 
All modules share a common runtime- and memory-optimized data structure 
and a standardized API, which allows users to combine them flexibly based on their 
required processing steps.

This section demonstrates how to use the Oligo Designer Toolsuite to develop 
a custom oligo design pipeline.

Setup
-----
Before getting started, you must install the Oligo Designer Toolsuite Python package. 
A step-by-step installation guide is available 
`here <https://oligo-designer-toolsuite.readthedocs.io/en/latest/_getting_started/installation.html>`_.

Oligo Sequence Generation
-------------------------
The first step is generating oligo sequences. 
`This notebook <_tutorials/1-oligo-sequences-generation>`_ describes how to produce oligo sequences either from a reference genomic sequence 
or at random.

Oligo Database
--------------
The oligo database is a central data structure in the Oligo Designer Toolsuite. 
It stores generated oligos and is used as both input and output for all pipeline steps. 
`This notebook <_tutorials/2-oligo-database>`_ shows how to create a database from previously generated sequences 
(in FASTA format) and demonstrates various functionalities, such as pre-filtering by attributes, 
loading, saving, and exporting.

Property Filters
----------------
The next step in a typical pipeline involves applying property filters. 
`This notebook <_tutorials/3-property-filters>`_ illustrates how to define and apply property filters on an oligo database.

Specificity Filters
-------------------
Specificity filters ensure that designed oligos do not bind to off-target sites. 
They operate with both the ``OligoDatabase`` and ``ReferenceDatabase`` objects. 
`This notebook <_tutorials/4-specificity-filters>`_ demonstrates how to create a reference database and apply specificity filters.

Oligo Sets Generation
---------------------
After applying property and specificity filters, the final phase involves selecting 
the best-performing sets of oligos. `This notebook <_tutorials/5-oligoset-generation>`_ shows how to score oligos and 
use a graph-based approach to identify non-overlapping, optimal oligo sets.