Introduction
========================

Oligonucleotides (abbrev. oligos) are short, synthetic strands of DNA or
RNA that have many application areas, ranging from research to disease
diagnosis or therapeutics. Oligos can be used as primers during DNA
amplification, as probes for *in situ* hybridization or as guide RNAs
for CRISPR-based gene editing. Based on the intended application and
experimental design, researchers can customize the length, sequence
composition, and thermodynamic properties of the designed oligos.

Various tools exist that provide custom design of oligo sequences
depending on the area of application. Interestingly, all those pipelines
have many common basic processing steps, ranging from the generation of
custom-length oligo sequences, the filtering of oligo sequences based on
thermodynamic properties as well as the selection of an optimal set of
oligos. Despite the fact that most tools apply the same basic processing
steps, each newly developed tool usually uses its own implementation and
different versions of package dependencies for those basic processing
steps. As a consequence, the comparability of tools that differ only in
certain steps is hampered, but also the development of new tools and the
update of existing tools is slowed down, because developers do not have
a common resource for basic functionalities to fall back on. We tackle
this issue by providing such a common resource in our *Oligo Designer
Toolsuite*. This Toolsuite is a collection of modules that provide all
basic functionalities for custom oligo design pipelines within a
flexible Python framework. All modules have a standardized I/O format
and can be combined individually depending on the required processing
steps.

|image0|


.. |image0| image:: ../figures/oligo_design.png
