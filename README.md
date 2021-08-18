# Probe Design

In a targeted spatial transcriptomics experiment, you need to know a priori what genes you would like to see. These genes are targeted by probes that target mRNA sequences by hybridization. However, we cannot design probes for all genes as some genes are just too similar in sequence to keep apart. A probe for one of these genes would also bind to the other, similar gene and therefore cannot give a specific readout.

This package provides a filter that removes the set of genes for which we cannot design a probe for a targeted spatial transcriptomics experiment.
