"""
The module provides a collection of comprehensive oligonucleotide design pipelines, each equipped with specific functionalities to facilitate efficient and targeted oligo synthesis for diverse applications. 

Classes:
- GenomicRegionGenerator: This class automates the generation of genomic regions from sequence data, streamlining the preparation process for oligonucleotide design.
- OligoSeqProbeDesigner: This class automates the creation of probes for the oligo-seq protocol, integrating various design criteria and computational strategies to optimize probe specificity and efficiency.
"""

from ._genomic_region_generator import GenomicRegionGenerator
from ._oligo_seq_probe_designer import OligoSeqProbeDesigner


__all__ = [
    "GenomicRegionGenerator",
    "OligoSeqProbeDesigner",
]
