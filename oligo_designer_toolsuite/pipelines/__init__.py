"""
Design pipelines for generating oligonucleotide sets for different molecular biology assays.

This package collects the ready-to-run pipeline classes and config-driven
workflow functions used by Oligo Designer Toolsuite. Each pipeline designs probe
sets for a specific assay, such as Oligo-seq, SCRINSHOT, MERFISH, seqFISH+,
CycleHCR, or HCR.

The classes expose the individual design steps for users who want more control.
The matching workflow functions run the full pipeline from a configuration file.
"""

from ._cycle_hcr_probe_designer import CycleHCRProbeDesigner, cycle_hcr_probe_designer
from ._genomic_region_generator import GenomicRegionGenerator
from ._hcr_probe_designer import HcrProbeDesigner, hcr_probe_designer
from ._merfish_probe_designer import MerfishProbeDesigner, merfish_probe_designer
from ._oligo_seq_probe_designer import OligoSeqProbeDesigner, oligo_seq_probe_designer
from ._scrinshot_probe_designer import ScrinshotProbeDesigner, scrinshot_probe_designer
from ._seqfish_plus_probe_designer import SeqFishPlusProbeDesigner, seqfish_plus_probe_designer

__all__ = [
    "GenomicRegionGenerator",
    "OligoSeqProbeDesigner",
    "oligo_seq_probe_designer",
    "ScrinshotProbeDesigner",
    "scrinshot_probe_designer",
    "SeqFishPlusProbeDesigner",
    "seqfish_plus_probe_designer",
    "MerfishProbeDesigner",
    "merfish_probe_designer",
    "CycleHCRProbeDesigner",
    "cycle_hcr_probe_designer",
    "HcrProbeDesigner",
    "hcr_probe_designer",
]
