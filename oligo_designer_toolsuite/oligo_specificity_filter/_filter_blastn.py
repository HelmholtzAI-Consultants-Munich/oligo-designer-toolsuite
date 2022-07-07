import os
from pathlib import Path

from ._filter_base import ProbeFilterBase


class ProbeFilterBlastn(ProbeFilterBase):
    """This class filters probes based on the blast alignment tool.

    :param"""

    def __init__(
        self, word_size, percent_identity, probe_length_min, probe_length_max, coverage
    ):
        super().__init__()

        self.word_size = word_size
        self.percent_identity = percent_identity
        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max
        self.coverage = coverage

        self.dir_blast = os.path.join(self.dir_output, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_probes = os.path.join(self.dir_output, "probes")
        Path(self.dir_probes).mkdir(parents=True, exist_ok=True)

        self.file_removed_genes = os.path.join(
            self.dir_output, "genes_with_insufficient_probes.txt"
        )
