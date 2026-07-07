############################################
# imports
############################################

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from abc import ABC, abstractmethod
from typing import Any

import yaml

from oligo_designer_toolsuite.config.pipelines.oligo_seq_probe_designer import OligoSeqProbeDesignerConfig

############################################
# Setup
############################################

SCRIPT_GENOMIC_REGION_GENERATOR = "oligo_designer_toolsuite/pipelines/_genomic_region_generator.py"
CONFIG_GENOMIC_REGION_GENERATOR = "data/configs/genomic_region_generator_custom.yaml"

SCRIPT_OLIGO_SEQ_PROBE_DESIGNER = "oligo_designer_toolsuite/pipelines/_oligo_seq_probe_designer.py"
CONFIG_OLIGO_SEQ_PROBE_DESIGNER = "data/configs/oligo_seq_probe_designer.yaml"

SCRIPT_SCRINSHOT_PROBE_DESIGNER = "oligo_designer_toolsuite/pipelines/_scrinshot_probe_designer.py"
CONFIG_SCRINSHOT_PROBE_DESIGNER = "data/configs/scrinshot_probe_designer.yaml"

SCRIPT_SEQFISHPLUS_PROBE_DESIGNER = "oligo_designer_toolsuite/pipelines/_seqfish_plus_probe_designer.py"
CONFIG_SEQFISHPLUS_PROBE_DESIGNER = "data/configs/seqfish_plus_probe_designer.yaml"

SCRIPT_MERFISH_PROBE_DESIGNER = "oligo_designer_toolsuite/pipelines/_merfish_probe_designer.py"
CONFIG_MERFISH_PROBE_DESIGNER = "data/configs/merfish_probe_designer.yaml"

SCRIPT_CYCLEHCR_PROBE_DESIGNER = "oligo_designer_toolsuite/pipelines/_cycle_hcr_probe_designer.py"
CONFIG_CYCLEHCR_PROBE_DESIGNER = "data/configs/cycle_hcr_probe_designer.yaml"

SCRIPT_HCR_PROBE_DESIGNER = "oligo_designer_toolsuite/pipelines/_hcr_probe_designer.py"
CONFIG_HCR_PROBE_DESIGNER = "data/configs/hcr_probe_designer.yaml"

############################################
# Tests
############################################


class PipelinesBase(unittest.TestCase, ABC):
    def setUp(self) -> None:
        self.tmp_path = self.setup_output_dir()
        self.script = self.setup_script()
        self.cmd_parameters = self.setup_cmd_parameters()

    def tearDown(self) -> None:
        if self.tmp_path is not None and os.path.exists(self.tmp_path):
            shutil.rmtree(self.tmp_path)

    @abstractmethod
    def setup_output_dir(self) -> str:
        pass

    @abstractmethod
    def setup_script(self) -> str:
        pass

    @abstractmethod
    def setup_cmd_parameters(self) -> list[str]:
        pass

    def test_main_function(self) -> None:
        # Run the script using subprocess
        result = subprocess.run(
            [
                sys.executable,
                self.script,
                *self.cmd_parameters,
            ],
            capture_output=True,
            text=True,
        )

        # Check terminal output
        print(result)

        # Check the return code to ensure the script ran successfully
        self.assertEqual(result.returncode, 0)


class TestGenomicRegionGenerator(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self) -> Any:
        with open(CONFIG_GENOMIC_REGION_GENERATOR, "r") as handle:
            config = yaml.safe_load(handle)
        dir_output = config.get("dir_output")
        if dir_output is None:
            return tempfile.mkdtemp()
        return os.path.abspath(dir_output)

    def setup_script(self) -> str:
        return os.path.abspath(SCRIPT_GENOMIC_REGION_GENERATOR)

    def setup_cmd_parameters(self) -> list[str]:
        return ["-c", os.path.abspath(CONFIG_GENOMIC_REGION_GENERATOR)]


class TestOligoSeqProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self) -> Any:
        with open(CONFIG_OLIGO_SEQ_PROBE_DESIGNER, "r") as handle:
            config = OligoSeqProbeDesignerConfig(**yaml.safe_load(handle))
        dir_output = config.general.dir_output
        if dir_output is None:
            return tempfile.mkdtemp()
        return os.path.abspath(dir_output)

    def setup_script(self) -> str:
        return os.path.abspath(SCRIPT_OLIGO_SEQ_PROBE_DESIGNER)

    def setup_cmd_parameters(self) -> list[str]:
        return ["-c", os.path.abspath(CONFIG_OLIGO_SEQ_PROBE_DESIGNER)]


class TestScrinshotProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self) -> Any:
        with open(CONFIG_SCRINSHOT_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        dir_output = config.get("general", {}).get("dir_output")
        if dir_output is None:
            return tempfile.mkdtemp()
        return os.path.abspath(dir_output)

    def setup_script(self) -> str:
        return os.path.abspath(SCRIPT_SCRINSHOT_PROBE_DESIGNER)

    def setup_cmd_parameters(self) -> list[str]:
        return ["-c", os.path.abspath(CONFIG_SCRINSHOT_PROBE_DESIGNER)]


class TestSeqfishplusProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self) -> Any:
        with open(CONFIG_SEQFISHPLUS_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        dir_output = config.get("dir_output")
        if dir_output is None:
            return tempfile.mkdtemp()
        return os.path.abspath(dir_output)

    def setup_script(self) -> str:
        return os.path.abspath(SCRIPT_SEQFISHPLUS_PROBE_DESIGNER)

    def setup_cmd_parameters(self) -> list[str]:
        return ["-c", os.path.abspath(CONFIG_SEQFISHPLUS_PROBE_DESIGNER)]


class TestMerfishProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self) -> Any:
        with open(CONFIG_MERFISH_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        dir_output = config.get("general", {}).get("dir_output")
        if dir_output is None:
            return tempfile.mkdtemp()
        return os.path.abspath(dir_output)

    def setup_script(self) -> str:
        return os.path.abspath(SCRIPT_MERFISH_PROBE_DESIGNER)

    def setup_cmd_parameters(self) -> list[str]:
        return ["-c", os.path.abspath(CONFIG_MERFISH_PROBE_DESIGNER)]


class TestCycleHCRProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self) -> Any:
        with open(CONFIG_CYCLEHCR_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        dir_output = config.get("general", {}).get("dir_output")
        if dir_output is None:
            return tempfile.mkdtemp()
        return os.path.abspath(dir_output)

    def setup_script(self) -> str:
        return os.path.abspath(SCRIPT_CYCLEHCR_PROBE_DESIGNER)

    def setup_cmd_parameters(self) -> list[str]:
        return ["-c", os.path.abspath(CONFIG_CYCLEHCR_PROBE_DESIGNER)]


class TestHCRProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self) -> Any:
        with open(CONFIG_HCR_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        dir_output = config.get("general", {}).get("dir_output")
        if dir_output is None:
            return tempfile.mkdtemp()
        return os.path.abspath(dir_output)

    def setup_script(self) -> str:
        return os.path.abspath(SCRIPT_HCR_PROBE_DESIGNER)

    def setup_cmd_parameters(self) -> list[str]:
        return ["-c", os.path.abspath(CONFIG_HCR_PROBE_DESIGNER)]
