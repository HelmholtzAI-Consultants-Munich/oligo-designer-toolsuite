############################################
# imports
############################################

import os
import shutil
import subprocess
import sys
import unittest
from abc import abstractmethod

import yaml

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

############################################
# Tests
############################################


class PipelinesBase:
    def setUp(self):
        self.tmp_path = self.setup_output_dir()
        self.script = self.setup_script()
        self.cmd_parameters = self.setup_cmd_parameters()

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    @abstractmethod
    def setup_output_dir(self):
        pass

    @abstractmethod
    def setup_script(self):
        pass

    @abstractmethod
    def setup_cmd_parameters(self):
        pass

    def test_main_function(self):
        # Run the script using subprocess
        result = subprocess.run(
            [
                sys.executable,
                self.script,
                self.cmd_parameters,
            ],
            capture_output=True,
            text=True,
        )

        # Check terminal output
        print(result)

        # Check the return code to ensure the script ran successfully
        self.assertEqual(result.returncode, 0)


class TestGenomicRegionGenerator(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self):
        with open(CONFIG_GENOMIC_REGION_GENERATOR, "r") as handle:
            config = yaml.safe_load(handle)
        return os.path.abspath(config["dir_output"])

    def setup_script(self):
        return os.path.abspath(SCRIPT_GENOMIC_REGION_GENERATOR)

    def setup_cmd_parameters(self):
        return f"-c{os.path.abspath(CONFIG_GENOMIC_REGION_GENERATOR)}"


class TestOligoSeqProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self):
        with open(CONFIG_OLIGO_SEQ_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        return os.path.abspath(config["dir_output"])

    def setup_script(self):
        return os.path.abspath(SCRIPT_OLIGO_SEQ_PROBE_DESIGNER)

    def setup_cmd_parameters(self):
        return f"-c{os.path.abspath(CONFIG_OLIGO_SEQ_PROBE_DESIGNER)}"


class TestScrinshotProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self):
        with open(CONFIG_SCRINSHOT_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        return os.path.abspath(config["dir_output"])

    def setup_script(self):
        return os.path.abspath(SCRIPT_SCRINSHOT_PROBE_DESIGNER)

    def setup_cmd_parameters(self):
        return f"-c{os.path.abspath(CONFIG_SCRINSHOT_PROBE_DESIGNER)}"


class TestSeqfishplusProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self):
        with open(CONFIG_SEQFISHPLUS_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        return os.path.abspath(config["dir_output"])

    def setup_script(self):
        return os.path.abspath(SCRIPT_SEQFISHPLUS_PROBE_DESIGNER)

    def setup_cmd_parameters(self):
        return f"-c{os.path.abspath(CONFIG_SEQFISHPLUS_PROBE_DESIGNER)}"


class TestMerfishProbeDesigner(PipelinesBase, unittest.TestCase):
    def setup_output_dir(self):
        with open(CONFIG_MERFISH_PROBE_DESIGNER, "r") as handle:
            config = yaml.safe_load(handle)
        return os.path.abspath(config["dir_output"])

    def setup_script(self):
        return os.path.abspath(SCRIPT_MERFISH_PROBE_DESIGNER)

    def setup_cmd_parameters(self):
        return f"-c{os.path.abspath(CONFIG_MERFISH_PROBE_DESIGNER)}"
