############################################
# imports
############################################

import unittest
from unittest.mock import patch

from oligo_designer_toolsuite._exceptions import (
    EmptyResultError,
    ExternalToolError,
    NetworkError,
    OligoDesignerError,
)
from oligo_designer_toolsuite.sequence_generator._ftp_loader import _connect
from oligo_designer_toolsuite.utils._checkers_and_helpers import run_external_tool

############################################
# Tests
############################################


class TestEmptyResultError(unittest.TestCase):
    """Tests the error that replaced ``sys.exit(1)`` when no oligos are left."""

    def test_it_exits_like_sys_exit_1_and_is_a_package_error(self):
        error = EmptyResultError("No oligos are left after filtering.")

        self.assertIsInstance(error, SystemExit)
        self.assertEqual(error.code, 1)
        self.assertIsInstance(error, OligoDesignerError)


class TestRunExternalTool(unittest.TestCase):
    """Tests the wrapper that reports why a command line tool failed."""

    def test_a_successful_tool_does_not_raise(self):
        run_external_tool(["true"])

    def test_the_error_names_the_tool_and_says_why(self):
        with self.assertRaises(ExternalToolError) as raised:
            run_external_tool(["sh", "-c", "echo 'BLAST Database error: No alias found' >&2; exit 2"])

        self.assertEqual(
            str(raised.exception), "sh failed with exit status 2: BLAST Database error: No alias found"
        )


class TestConnect(unittest.TestCase):
    """Tests the FTP connection that the NCBI and Ensembl loaders share."""

    def test_an_unreachable_server_raises_a_network_error(self):
        with patch(
            "oligo_designer_toolsuite.sequence_generator._ftp_loader.FTP", side_effect=OSError("timed out")
        ):
            with self.assertRaises(NetworkError) as raised:
                _connect("ftp.ncbi.nlm.nih.gov")

        self.assertEqual(
            str(raised.exception), "Could not connect to ftp.ncbi.nlm.nih.gov. Please try again later."
        )


if __name__ == "__main__":
    unittest.main()
