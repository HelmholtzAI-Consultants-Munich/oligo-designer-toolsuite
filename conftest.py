import pytest
import warnings
import os
import shutil
import glob

@pytest.fixture(autouse=True, scope='session')
def suppress_effidict_warnings():
    # Suppress specific warnings related to the FileNotFoundError in effidict
    warnings.filterwarnings("ignore", category=ResourceWarning, message=r".*effidict.*")

    yield

    # Additional code to handle the cleanup more gracefully if necessary
    def handle_effidict_cleanup():
        base_path = '/home/runner/work/oligo-designer-toolsuite/oligo-designer-toolsuite/tmp_oligo_database/test_oligo_database/'
        cache_files_pattern = os.path.join(base_path, 'cache_files*')

        for cache_dir in glob.glob(cache_files_pattern):
            try:
                if os.path.exists(cache_dir):
                    shutil.rmtree(cache_dir)
            except FileNotFoundError:
                pass

    handle_effidict_cleanup()
