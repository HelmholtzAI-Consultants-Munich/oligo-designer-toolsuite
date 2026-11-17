############################################
# imports
############################################

import logging
import os
from datetime import datetime
from pathlib import Path

############################################
# Logging utils
############################################

# Library-wide logger
logger = logging.getLogger("oligo_designer_toolsuite")


def configure_root_logger(
    dir_output: str,
    pipeline_name: str,
    log_level: int = logging.NOTSET,
    include_console: bool = False,
) -> None:
    """
    Set up logging configuration for a pipeline.

    This function creates a consistent logging setup across all pipelines, creating a log file
    in the output directory and optionally writing to the console.

    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    :param pipeline_name: Name of the pipeline (used in log file name).
    :type pipeline_name: str
    :param log_level: Logging level (default: logging.NOTSET).
    :type log_level: int
    :param include_console: Whether to also log to console (default: False).
    :type include_console: bool
    """
    timestamp = datetime.now()

    # ensure output directory exists
    dir_output = os.path.abspath(dir_output)
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    file_logger = os.path.join(
        dir_output,
        f"log_{pipeline_name}_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
    )

    handlers: list[logging.Handler] = [logging.FileHandler(file_logger)]
    if include_console:
        handlers.append(logging.StreamHandler())

    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=log_level,
        handlers=handlers,
    )
