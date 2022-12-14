"""
Retrieve and extract sector data bundle.
"""

import logging

logger = logging.getLogger(__name__)

import os
import sys
import tarfile
from pathlib import Path

from _helpers import progress_retrieve, configure_logging


if __name__ == "__main__":
    configure_logging(snakemake)

    url = "https://zenodo.org/record/6412255/files/pypsa-eur-sec-data-bundle.tar.gz"

    # Should be [.../]pypsa-eur-sec/data:
    to_fn = Path(snakemake.output[0]).parent.parent

    tarball_fn = to_fn.joinpath("sector-bundle.tar.gz")

    logger.info(f"Downloading databundle from '{url}'.")
    progress_retrieve(url, tarball_fn)

    logger.info(f"Extracting databundle.")
    tarfile.open(tarball_fn).extractall(to_fn)

    tarball_fn.unlink()

    logger.info(f"Databundle available in '{to_fn}'.")
