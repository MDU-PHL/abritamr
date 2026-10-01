"""Configure the package-wide logger."""

import logging

logging.basicConfig(
    format="[%(levelname)s:%(asctime)s] %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO,
)

log = logging.getLogger(__name__)
