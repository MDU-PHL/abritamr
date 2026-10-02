"""Check thay all the depednencies for abritamr are installed"""

from abritamr.logger import log
from abritamr.utils import check_amrfinder, check_sourmash, check_any2fasta


def check_deps() -> bool:

    check_amrfinder()
    check_sourmash()
    check_any2fasta()

    return True
