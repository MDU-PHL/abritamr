from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class Criteria:
    drugclass: Optional[str]
    drugsubclass: Optional[str]
    species: Optional[str]
    action: str   # "report" or "not_report"
    drug: Optional[str]
    abx: Optional[str]
    gene: Optional[str] = None


def default_criteria()


# abritamr = Abritamr(RULES)
# abritamr.run(data, output, etc)