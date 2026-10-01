"""Data models and loaders for AMR classification and inference criteria."""

from dataclasses import dataclass
import pandas as pd


@dataclass(frozen=True)
class InferRules:
    """Represent a rule that infers susceptibility for a drug and species."""

    drugname: str
    species: str
    amrrules_species: str
    rule_id: str
    rule_version: str
    rule: str
    inferred: str
    original_rule: str = None
    source: str = "not supplied"


@dataclass(frozen=True)
class ClassDefintion:
    """Represent a rule assigning AMR classes to reference catalog entries."""

    class_id: int
    class_curation_id: str
    definition: str
    abritamr_class: str = None
    abritamr_subclass: str = None


@dataclass(frozen=True)
class Criteria:
    """Represent configured reporting, typing, or susceptibility criteria."""

    criteria_id: str
    criteria_version: str
    criteria: str
    additional_status_criteria: str = None
    additional_type_criteria: str = None
    status: str = None
    amrtype: str = None
    species: str = None
    drugname: str = None
    inferred: str = None
    translation: str = None

    def __post_init__(self):
        """Require a drug name for inference and either status or inference."""
        if not self.drugname and self.inferred:
            raise ValueError(f"You must supply a 'drugname'")
        if not self.status and not self.inferred:
            raise ValueError(
                f"You must supply one of status (reportable or not-reportable) or inferred (example resistant, R, susceptible)"
            )


def get_abritamr_reporting(cfgpath: str = "") -> list:
    """Load reporting and typing criteria from a CSV configuration file."""
    abritamr_rep = []
    try:
        reps = pd.read_csv(f"{cfgpath}")
    except Exception as e:
        raise RuntimeError(f"Error reading abritamr reporting rules: {e}")

    for criteria in reps.to_dict(orient="records"):
        # # print(criteria)
        rep = Criteria(**criteria)
        abritamr_rep.append(rep)

    return abritamr_rep


def get_abritamr_defs(cfgpath: str = "") -> list:
    """Load class definitions from a CSV configuration file."""
    abritamr_defs = []
    try:
        defs = pd.read_csv(f"{cfgpath}")
    except Exception as e:
        raise RuntimeError(f"Error reading abritamr class definitions: {e}")

    for d in defs.to_dict(orient="records"):
        dfs = ClassDefintion(**d)
        abritamr_defs.append(dfs)

    return abritamr_defs
