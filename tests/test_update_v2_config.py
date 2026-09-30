import json
import pathlib

import pytest


CONFIG_PATH = (
    pathlib.Path(__file__).parent.parent
    / "abritamr"
    / "db"
    / "update_vars_v2.json"
)


@pytest.fixture
def config():
    with CONFIG_PATH.open() as config_file:
        return json.load(config_file)


def test_v2_config_has_expected_sections(config):
    assert set(config) == {
        "for_curation",
        "other_antimicrobials",
        "other_non_amr",
    }
    assert set(config["for_curation"]["complex"]) == {
        "columns",
        "class",
        "type",
        "rules",
    }


def test_v2_complex_rules_have_valid_conditions(config):
    complex_rules = config["for_curation"]["complex"]
    supported_operators = {"eq", "in", "not in", "not eq", "r in"}

    for rule_group in complex_rules["rules"].values():
        for rule in rule_group:
            assert rule["rule_join"] in {"and", "or"}
            assert rule["class"]
            assert rule["subclass"]
            assert rule["rule"]

            for condition in rule["rule"]:
                assert set(condition) == {"column", "value", "operator"}
                assert condition["column"] in {
                    "allele",
                    "class",
                    "gene_family",
                    "product_name",
                    "subclass",
                    "subtype",
                }
                assert condition["operator"] in supported_operators


@pytest.mark.parametrize(
    ("section", "expected_class"),
    [
        ("other_antimicrobials", "Other antimicrobial"),
        ("other_non_amr", "Other non-antimicrobial"),
    ],
)
def test_v2_simple_rules_have_class_and_subclass(section, expected_class, config):
    classification = config[section]

    assert classification["class"]
    assert classification["rules"] == {
        "class": expected_class,
        "subclass": "SUBCLASS",
    }


def test_v2_complex_rule_groups_cover_supported_classes(config):
    complex_rules = config["for_curation"]["complex"]

    assert set(complex_rules["rules"]) == {
        "AMINOGLYCOSIDE",
        "BETA-LACTAM",
        "VIRULENCE",
    }
    assert set(complex_rules["class"]) == {"AMINOGLYCOSIDE", "BETA-LACTAM"}
    assert complex_rules["type"] == ["VIRULENCE"]


def test_v2_beta_lactam_rules_include_required_subclasses(config):
    rules = config["for_curation"]["complex"]["rules"]["BETA-LACTAM"]

    assert {rule["subclass"] for rule in rules} == {
        "Beta-lactam (Susceptible)",
        "Carbapenemase",
        "Carbapenemase (MBL)",
        "Carbapenemase (OXA-51 family)",
        "Penicillin resistance (Staphylococcus aureus)",
        "AmpC",
        "ESBL (KPC variant)",
        "ESBL",
    }
