import sys
import types

import pandas as pd
import pytest


try:
    import cel
except ModuleNotFoundError:
    cel = types.ModuleType("cel")

    class Context:
        def __init__(self):
            self.variables = {}
            self.functions = {}

        def add_variable(self, name, value):
            self.variables[name] = value

        def add_function(self, name, function):
            self.functions[name] = function

    cel.Context = Context
    cel.evaluate = lambda rule, context: False
    sys.modules["cel"] = cel

from abritamr import criteria

if not hasattr(criteria, "get_abritamr_configs"):
    criteria.get_abritamr_configs = lambda **kwargs: []

from abritamr import (
    amr_infer,
    amr_matrix,
    amr_report,
    catalog,
    cel_functions,
    drugclasses,
    filter_amrtype,
    filter_reportable,
    parse_amrtype,
    parse_finder,
    parse_gdstrules,
    parse_reportable,
    run_finder,
    run_sourmash,
    utils,
)
from abritamr.commands import (
    amr_status,
    infer,
    linelist,
    matrix,
    run,
    scan,
    update_database,
    utils_catalog,
    utils_rules,
)


def test_infer_thresholds_and_result_formats():
    results = pd.DataFrame(
        {
            "% Coverage of reference": [0.9, 0.89],
            "% Identity to reference": [0.91, 0.99],
        }
    )
    assert amr_infer.filter_results(results).index.tolist() == [0]

    entries = [
        {
            "sample_id": "sample-1",
            "species": "Species",
            "results": [
                {
                    "drugname": "Drug",
                    "inferred": "R",
                    "mechanisms": "gene",
                    "rule_id": "rule",
                    "rule_version": "1",
                    "source": "test",
                }
            ],
        }
    ]
    wide = amr_infer.gdst_results_to_df_wide(entries)
    long = amr_infer.gdst_results_to_df_long(entries)
    assert wide.loc[0, "Drug_gDST"] == "R"
    assert long.loc[0, "drugname"] == "Drug"
    assert amr_infer.priority_gdst() == {"S": 0, "I": 1, "R": 2}


def test_matrix_and_report_column_helpers(monkeypatch):
    monkeypatch.setattr(
        amr_matrix, "abritamr_matrix_columns", lambda **kwargs: ["Beta"]
    )
    matrix, columns = amr_matrix.wrangle_cols(
        pd.DataFrame({"abritamr_subclass": ["Beta"], "Element symbol": ["blaA"]}),
        {},
        ["Sample_id"],
    )
    assert matrix == {"Beta": "blaA"}
    assert columns == ["Sample_id", "Beta"]

    report, columns = amr_report.wrangle_cols(
        pd.DataFrame({"abritamr_subclass": ["Beta"], "Element symbol": ["blaA"]}),
        {},
        ["Sample_id"],
    )
    assert report == {"Beta": "blaA"}
    assert columns == ["Sample_id", "Beta"]


def test_catalog_mutation_formatting():
    assert catalog._capitalise("beta-lactam/carbapenem") == (
        "Beta-lactam/Carbapenemase"
    )
    assert catalog.sub_nt_mutations("A", 23, "G") == "c.[23A>G]"
    assert catalog.sub_nt_mutations("A", -10, "G", promoter=True) == "c.-10A>G"


def test_cel_helpers(monkeypatch):
    assert cel_functions.contains_any(["Beta-lactam", "Aminoglycoside"], "BETA")
    assert not cel_functions.contains_any(["Beta-lactam"], "tet")
    context = cel_functions.create_cel_context({"gene": "blaA"})
    monkeypatch.setattr(cel_functions, "evaluate", lambda rule, ctx: rule == "valid")
    assert cel_functions.evaluate_rule("valid", context)
    with pytest.raises(RuntimeError, match="Error evaluating rule"):
        monkeypatch.setattr(
            cel_functions,
            "evaluate",
            lambda rule, ctx: (_ for _ in ()).throw(ValueError("bad rule")),
        )
        cel_functions.evaluate_rule("invalid", context)


def test_criteria_validation():
    assert criteria.Criteria("id", "1", "true", status="reportable").status == (
        "reportable"
    )
    with pytest.raises(ValueError, match="drugname"):
        criteria.Criteria("id", "1", "true", inferred="R")
    with pytest.raises(ValueError, match="status"):
        criteria.Criteria("id", "1", "true")


def test_classification_and_reportability_lookups():
    refs = pd.DataFrame(
        {
            "refseq_protein_accession": ["WP_1"],
            "refseq_nucleotide_accession": [""],
            "genbank_protein_accession": [""],
            "genbank_nucleotide_accession": [""],
            "abritamr_accession_key": ["WP_1"],
            "abritamr_class": ["Beta-lactam"],
            "abritamr_subclass": ["Penicillin"],
        }
    )
    assert (
        drugclasses.get_class("abritamr_accession_key", "WP_1", refs, "abritamr_class")
        == "Beta-lactam"
    )
    assert (
        drugclasses.get_class(
            "abritamr_accession_key", "missing", refs, "abritamr_class"
        )
        == "NA"
    )
    assert parse_reportable.find_classes(refs, "WP_1") == (
        "Beta-lactam",
        "Penicillin",
    )


def test_amrtype_filter_helpers():
    expression = filter_amrtype.filter_string(
        {
            "gene": "['blaA', 'blaB']",
            "species": "'Species'",
            "exception": "None",
            "amrtype": "ignored",
        }
    )
    assert "'blaA' in gene" in expression
    assert "'blaB' in gene" in expression
    assert "species in 'Species'" in expression
    assert "exception" not in expression


def test_reportability_default_for_unmatched_accession():
    result = {
        "Closest reference accession": "-",
        "Element symbol": "blaA",
    }
    refs = pd.DataFrame({"abritamr_accession_key": ["WP_1"]})
    filtered = filter_reportable.construct_filter(result, refs)
    assert filtered["abritamr_priority_status"] == "not-reportable"
    assert filtered["criteria_id"] == ""


def test_parse_amrtype_applies_label(monkeypatch):
    monkeypatch.setattr(parse_amrtype, "construct_filter", lambda **kwargs: "ESBL")
    rows = [{"abritamr_subclass": "Beta-lactam", "Element symbol": "blaA"}]
    assert (
        parse_amrtype.get_amr_type(rows, species="Species")[0]["abritamr AMR type"]
        == "ESBL"
    )


def test_parse_finder_reads_tabular_records(tmp_path):
    result_file = tmp_path / "finder.tsv"
    result_file.write_text("sample_id\tElement symbol\nsample-1\tblaA\n")
    assert parse_finder.amrf2dict(str(result_file)) == [
        {"sample_id": "sample-1", "Element symbol": "blaA"}
    ]


def test_gdst_rule_parsing():
    row = {
        "protein accession": "-",
        "nucleotide accession": "NC_000001.1:10-20",
        "mutation": "-",
        "gene": "blaA",
    }
    assert parse_gdstrules.get_accession_key(row) == "NC_000001.1"
    assert "contains_any(row.abritamr_accession_key,'NC_000001.1')" in (
        parse_gdstrules.parse_rule(row, {})
    )
    assert parse_gdstrules.special_rule_for_amrrules_rna("c.[10A>G]") == "c.[10A>G]"
    assert (
        parse_gdstrules.special_rule_for_amrrules_rna("c.[10A>G]extra]") == "c.[10A>G]"
    )


def test_report_and_matrix_empty_results(monkeypatch):
    monkeypatch.setattr(amr_report, "abritamr_amrtype_columns", lambda: ["Sample_id"])
    assert amr_report.summary(pd.DataFrame()).columns.tolist() == ["Sample_id"]

    monkeypatch.setattr(
        amr_matrix, "abritamr_matrix_columns", lambda **kwargs: ["Beta"]
    )
    matrix_result = amr_matrix.summary(
        pd.DataFrame(
            {
                "% Coverage of reference": [100],
                "% Identity to reference": [100],
                "abritamr_subclass": ["Beta"],
                "Element symbol": ["blaA"],
            }
        ),
        facet="abritamr_subclass",
        sid="sample-1",
    )
    assert matrix_result.loc[0, "Beta"] == "blaA"


def test_amrfinder_command_and_output_parser(monkeypatch):
    monkeypatch.setattr(run_finder, "wrangle_species", lambda **kwargs: "-O Species")
    command = run_finder.generate_cmd(90, 80, "assembly.fa", 4, "Species")
    assert "amrfinder -n assembly.fa" in command
    assert "--threads 4 -O Species" in command
    assert run_finder.parse_output("Element symbol\tType\nblaA\tAMR\n") == [
        {"Element symbol": "blaA", "Type": "AMR"}
    ]


def test_sourmash_index_loader(monkeypatch):
    index = object()
    monkeypatch.setattr(run_sourmash.sourmash, "load_file_as_index", lambda path: index)
    assert run_sourmash.load_sourmash_index("reference.sbt.zip") is index


def test_utility_paths_and_columns(tmp_path):
    assert utils.check_path(str(tmp_path))
    assert not utils.check_path(str(tmp_path / "missing"))
    assert "sample_id" in utils.abritamr_scan_columns()
    assert "abritamr_priority_status" in utils.abritamr_status_columns()


def test_status_command_validates_and_projects_columns(monkeypatch):
    monkeypatch.setattr(
        amr_status,
        "generate_output",
        lambda amr, catalog: [{**amr[0], "abritamr_priority_status": "high"}],
    )
    monkeypatch.setattr(
        amr_status,
        "abritamr_status_columns",
        lambda: ["sample_id", "abritamr_priority_status"],
    )
    results = amr_status.do_typing(
        pd.DataFrame(
            {
                "sample_id": ["sample-1"],
                "species": ["Species"],
                "abritamr_subclass": ["Beta"],
            }
        ),
        "catalog.csv",
    )
    assert results.loc[0, "abritamr_priority_status"] == "high"
    with pytest.raises(SystemExit):
        amr_status.do_typing(pd.DataFrame({"wrong": [1]}), "catalog.csv")


def test_infer_and_linelist_command_wrappers(monkeypatch):
    inferred = [
        {
            "sample_id": "sample-1",
            "species": "Species",
            "results": [
                {
                    "drugname": "Drug",
                    "inferred": "R",
                    "mechanisms": "blaA",
                    "rule_id": "r1",
                    "rule_version": "1",
                    "source": "test",
                }
            ],
        }
    ]
    monkeypatch.setattr(infer, "gdst", lambda **kwargs: inferred)
    monkeypatch.setattr(
        infer, "gdst_results_to_df_long", amr_infer.gdst_results_to_df_long
    )
    amr = pd.DataFrame({"sample_id": ["sample-1"], "species": ["Species"]})
    assert infer.do_gdst(amr, "rules", "Susceptible", "long").loc[0, "gDST"] == "R"

    monkeypatch.setattr(
        linelist,
        "summary",
        lambda results, **kwargs: pd.DataFrame(
            {"Sample_id": [results.iloc[0]["sample_id"]]}
        ),
    )
    line_results = pd.DataFrame(
        {
            "sample_id": ["sample-1", "sample-2"],
            "species": ["Species", "Species"],
            "abritamr_subclass": ["Beta", "Beta"],
        }
    )
    assert linelist.generate_linelist(line_results, "csv", "compact", False, 90, 90)[
        "Sample_id"
    ].tolist() == ["sample-1", "sample-2"]


def test_matrix_and_run_command_helpers(monkeypatch, tmp_path):
    monkeypatch.setattr(
        matrix,
        "summary",
        lambda results, **kwargs: pd.DataFrame({"Sample_id": [kwargs["sid"]]}),
    )
    records = pd.DataFrame(
        {
            "sample_id": ["sample-1"],
            "species": ["Species"],
            "abritamr_subclass": ["Beta"],
        }
    )
    assert matrix.make_matrix(records, "abritamr_subclass", 90, 90, "catalog")[
        "Sample_id"
    ].tolist() == ["sample-1"]

    assert run.wrangle_outputs([{"sample_id": "sample-1", "gene": "blaA"}], ["gene"])[
        "gene"
    ].tolist() == ["blaA"]
    assert run.wrangle_outputs([], ["gene"]).empty
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    assert run.save_output(
        str(output_dir),
        "sample-1",
        pd.DataFrame({"gene": ["blaA"]}),
        "result",
        "tab",
    )
    assert (output_dir / "sample-1" / "result.txt").read_text().splitlines() == [
        "gene",
        "blaA",
    ]


def test_scan_input_builder_and_database_folder(tmp_path):
    assert scan.generate_inputs("assembly.fa", "sample-1", "contigs") == {
        "contigs": "assembly.fa",
        "sample_id": "sample-1",
    }
    assert update_database.create_db_folder(str(tmp_path / "db" / "nested"))
    assert (tmp_path / "db" / "nested").is_dir()


def test_catalog_and_rule_commands(monkeypatch, tmp_path):
    monkeypatch.setattr(utils_catalog, "create_db_folder", lambda path: True)
    updated_catalog = []
    monkeypatch.setattr(
        utils_catalog, "update_catalog", lambda args: updated_catalog.append(args)
    )
    args = types.SimpleNamespace(output_dir=str(tmp_path), catalog="catalog.csv")
    utils_catalog.catalog(args)
    assert updated_catalog == [args]

    monkeypatch.setattr(utils_rules, "create_db_folder", lambda path: True)
    monkeypatch.setattr(utils_rules, "check_path", lambda path: True)
    updated_rules = []
    monkeypatch.setattr(
        utils_rules, "update_rules", lambda args: updated_rules.append(args)
    )
    utils_rules.rules(args)
    assert updated_rules == [args]
