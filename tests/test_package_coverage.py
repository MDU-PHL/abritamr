import importlib
import logging
import re
import sys
from types import SimpleNamespace
from unittest.mock import Mock
from xml.etree import ElementTree
from zipfile import ZipFile

import numpy
import pandas
import pytest

from abritamr.Collate import Collate, MduCollate
from abritamr.CustomLog import CustomFormatter
from abritamr.RunFinder import RunFinder


@pytest.fixture
def update_module(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    return importlib.import_module("abritamr.Update")


@pytest.fixture
def cli_module(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    return importlib.import_module("abritamr.abritamr")


def bare_run_finder():
    finder = RunFinder.__new__(RunFinder)
    finder.logger = logging.getLogger(__name__)
    finder.db = "2026-03-24.1"
    finder.amrfinder_db = finder.db
    finder.organism = ""
    finder.input = "assembly.fa"
    finder.run_type = "assembly"
    finder.jobs = 4
    finder.prefix = "sample"
    finder.identity = ""
    return finder


def bare_mdu_collate():
    collate = MduCollate.__new__(MduCollate)
    collate.logger = logging.getLogger(__name__)
    collate.NONE_CODES = {
        "Salmonella": "CPase_ESBL_AmpC_16S_NEG",
        "Shigella": "CPase_ESBL_AmpC_16S_NEG",
        "Staphylococcus": "Mec_VanAB_Linez_NEG",
        "Enterococcus": "Van_Linez_NEG",
        "Other": "CPase_16S_mcr_NEG",
    }
    return collate


def test_custom_formatter_uses_level_color():
    record = logging.LogRecord("test", logging.WARNING, "", 1, "warning", (), None)

    assert CustomFormatter().format(record).startswith("\x1b[33;21m[WARNING:")
    assert CustomFormatter().format(record).endswith("\x1b[0m")


@pytest.mark.parametrize(
    ("database", "stderr", "returncode", "expected"),
    [
        ("2026-03-24.1", "", 0, True),
        ("other-version", "", 0, False),
        ("", "AMRFinderPlus 2026-03-24\n", 0, True),
        ("", "version unavailable\n", 0, False),
    ],
)
def test_check_amrfinder(database, stderr, returncode, expected, monkeypatch):
    finder = bare_run_finder()
    finder.amrfinder_db = database
    monkeypatch.setattr(
        "abritamr.RunFinder.subprocess.run",
        Mock(return_value=SimpleNamespace(stderr=stderr, returncode=returncode)),
    )

    assert finder._check_amrfinder() is expected


@pytest.mark.parametrize(("returncode", "expected"), [(0, True), (1, None)])
def test_run_cmd_reports_subprocess_result(returncode, expected, monkeypatch):
    finder = bare_run_finder()
    monkeypatch.setattr(
        "abritamr.RunFinder.subprocess.run",
        Mock(return_value=SimpleNamespace(returncode=returncode, stderr="failure")),
    )

    assert finder._run_cmd("amrfinder") is expected


def test_run_returns_run_data_after_validating_output(monkeypatch):
    finder = bare_run_finder()
    finder._check_amrfinder = Mock(return_value=True)
    finder._run_cmd = Mock()
    finder._check_outputs = Mock()

    result = finder.run()

    assert result == ("assembly", "assembly.fa", "sample")
    finder._run_cmd.assert_called_once_with(finder._generate_cmd())
    finder._check_outputs.assert_called_once_with()


def test_collate_joins_other_than_isolate_and_deduplicates():
    collate = Collate.__new__(Collate)
    values = {"Isolate": "sample", "Beta-lactam": ["blaA", "blaA", "blaB"]}

    assert collate.joins(values) == {"Isolate": "sample", "Beta-lactam": "blaA,blaB"}


def test_collate_adds_caret_to_nonempty_partial_values():
    collate = Collate.__new__(Collate)
    frame = pandas.DataFrame({"Isolate": ["sample"], "ESBL": ["blaCTX*"], "Other": [""]})

    result = collate._add_caret(frame, ["ESBL", "Other"])

    assert result.loc[0, "ESBL"] == "blaCTX^"
    assert result.loc[0, "Other"] == ""


def test_collate_batches_and_combines_results(tmp_path, monkeypatch):
    collate = Collate.__new__(Collate)
    collate.logger = logging.getLogger(__name__)
    batch = tmp_path / "batch.tsv"
    batch.write_text("sample-a\tassembly-a.fa\nsample-b\tassembly-b.fa\n")
    match = pandas.DataFrame({"Isolate": ["one"], "Beta-lactam": ["blaA"]})
    partial = pandas.DataFrame({"Isolate": ["one"]})
    virulence = pandas.DataFrame({"Isolate": ["one"], "Other": ["gene"]})
    calls = []

    def collate_sample(prefix):
        calls.append(prefix)
        return match, partial, virulence

    monkeypatch.setattr(collate, "collate", collate_sample)
    results = collate._batch_collate(batch)

    assert calls == ["sample-a", "sample-b"]
    assert [len(result) for result in results] == [2, 2, 2]
    assert collate._combine_df(pandas.DataFrame(), match).equals(match)


def test_mdu_qc_validates_columns_and_adds_sentinel(tmp_path):
    collate = bare_mdu_collate()
    qc = tmp_path / "qc.csv"
    qc.write_text("sample,SPECIES_EXP,SPECIES_OBS,TEST_QC\ns1,E. coli,E. coli,PASS\n")
    collate.mduqc = qc

    result = collate.mdu_qc_tab()

    assert list(result.columns) == ["ISOLATE", "SPECIES_EXP", "SPECIES_OBS", "TEST_QC"]
    assert result.iloc[-1]["ISOLATE"] == "9999-99888"

    qc.write_text("sample,SPECIES_EXP\ns1,E. coli\n")
    with pytest.raises(SystemExit):
        collate.mdu_qc_tab()


@pytest.mark.parametrize(
    ("gene", "expected"),
    [("blaCTX-M-15", "CTX-M-15"), ("blaCTX-M-15*", "CTX-M-15*"), ("blaZ", "blaZ")],
)
def test_mdu_strip_bla(gene, expected):
    assert bare_mdu_collate().strip_bla(gene) == expected


def test_mdu_ids_and_negative_codes():
    collate = bare_mdu_collate()
    pattern = re.compile(r"(?P<id>[0-9]{4}-[0-9]{5,6})-?(?P<itemcode>.{1,})?")

    assert collate.assign_mduid("1234-56789-ABC", pattern) == "1234-56789"
    assert collate.assign_itemcode("1234-56789-ABC", pattern) == "ABC"
    assert collate.assign_mduid("sample/path", pattern) == "path"
    assert collate.none_replacement_code("Salmonella") == "CPase_ESBL_AmpC_16S_NEG"
    assert collate.none_replacement_code("Unknown") == "CPase_16S_mcr_NEG"


@pytest.mark.parametrize(
    ("method", "matches", "nonmatches"),
    [
        ("_ampicillin_res_sal", ["Beta-lactam", "Ampicillin"], ["Other"]),
        ("_chloramphenicol_res_sal", ["Chloramphenicol"], ["Tetracycline"]),
        ("_cefo_esbl_res_sal", ["ESBL"], ["AmpC"]),
        ("_cefo_ampc_res_sal", ["AmpC"], ["ESBL"]),
        ("_carbapenem_res_salmo", ["Carbapenemase"], ["KPC"]),
        ("_azi_res_salmo", ["Azithromycin", "Macrolide"], ["Other"]),
        ("_gentamicin_res_salm", ["Gentamicin"], ["Kanamycin"]),
        ("_kanamycin_res_salm", ["Kanamycin"], ["Gentamicin"]),
        ("_streptomycin_res_salm", ["Streptomycin"], ["Tetracycline"]),
        ("_spectinomycin_res_salm", ["Streptomycin"], ["Tetracycline"]),
        ("_tetra_res_salmo", ["Tetracycline"], ["Other"]),
        ("_cipro_res_salmo", ["Quinolone"], ["Other"]),
        ("_sulf_res_salmo", ["Sulfonamide"], ["Other"]),
        ("_trimet_res_salmo", ["Trimethoprim"], ["Other"]),
        ("_rmt_res_salmo", ["Aminoglycosides (Ribosomal methyltransferase)"], ["Other"]),
        ("_colistin_res_salmo", ["Colistin"], ["Other"]),
    ],
)
def test_salmonella_drug_column_filters(method, matches, nonmatches):
    collate = bare_mdu_collate()
    filter_column = getattr(collate, method)

    assert all(filter_column(column, "gene") == "gene" for column in matches)
    assert all(filter_column(column, "gene") == "" for column in nonmatches)


def test_reporting_logic_general_applies_species_exclusions():
    collate = bare_mdu_collate()
    frame = pandas.DataFrame(
        {
            "Isolate": ["sample"],
            "Carbapenemase (MBL)": ["blaL1"],
            "ESBL": ["blaEC-1"],
            "Vancomycin": ["vanA,other"],
            "Methicillin": ["mecA,other"],
        }
    )
    row = next(frame.iterrows())

    reported, not_reported = collate.reporting_logic_general(
        row, "Stenotrophomonas maltophilia"
    )

    assert reported == ["vanA", "mecA"]
    assert not_reported == ["blaL1", "blaEC-1", "other", "other"]


def test_mdu_collects_genes_and_extracts_qc_passed_isolates(tmp_path):
    collate = bare_mdu_collate()
    row = (
        0,
        pandas.Series({"Isolate": "sample", "ESBL": "blaA,blaB", "Other": numpy.nan}),
    )
    assert collate.get_all_genes(row) == ["sample", "blaA", "blaB"]

    qc = tmp_path / "qc.csv"
    qc.write_text(
        "ISOLATE,SPECIES_EXP,SPECIES_OBS,TEST_QC\n"
        "pass,E. coli,E. coli,PASS\n"
        "fail,E. coli,E. coli,FAIL\n"
    )
    collate.mduqc = qc
    assert collate._extract_plus_isolates("E. coli") == ["pass"]


def test_mdu_reporting_logic_salmonella_classifies_resistance():
    collate = bare_mdu_collate()
    frame = pandas.DataFrame(
        {
            "Isolate": ["1234-56789-SALM"],
            "Ampicillin": ["blaTEM"],
            "Sulfonamide": ["sul1"],
            "Trimethoprim": ["dfrA"],
            "Ciprofloxacin": [""],
        }
    )

    result = collate.reporting_logic_salmonella(next(frame.iterrows()))

    assert result["MDU Sample ID"] == "1234-56789"
    assert result["Item code"] == "SALM"
    assert result["Ampicillin - ResMech"] == "blaTEM"
    assert result["Ampicillin - Interpretation"] == "Resistant"
    assert set(result["Trim-Sulpha - ResMech"].split(";")) == {"dfrA", "sul1"}
    assert result["Ciprofloxacin - Interpretation"] == "Susceptible"


def test_mdu_reporting_salmonella_selects_requested_isolates(tmp_path):
    collate = bare_mdu_collate()
    match = tmp_path / "matches.tsv"
    pandas.DataFrame({"Isolate": ["keep", "skip"], "Ampicillin": ["blaA", "blaB"]}).to_csv(
        match, sep="\t", index=False
    )
    collate.reporting_logic_salmonella = Mock(
        side_effect=lambda row: {
            "MDU Sample ID": row[1]["Isolate"],
            "Item code": "",
            **{
                f"{name} - {suffix}": "None detected"
                for name in [
                    "Ampicillin",
                    "Cefotaxime (ESBL)",
                    "Cefotaxime (AmpC)",
                    "Tetracycline",
                    "Gentamicin",
                    "Kanamycin",
                    "Streptomycin",
                    "Sulfathiazole",
                    "Trimethoprim",
                    "Trim-Sulpha",
                    "Chloramphenicol",
                    "Ciprofloxacin",
                    "Meropenem",
                    "Azithromycin",
                    "Aminoglycosides (RMT)",
                    "Colistin",
                    "Other",
                ]
                for suffix in ["ResMech", "Interpretation"]
            },
        }
    )

    result = collate.mdu_reporting_salmonella(match, ["keep"])

    assert list(result["MDU Sample ID"]) == ["keep"]


def test_mdu_spreadsheet_saves_general_and_interpreted_results(tmp_path, monkeypatch):
    collate = bare_mdu_collate()
    collate.sop_name = "report"
    collate.runid = "RUN"
    monkeypatch.chdir(tmp_path)
    matches = pandas.DataFrame({"value": ["match"]})
    partials = pandas.DataFrame({"value": ["partial"]})

    collate.save_spreadsheet_general(matches, partials)

    output = tmp_path / "RUN_report.xlsx"
    assert output.exists()
    with ZipFile(output) as workbook:
        contents = ElementTree.fromstring(workbook.read("xl/workbook.xml"))
    namespace = {"xlsx": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    assert [sheet.attrib["name"] for sheet in contents.findall(".//xlsx:sheet", namespace)] == [
        "report",
        "Passed QC partial",
    ]

    collate.save_spreadsheet_interpreted([("Salmonella enterica", matches)])
    with ZipFile(output) as workbook:
        contents = ElementTree.fromstring(workbook.read("xl/workbook.xml"))
    assert [sheet.attrib["name"] for sheet in contents.findall(".//xlsx:sheet", namespace)] == [
        "report-01"
    ]


def test_update_transforms_and_config(update_module):
    update = update_module
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", update._get_date())
    rename, other_amr, other_non_amr, oxa, address = update._get_vars()
    assert rename["FLUOROQUINOLONE"] == "Quinolone"
    assert "BACITRACIN" in other_amr
    assert "ARSENIC" in other_non_amr
    assert "optrA" in oxa
    assert "@" in address

    assert update._capitalise("CARBAPENEM/OTHER") == "Carbapenemase/Other"
    assert update._oxa_phen({"subclass": "FLORFENICOL"}) == (
        "Oxazolidinone/Phenicol",
        "Florfenicol",
    )
    assert update._other_antimicrobials({"subclass": "BACITRACIN"}) == (
        "Other antimicrobial",
        "Bacitracin",
    )
    assert update._other_non_antimicrobials({"subclass": "ARSENIC"}) == (
        "Other non-antimicrobial",
        "Arsenic",
    )
    assert update.virulence({"class": "", "subclass": ""}) == ("Virulence", "Other")


def test_update_classification_helpers(update_module):
    update = update_module
    assert update._get_keys(
        {
            "rename_key": {"A": "B"},
            "other_amr": ["AMR"],
            "other_non_amr": ["NON-AMR"],
            "oxa_phen_list": ["gene"],
            "email_address": "review@example.org",
        }
    ) == ({"A": "B"}, ["AMR"], ["NON-AMR"], ["gene"], "review@example.org")

    row = {"gene_family": "cfr", "class": "", "subclass": "OXAZOLIDINONE"}
    assert update.cfr(row) == ("Multidrug", "Oxazolidinone")
    aminoglycoside = {
        "product_name": "16S rRNA methyltransferase",
        "class": "AMINOGLYCOSIDE",
        "subclass": "AMINOGLYCOSIDE",
    }
    assert update._aminoglycosides(aminoglycoside) == (
        "Aminoglycoside",
        "Aminoglycosides (Ribosomal methyltransferase)",
    )
    assert update._rename(
        {"FLUOROQUINOLONE": "Quinolone", "CARBAPENEM": "Carbapenemase"},
        {"class": "FLUOROQUINOLONE", "subclass": "CARBAPENEM"},
    ) == ("Quinolone", "Carbapenemase")
    assert update.virulence({"class": "INTIMIN", "subclass": "ECOLI"}) == (
        "Virulence",
        "Intimin_ecoli",
    )
    assert update.virulence({"class": "STX2", "subclass": "STX"}) == ("Virulence", "Stx")


def test_update_existing_catalog_helpers(update_module, monkeypatch):
    update = update_module
    expected_path = update.pathlib.Path(update.__file__).parent / "db" / "refgenes_latest.csv"
    assert update._check_existing() == expected_path

    monkeypatch.setattr(update, "_check_existing", lambda: False)
    assert update._get_previous_refgenes() is False

    previous = pandas.DataFrame(
        [{"key": "existing", "class_new": "AMR", "subclass_new": "ESBL"}]
    )
    monkeypatch.setattr(update, "_check_existing", lambda: expected_path)
    monkeypatch.setattr(update.pandas, "read_csv", lambda path: previous.copy())
    assert update._get_previous_refgenes().equals(previous)


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (
            {
                "subtype": "AMR-SUSCEPTIBLE",
                "product_name": "",
                "subclass": "",
                "allele": "",
                "gene_family": "",
                "class": "BETA-LACTAM",
            },
            ("Beta-lactam", "Beta-lactam (Susceptible)"),
        ),
        (
            {
                "subtype": "AMR",
                "product_name": "carbapenem-hydrolyzing enzyme",
                "subclass": "CARBAPENEM",
                "allele": "",
                "gene_family": "",
                "class": "BETA-LACTAM",
            },
            ("Beta-lactam", "Carbapenemase"),
        ),
        (
            {
                "subtype": "AMR",
                "product_name": "",
                "subclass": "PENICILLIN",
                "allele": "",
                "gene_family": "blaZ",
                "class": "BETA-LACTAM",
            },
            ("Beta-lactam", "Penicillin resistance (Staphylococcus aureus)"),
        ),
    ],
)
def test_update_beta_lactam_classification(update_module, row, expected):
    assert update_module._beta_lactams(row) == expected


def test_update_makes_keys_and_classifies_rows(update_module):
    update = update_module
    source = pandas.DataFrame(
        [
            {
                "allele": "gyrA",
                "whitelisted_taxa": "Escherichia",
                "refseq_nucleotide_accession": "NC_1",
                "refseq_protein_accession": "",
                "genbank_protein_accession": "WP_1",
                "subtype": "POINT",
            },
            {
                "allele": "blaZ",
                "whitelisted_taxa": "",
                "refseq_nucleotide_accession": "NC_2",
                "refseq_protein_accession": "",
                "genbank_protein_accession": "WP_2",
                "subtype": "AMR",
            },
        ]
    )
    keyed = update._make_key(source)
    assert keyed["key"].tolist() == ["gyrA_Escherichia_NC_1", "WP_2"]

    classified = update._make_dict(
        pandas.DataFrame(
            [
                {
                    "class": "BETA-LACTAM",
                    "subclass": "CARBAPENEM",
                    "subtype": "AMR",
                    "product_name": "carbapenem-hydrolyzing enzyme",
                    "allele": "blaKPC",
                    "gene_family": "blaKPC",
                    "type": "AMR",
                },
                {
                    "class": "",
                    "subclass": "",
                    "subtype": "",
                    "product_name": "multidrug efflux",
                    "allele": "",
                    "gene_family": "",
                    "type": "AMR",
                },
            ]
        ),
        other_amr=[],
        other_non_amr=[],
        rename_key={},
    )
    assert classified[0]["enhanced_subclass"] == "Carbapenemase"
    assert (classified[1]["enhanced_class"], classified[1]["enhanced_subclass"]) == (
        "Multidrug",
        "Other",
    )


def test_update_catalog_download_and_archive_are_mocked(update_module, monkeypatch, tmp_path):
    update = update_module
    run = Mock(return_value=SimpleNamespace(returncode=0, stderr=""))
    monkeypatch.setattr(update.subprocess, "run", run)
    assert update._get_new_catalog()
    assert run.call_args.args[0].startswith("wget -O ReferenceGeneCatalog.txt ")
    assert run.call_args.kwargs["shell"] is True

    monkeypatch.setattr(update, "_get_date", lambda: "2026-01-02")
    update._archive_old_ref_catalog()
    assert run.call_args.args[0] == [
        "cp",
        str(update.pathlib.Path(update.__file__).parent / "db" / "refgenes_latest.csv"),
        str(update.pathlib.Path(update.__file__).parent / "db" / "refgenes_latest.csv.2026-01-02"),
    ]
    monkeypatch.setattr(
        update.subprocess,
        "run",
        Mock(return_value=SimpleNamespace(returncode=1, stderr="network error")),
    )
    with pytest.raises(SystemExit):
        update._get_new_catalog()


def test_update_open_catalog_fills_missing_values_and_builds_key(update_module, monkeypatch):
    update = update_module
    monkeypatch.setattr(update, "_get_new_catalog", Mock())
    monkeypatch.setattr(
        update.pandas,
        "read_csv",
        lambda *args, **kwargs: pandas.DataFrame(
            [
                {
                    "allele": "gene",
                    "whitelisted_taxa": numpy.nan,
                    "refseq_nucleotide_accession": "NC_1",
                    "refseq_protein_accession": "",
                    "genbank_protein_accession": "WP_1",
                    "subtype": "POINT",
                }
            ]
        ),
    )

    result = update._open_catalog()

    assert result.loc[0, "whitelisted_taxa"] == ""
    assert result.loc[0, "key"] == "gene__NC_1"


def test_update_entry_comparison_marks_new_and_changed_rows(update_module):
    update = update_module
    new = pandas.DataFrame(
        [
            {"key": "same", "class": "BETA-LACTAM", "subclass": "ESBL", "Status": ""},
            {"key": "new", "class": "AMR", "subclass": "Other", "Status": ""},
        ]
    )
    old = pandas.DataFrame(
        [
            {"key": "same", "class_new": "BETA-LACTAM", "subclass_new": "AmpC"},
        ]
    )

    result = update._new_entries(new, old).set_index("key")

    assert result.loc["same", "Status"] == "updated"
    assert result.loc["same", "Previous_subclass"] == "AmpC"
    assert result.loc["new", "Status"] == "new"
    assert update._compare_to_existing(new.to_dict("records"), False).equals(
        pandas.DataFrame(new.to_dict("records"))
    )


def test_update_save_and_email_use_expected_outputs(update_module, monkeypatch, tmp_path):
    update = update_module
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(update, "_get_date", lambda: "2026-01-02")
    frame = pandas.DataFrame({"key": ["gene"]})
    saved = update._save_df(frame)
    assert saved == "refgenes_2026-01-02.csv"
    assert pandas.read_csv(saved).loc[0, "key"] == "gene"

    run = Mock(return_value=SimpleNamespace(returncode=0, stderr=""))
    monkeypatch.setattr(update.subprocess, "run", run)
    update._email("review@example.org", saved)
    assert "review@example.org" in run.call_args.args[0]
    assert saved in run.call_args.args[0]


def test_update_create_refgenes_orchestrates_without_network(update_module, monkeypatch):
    update = update_module
    monkeypatch.setattr(update, "_get_vars", Mock(return_value=({}, [], [], [], "review@example.org")))
    monkeypatch.setattr(update, "_open_catalog", Mock(return_value=pandas.DataFrame()))
    monkeypatch.setattr(update, "_make_dict", Mock(return_value=[]))
    monkeypatch.setattr(update, "_get_previous_refgenes", Mock(return_value=False))
    monkeypatch.setattr(update, "_compare_to_existing", Mock(return_value=[]))
    monkeypatch.setattr(update, "_save_df", Mock(return_value="output.csv"))
    monkeypatch.setattr(update, "_email", Mock())

    update.create_refgenes()

    update._email.assert_called_once_with(adrs="review@example.org", pth="output.csv")


def test_cli_subcommands_dispatch_to_pipeline_functions(cli_module, monkeypatch):
    cli = cli_module
    for command, function in (
        (["run", "--contigs", "assembly.fa"], "run_pipeline"),
        (["report", "--qc", "qc.csv"], "mdu"),
        (["update_db"], "update_db"),
    ):
        handler = Mock()
        monkeypatch.setattr(cli, function, handler)
        monkeypatch.setattr(sys, "argv", ["abritamr", *command])
        cli.main()
        handler.assert_called_once()


def test_cli_help_and_version(cli_module, monkeypatch, capsys):
    cli = cli_module
    monkeypatch.setattr(sys, "argv", ["abritamr"])
    cli.main()
    assert "usage:" in capsys.readouterr().err

    monkeypatch.setattr(sys, "argv", ["abritamr", "--version"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 0
    assert cli.__version__ in capsys.readouterr().out


def test_pipeline_entrypoints_orchestrate_components(cli_module, monkeypatch):
    cli = cli_module
    setup = Mock()
    setup.return_value.setup.return_value = "input"
    finder = Mock()
    finder.return_value.run.return_value = "amr"
    collate = Mock()
    monkeypatch.setattr(cli, "SetupAMR", setup)
    monkeypatch.setattr(cli, "RunFinder", finder)
    monkeypatch.setattr(cli, "Collate", collate)

    cli.run_pipeline(SimpleNamespace())

    setup.assert_called_once()
    finder.assert_called_once_with("input")
    collate.assert_called_once_with("amr")
    collate.return_value.run.assert_called_once_with()


def test_report_and_update_entrypoints_orchestrate_components(cli_module, monkeypatch):
    cli = cli_module
    setup = Mock()
    setup.return_value.setup.return_value = "report inputs"
    report = Mock()
    monkeypatch.setattr(cli, "SetupMDU", setup)
    monkeypatch.setattr(cli, "MduCollate", report)

    cli.mdu(SimpleNamespace())

    setup.assert_called_once()
    report.assert_called_once_with("report inputs")
    report.return_value.run.assert_called_once_with()

    monkeypatch.setattr(cli, "create_refgenes", Mock())
    cli.update_db(SimpleNamespace())
    cli.create_refgenes.assert_called_once_with()
