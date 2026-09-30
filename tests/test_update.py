import json
import pathlib
import datetime
from types import SimpleNamespace

import pandas
import pytest

from abritamr import Update


def test_get_date_uses_iso_format():
    assert Update._get_date() == datetime.date.today().strftime("%Y-%m-%d")


def _row(**overrides):
    row = {
        "class": "OTHER",
        "subclass": "SUBCLASS",
        "gene_family": "gene",
        "product_name": "product",
        "subtype": "AMR",
        "type": "AMR",
        "allele": "gene",
    }
    row.update(overrides)
    return row


def test_get_keys_returns_config_values():
    config = {
        "rename_key": {"A": "a"},
        "other_amr": ["AMR"],
        "other_non_amr": ["NON_AMR"],
        "oxa_phen_list": ["gene"],
        "email_address": "test@example.org",
    }

    assert Update._get_keys(config) == (
        config["rename_key"],
        config["other_amr"],
        config["other_non_amr"],
        config["oxa_phen_list"],
        config["email_address"],
    )


def test_get_vars_reads_the_update_configuration(tmp_path, monkeypatch):
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    config = {
        "rename_key": {},
        "other_amr": [],
        "other_non_amr": [],
        "oxa_phen_list": [],
        "email_address": "test@example.org",
    }
    (db_dir / "update_vars.json").write_text(json.dumps(config))
    monkeypatch.setattr(Update, "__file__", str(tmp_path / "Update.py"))

    assert Update._get_vars() == Update._get_keys(config)


def test_get_vars_exits_when_configuration_is_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(Update, "__file__", str(tmp_path / "Update.py"))

    with pytest.raises(SystemExit):
        Update._get_vars()


def test_make_key_uses_protein_accession_or_point_mutation_key():
    catalog = pandas.DataFrame(
        [
            {
                "allele": "geneA",
                "whitelisted_taxa": "taxon",
                "refseq_nucleotide_accession": "NC_1",
                "refseq_protein_accession": "WP_1",
                "genbank_protein_accession": "GB_1",
                "subtype": "AMR",
            },
            {
                "allele": "geneB",
                "whitelisted_taxa": "taxon",
                "refseq_nucleotide_accession": "NC_2",
                "refseq_protein_accession": "",
                "genbank_protein_accession": "GB_2",
                "subtype": "POINT",
            },
        ]
    )

    result = Update._make_key(catalog)

    assert result["key"].tolist() == ["WP_1", "geneB_taxon_NC_2"]
    assert result["mut_acc"].tolist() == ["geneA_taxon_NC_1", "geneB_taxon_NC_2"]


def test_make_key_preserves_an_existing_key():
    catalog = pandas.DataFrame({"key": ["existing"]})

    assert Update._make_key(catalog) is catalog
    assert catalog["key"].tolist() == ["existing"]


def test_check_existing_returns_catalog_path_or_false(tmp_path, monkeypatch):
    monkeypatch.setattr(Update, "__file__", str(tmp_path / "Update.py"))
    db_dir = tmp_path / "db"
    db_dir.mkdir()

    assert Update._check_existing() is False

    catalog_path = db_dir / "refgenes_latest.csv"
    catalog_path.touch()
    assert Update._check_existing() == catalog_path


def test_get_previous_refgenes_loads_or_returns_false(tmp_path, monkeypatch):
    catalog_path = tmp_path / "refgenes_latest.csv"
    catalog_path.touch()
    catalog = pandas.DataFrame({"key": ["WP_1"]})
    monkeypatch.setattr(Update, "_check_existing", lambda: catalog_path)
    monkeypatch.setattr(Update.pandas, "read_csv", lambda path: catalog)

    assert Update._get_previous_refgenes().equals(catalog)

    monkeypatch.setattr(Update, "_check_existing", lambda: False)
    assert Update._get_previous_refgenes() is False


def test_simple_classification_helpers():
    row = _row(**{"class": "PHENICOL/OXAZOLIDINONE", "subclass": "CHLORAMPHENICOL"})

    assert Update._capitalise("BETA-LACTAM/CARBAPENEM") == "Beta-lactam/Carbapenemase"
    assert Update._oxa_phen(row) == ("Oxazolidinone/Phenicol", "Chloramphenicol")
    assert Update._other_antimicrobials(row) == (
        "Other antimicrobial",
        "Chloramphenicol",
    )
    assert Update._other_non_antimicrobials(row) == (
        "Other non-antimicrobial",
        "Chloramphenicol",
    )
    assert Update._rename(
        {"FLUOROQUINOLONE": "Quinolone", "CHLORAMPHENICOL": "Phenicol"},
        row,
    ) == ("Phenicol/Oxazolidinone", "Phenicol")


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (
            _row(subtype="AMR-SUSCEPTIBLE"),
            ("Beta-lactam", "Beta-lactam (Susceptible)"),
        ),
        (
            _row(product_name="carbapenem-hydrolyzing enzyme"),
            ("Beta-lactam", "Carbapenemase"),
        ),
        (
            _row(
                product_name="metallo-beta-lactamase",
                subclass="CARBAPENEM",
            ),
            ("Beta-lactam", "Carbapenemase (MBL)"),
        ),
        (
            _row(product_name="OXA-51 family protein"),
            ("Beta-lactam", "Carbapenemase (OXA-51 family)"),
        ),
        (
            _row(gene_family="blaZ"),
            ("Beta-lactam", "Penicillin resistance (Staphylococcus aureus)"),
        ),
        (
            _row(product_name="class C beta-lactamase"),
            ("Beta-lactam", "AmpC"),
        ),
        (
            _row(gene_family="blaKPC", subclass="ESBL"),
            ("Beta-lactam", "ESBL (KPC variant)"),
        ),
        (
            _row(product_name="extended-spectrum beta-lactamase"),
            ("Beta-lactam", "ESBL"),
        ),
        (
            _row(**{"class": "BETA-LACTAM", "subclass": "CARBAPENEM"}),
            ("Beta-lactam", "Carbapenemase"),
        ),
    ],
)
def test_beta_lactam_classification(row, expected):
    assert Update._beta_lactams(row) == expected


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (
            _row(
                **{
                    "class": "AMINOGLYCOSIDE",
                    "product_name": "16S rRNA methyltransferase",
                }
            ),
            ("Aminoglycoside", "Aminoglycosides (Ribosomal methyltransferase)"),
        ),
        (
            _row(
                **{
                    "class": "AMINOGLYCOSIDE",
                    "product_name": "aminoglycoside resistance",
                }
            ),
            ("Aminoglycoside", "Subclass"),
        ),
        (
            _row(gene_family="cfr", **{"class": "", "subclass": "PHENICOL"}),
            ("Multidrug", "Phenicol"),
        ),
        (
            _row(
                **{
                    "class": "INTIMIN",
                    "subclass": "EAE",
                    "type": "VIRULENCE",
                }
            ),
            ("Virulence", "Intimin_eae"),
        ),
        (
            _row(
                **{
                    "class": "STX1",
                    "subclass": "STX1A",
                    "type": "VIRULENCE",
                }
            ),
            ("Virulence", "Stx1a"),
        ),
        (
            _row(**{"class": "", "subclass": "", "type": "VIRULENCE"}),
            ("Virulence", "Other"),
        ),
    ],
)
def test_special_classification_helpers(row, expected):
    if row["class"] == "AMINOGLYCOSIDE":
        actual = Update._aminoglycosides(row)
    elif row["gene_family"] == "cfr":
        actual = Update.cfr(row)
    else:
        actual = Update.virulence(row)

    assert actual == expected


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (
            _row(
                **{
                    "class": "AMINOGLYCOSIDE",
                    "product_name": "rRNA methyltransferase",
                }
            ),
            ("Aminoglycoside", "Aminoglycosides (Ribosomal methyltransferase)"),
        ),
        (
            _row(gene_family="cfr", **{"class": "", "subclass": "PHENICOL"}),
            ("Multidrug", "Phenicol"),
        ),
        (
            _row(**{"class": "BETA-LACTAM", "subtype": "AMR-SUSCEPTIBLE"}),
            ("Beta-lactam", "Beta-lactam (Susceptible)"),
        ),
        (
            _row(**{"class": "OTHER_AMR"}),
            ("Other antimicrobial", "Subclass"),
        ),
        (
            _row(
                **{
                    "class": "INTIMIN",
                    "subclass": "EAE",
                    "type": "VIRULENCE",
                }
            ),
            ("Virulence", "Intimin_eae"),
        ),
        (
            _row(**{"class": "OTHER_NON_AMR"}),
            ("Other non-antimicrobial", "Subclass"),
        ),
        (
            _row(**{"class": "FLUOROQUINOLONE"}),
            ("Quinolone", "Mapped subclass"),
        ),
        (
            _row(**{"class": "MULTIDRUG"}),
            ("Multidrug", "Subclass"),
        ),
        (
            _row(
                **{
                    "class": "",
                    "subclass": "",
                    "product_name": "multidrug efflux pump",
                }
            ),
            ("Multidrug", "Other"),
        ),
        (
            _row(**{"class": "", "subclass": ""}),
            ("Other", "Other"),
        ),
        (_row(**{"class": "UNMAPPED_CLASS"}), ("Unmapped_class", "Subclass")),
    ],
)
def test_logic_classifies_each_catalog_category(row, expected):
    result = Update._logic(
        [row],
        other_amr=["OTHER_AMR"],
        other_non_amr=["OTHER_NON_AMR"],
        rename_key={
            "FLUOROQUINOLONE": "Quinolone",
            "SUBCLASS": "Mapped subclass",
        },
    )

    assert (result[0]["enhanced_class"], result[0]["enhanced_subclass"]) == expected


def test_make_dict_applies_classification_to_dataframe():
    catalog = pandas.DataFrame([_row(**{"class": "MULTIDRUG"})])

    result = Update._make_dict(catalog, [], [], {})

    assert result[0]["enhanced_class"] == "Multidrug"
    assert result[0]["enhanced_subclass"] == "Subclass"


def test_updated_entries_marks_changed_catalog_rows():
    current = pandas.DataFrame(
        [
            {"key": "same", "class": "AMR", "subclass": "Beta", "Status": ""},
            {"key": "changed", "class": "AMR", "subclass": "Updated", "Status": ""},
            {"key": "new", "class": "AMR", "subclass": "New", "Status": ""},
        ]
    )
    previous = pandas.DataFrame(
        [
            {"key": "same", "class_new": "AMR", "subclass_new": "Beta"},
            {"key": "changed", "class_new": "AMR", "subclass_new": "Old"},
        ]
    )

    result = Update._compare_to_existing(current, previous).set_index("key")

    assert result.loc["same", "Status"] == "existing"
    assert result.loc["changed", "Status"] == "updated"
    assert result.loc["changed", "Previous_subclass"] == "Old"
    assert result.loc["new", "Status"] == "new"


def test_compare_to_existing_without_previous_catalog_returns_current_rows():
    current = [{"key": "new", "class": "AMR"}]

    assert Update._compare_to_existing(current, False).to_dict("records") == current


def test_get_new_catalog_returns_success_and_exits_on_failure(monkeypatch):
    monkeypatch.setattr(
        Update.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0),
    )
    assert Update._get_new_catalog()

    monkeypatch.setattr(
        Update.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=1, stderr="download failed"),
    )
    with pytest.raises(SystemExit):
        Update._get_new_catalog()


def test_open_catalog_downloads_and_builds_keys(monkeypatch):
    catalog = pandas.DataFrame(
        [
            {
                "allele": "gene",
                "whitelisted_taxa": "taxon",
                "refseq_nucleotide_accession": "NC_1",
                "refseq_protein_accession": "WP_1",
                "genbank_protein_accession": "",
                "subtype": "AMR",
            }
        ]
    )
    monkeypatch.setattr(Update, "_get_new_catalog", lambda: True)
    monkeypatch.setattr(Update.pandas, "read_csv", lambda *args, **kwargs: catalog)

    result = Update._open_catalog()

    assert result["key"].tolist() == ["WP_1"]


def test_email_and_archive_delegate_to_subprocess(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(Update, "_get_date", lambda: "2026-01-02")
    monkeypatch.setattr(
        Update.subprocess,
        "run",
        lambda *args, **kwargs: calls.append((args, kwargs))
        or SimpleNamespace(returncode=0, stderr=""),
    )

    Update._archive_old_ref_catalog()
    Update._email("test@example.org", "catalog.csv")

    assert calls[0][0][0] == [
        "cp",
        str(pathlib.Path(Update.__file__).parent / "db" / "refgenes_latest.csv"),
        str(pathlib.Path(Update.__file__).parent / "db" / "refgenes_latest.csv.2026-01-02"),
    ]
    assert "mailx" in calls[1][0][0]


def test_save_df_writes_dated_csv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(Update, "_get_date", lambda: "2026-01-02")
    catalog = pandas.DataFrame({"gene": ["geneA"]})

    result = Update._save_df(catalog)

    assert result == "refgenes_2026-01-02.csv"
    assert pandas.read_csv(tmp_path / result)["gene"].tolist() == ["geneA"]


def test_create_refgenes_runs_update_steps(monkeypatch):
    calls = []
    monkeypatch.setattr(Update, "_get_vars", lambda: ("rename", "amr", "non-amr", [], "email"))
    monkeypatch.setattr(Update, "_open_catalog", lambda: "downloaded")
    monkeypatch.setattr(Update, "_make_dict", lambda **kwargs: calls.append(("make", kwargs)) or "classified")
    monkeypatch.setattr(Update, "_get_previous_refgenes", lambda: "previous")
    monkeypatch.setattr(
        Update,
        "_compare_to_existing",
        lambda **kwargs: calls.append(("compare", kwargs)) or "compared",
    )
    monkeypatch.setattr(Update, "_save_df", lambda df: calls.append(("save", df)) or "catalog.csv")
    monkeypatch.setattr(Update, "_email", lambda **kwargs: calls.append(("email", kwargs)))

    Update.create_refgenes()

    assert [call[0] for call in calls] == ["make", "compare", "save", "email"]
    assert calls[2] == ("save", "compared")
    assert calls[3][1] == {"adrs": "email", "pth": "catalog.csv"}
