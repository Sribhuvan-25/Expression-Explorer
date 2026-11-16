"""GSE66586 loader: label parsing, the 5/3/2 guard, symbol normalisation,
and probe collapse -- all against small synthetic GEO files in tmp_path,
so nothing here touches the network."""
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from app.ingest import gse66586

# The real series matrix's sample header, in its column order.
REAL_SAMPLES = [
    # gsm, title, source_name, cell type, ATCC provider
    ("GSM1625935", "hBM-MSCs", "Control_1", "bone marrow mesenchymal stem cells", "ATCC:PCS-500-012"),
    ("GSM1625937", "HCN1-A", "Control_2", "cortical neuronal", "ATCC:CRL\u00ad10442"),
    ("GSM1625938", "BE(2)c", "MYCN_Amplified_1", "neuroblastoma", "ATCC:CRL-2268"),
    ("GSM1625939", "SK-N-DZ", "MYCN_Amplified_2", "neuroblastoma", "ATCC:CRL-2149"),
    ("GSM1625940", "IMR-32", "MYCN_Amplified_3", "neuroblastoma", "ATCC:CCL-127"),
    ("GSM1625941", "SK-N-BE(2)c", "MYCN_Amplified_4", "neuroblastoma", "ATCC:CRL-2271"),
    ("GSM1625942", "CHP-212", "MYCN_Amplified_5", "neuroblastoma", "ATCC:CRL-2273"),
    ("GSM1625944", "SK-N-AS", "MYCN_NotAmplified_1", "neuroblastoma", "ATCC:CRL-2137"),
    ("GSM1625945", "SH-SY5Y", "MYCN_NotAmplified_2", "neuroblastoma", "ATCC:CRL-2266"),
    ("GSM1625946", "SK-N-SH", "MYCN_NotAmplified_3", "neuroblastoma", "ATCC:HTB-11"),
]


def _q(values):
    return "\t".join(f'"{v}"' for v in values)


def _write_series_matrix(path: Path, samples, probes: dict[str, list[float]]) -> None:
    gsm = [s[0] for s in samples]
    lines = [
        '!Series_title\t"Neuroblastoma cell line gene expression profiling"',
        '!Series_geo_accession\t"GSE66586"',
        f'!Series_sample_id\t"{" ".join(gsm)} "',
        "",
        "!Sample_title\t" + _q(s[1] for s in samples),
        "!Sample_geo_accession\t" + _q(gsm),
        "!Sample_source_name_ch1\t" + _q(s[2] for s in samples),
        "!Sample_organism_ch1\t" + _q("Homo sapiens" for _ in samples),
        "!Sample_characteristics_ch1\t" + _q(f"cell type: {s[3]}" for s in samples),
        "!Sample_biomaterial_provider_ch1\t" + _q(s[4] for s in samples),
        "!Sample_platform_id\t" + _q("GPL6244" for _ in samples),
        "!series_matrix_table_begin",
        '"ID_REF"\t' + _q(gsm),
    ]
    for probe, values in probes.items():
        lines.append(probe + "\t" + "\t".join(str(v) for v in values))
    lines.append("!series_matrix_table_end")
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def _write_annot(path: Path, probe_symbols: dict[str, str]) -> None:
    lines = [
        "^Annotation",
        "!Annotation_platform = GPL6244",
        "#ID = ID from Platform data table",
        "!platform_table_begin",
        "ID\tGene title\tGene symbol\tGene ID",
    ]
    for probe, symbol in probe_symbols.items():
        lines.append(f"{probe}\ttitle\t{symbol}\t1")
    lines.append("!platform_table_end")
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def _flat(level: float) -> list[float]:
    return [level] * len(REAL_SAMPLES)


@pytest.fixture
def geo_files(tmp_path, monkeypatch):
    """Point the loader at tmp_path. The files already exist there, so
    `download()` short-circuits; the URLs are also made unreachable so a
    regression that tried to fetch would fail rather than hit GEO."""
    monkeypatch.setattr(gse66586, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(gse66586, "SERIES_MATRIX_URL", "http://invalid.test/matrix.gz")
    monkeypatch.setattr(gse66586, "PLATFORM_ANNOT_URL", "http://invalid.test/annot.gz")

    def write(samples=REAL_SAMPLES, probes=None, annot=None):
        probes = probes or {"7892501": _flat(7.0)}
        annot = annot or {"7892501": "MYCN"}
        _write_series_matrix(tmp_path / "GSE66586_series_matrix.txt.gz", samples, probes)
        _write_annot(tmp_path / "GPL6244.annot.gz", annot)

    return write


def test_labels_parse_to_the_three_mycn_groups(geo_files):
    geo_files()
    ds = gse66586.load()

    groups = {r.sample_id: r.group_columns for r in ds.samples.itertuples()}
    status = {s: g["mycn_status"] for s, g in groups.items()}
    assert sorted(status.values()).count("MYCN-amplified") == 5
    assert sorted(status.values()).count("MYCN non-amplified") == 3
    assert sorted(status.values()).count("non-malignant") == 2
    assert status["GSM1625935"] == "non-malignant"
    assert status["GSM1625940"] == "MYCN-amplified"
    assert status["GSM1625945"] == "MYCN non-amplified"

    assert groups["GSM1625940"]["cell_line"] == "IMR-32"
    assert groups["GSM1625937"]["cell_type"] == "cortical neuronal"
    # Identity comes from the ATCC catalogue, not the title: GSM1625941 is
    # titled "SK-N-BE(2)c" but is CRL-2271, the parental SK-N-BE(2).
    assert groups["GSM1625938"]["cell_line_identity"] == "SK-N-BE(2)-C"
    assert groups["GSM1625941"]["cell_line_identity"] == "SK-N-BE(2)"
    # The soft hyphen GEO stores in HCN-1A's catalogue number still resolves.
    assert groups["GSM1625937"]["cell_line_identity"] == "HCN-1A"
    assert list(ds.matrix.columns) == [s[0] for s in REAL_SAMPLES]


def test_regrouped_matrix_fails_loudly(geo_files):
    """A re-issued matrix that labels SK-N-AS as amplified (6/2/2) must not
    load as a silently regrouped cohort."""
    samples = list(REAL_SAMPLES)
    gsm, title, _, ctype, atcc = samples[7]
    samples[7] = (gsm, title, "MYCN_Amplified_6", ctype, atcc)
    geo_files(samples=samples)

    with pytest.raises(ValueError, match="expected"):
        gse66586.load()


def test_unrecognised_source_name_fails_loudly(geo_files):
    samples = list(REAL_SAMPLES)
    gsm, title, _, ctype, atcc = samples[0]
    samples[0] = (gsm, title, "Normal_1", ctype, atcc)
    geo_files(samples=samples)

    with pytest.raises(ValueError, match="Unrecognised GSE66586 source name"):
        gse66586.load()


def test_legacy_symbols_are_renamed_and_recorded_as_aliases(geo_files):
    """GPL6244's annotation still names CIP2A 'KIAA1524' and PABIR1
    'FAM122A'. The matrix must carry the current symbol (or a CIP2A query
    returns 'gene not found') and record the legacy name it came from."""
    geo_files(
        probes={"p1": _flat(8.8), "p2": _flat(8.9), "p3": _flat(10.0)},
        annot={"p1": "KIAA1524", "p2": "FAM122A///FAM122A", "p3": "MYCN"},
    )
    ds = gse66586.load()

    assert {"CIP2A", "PABIR1", "MYCN"} == set(ds.matrix.index)
    assert "KIAA1524" not in ds.matrix.index
    assert "FAM122A" not in ds.matrix.index
    assert ds.features.loc["CIP2A", "aliases"] == ["KIAA1524"]
    assert ds.features.loc["PABIR1", "aliases"] == ["FAM122A"]
    assert ds.features.loc["MYCN", "aliases"] == []


@pytest.mark.parametrize("legacy_level, current_level", [(9.0, 4.0), (4.0, 9.0)])
def test_legacy_current_collision_keeps_highest_mean_probe(geo_files, legacy_level, current_level):
    """If a platform carried a probe under both KIAA1524 and CIP2A, they are
    one gene: exactly one CIP2A row, from whichever probe has the higher
    mean -- whichever name it was annotated under."""
    geo_files(
        probes={"legacy": _flat(legacy_level), "current": _flat(current_level)},
        annot={"legacy": "KIAA1524", "current": "CIP2A"},
    )
    ds = gse66586.load()

    assert list(ds.matrix.index) == ["CIP2A"]
    assert (ds.matrix.loc["CIP2A"] == max(legacy_level, current_level)).all()
    assert ds.features.loc["CIP2A", "aliases"] == ["KIAA1524"]


def test_multi_probe_gene_collapses_to_highest_mean_not_group_separation(geo_files):
    """ARPP19 has two transcript clusters on GPL6244. The kept probe must be
    the highest-mean one even when the other separates MYCN groups
    perfectly -- choosing on separation would be circular (METHODS 7.1)."""
    amp = [s[2].startswith("MYCN_Amplified") for s in REAL_SAMPLES]
    separating = [6.0 if a else 2.0 for a in amp]  # mean 3.6, perfect split
    geo_files(
        probes={"high": _flat(11.0), "separating": separating},
        annot={"high": "ARPP19", "separating": "ARPP19"},
    )
    ds = gse66586.load()

    assert list(ds.matrix.index) == ["ARPP19"]
    assert (ds.matrix.loc["ARPP19"] == 11.0).all()
