"""
GSE324984 loader: the properties that keep its numbers attributable to
the right samples and genes. Synthetic GEO files only -- no network.

Each test corresponds to something found in the real submission:
  - matrix columns are library names, not GSMs, and must be mapped;
  - the TPM file is a truncated prefix of the counts file (Excel limit);
  - one Ensembl id is named differently in the two files (TP63/KDM5A);
  - symbols repeat across Ensembl ids (e.g. CAST);
  - a swapped column would flow silently into every result.
"""
import gzip

import numpy as np
import pandas as pd
import pytest

from app.ingest import gse324984

LIBS = {
    "GSM9589988": ("C4-2B control biological replicate 1", "C4-2B", "empty vector", "C4-2B_Ctrl_1"),
    "GSM9589989": ("C4-2B control biological replicate 2", "C4-2B", "empty vector", "C4-2B_Ctrl_2"),
    "GSM9589990": ("C4-2B MYCL overexpression replicate 1", "C4-2B", "MYCL overexpression", "C4-2B_MYCL_1"),
    "GSM9589991": ("C4-2B MYCL overexpression replicate 2", "C4-2B", "MYCL overexpression", "C4-2B_MYCL_2"),
    "GSM9589992": ("H660 shControl replicate 1", "NCI-H660", "shCtrl", "shCtrl_H660_1"),
    "GSM9589993": ("H660 shControl replicate 2", "NCI-H660", "shCtrl", "shCtrl_H660_2"),
    "GSM9589994": ("H660 MYCL knockdown replicate 1", "NCI-H660", "shMYCL", "shMYCL_H660_1"),
    "GSM9589995": ("H660 MYCL knockdown replicate 2", "NCI-H660", "shMYCL", "shMYCL_H660_2"),
}
# The real files order columns H660 first; the series matrix lists C4-2B first.
FILE_COLUMN_ORDER = [v[3] for v in list(LIBS.values())[4:]] + [v[3] for v in list(LIBS.values())[:4]]
MYCL_ID = "ENSG00000116990"
N_GENES = 600
N_TPM = 500  # TPM file is a prefix of the counts file


def _write_series_matrix(path):
    gsms = list(LIBS)
    def row(key, values):
        return key + "\t" + "\t".join(f'"{v}"' for v in values) + "\n"
    with gzip.open(path, "wt") as fh:
        fh.write('!Series_title\t"test"\n')
        fh.write(row("!Sample_title", [LIBS[g][0] for g in gsms]))
        fh.write(row("!Sample_geo_accession", gsms))
        fh.write(row("!Sample_characteristics_ch1", [f"cell line: {LIBS[g][1]}" for g in gsms]))
        fh.write(row("!Sample_characteristics_ch1", ["cell type: Prostate cancer cells"] * 8))
        fh.write(row("!Sample_characteristics_ch1", [f"treatment: {LIBS[g][2]}" for g in gsms]))
        fh.write(row("!Sample_description", [f"Library name: {LIBS[g][3]}" for g in gsms]))
        fh.write('!series_matrix_table_begin\n"ID_REF"\n!series_matrix_table_end\n')


def _synthetic(seed=0):
    rng = np.random.default_rng(seed)
    ids = [f"ENSG{i:011d}" for i in range(1, N_GENES + 1)]
    ids[10] = MYCL_ID
    names = [f"GENE{i}" for i in range(N_GENES)]
    names[10] = "MYCL"
    base = rng.gamma(1.0, 200.0, size=N_GENES)
    counts = pd.DataFrame(
        {lib: rng.poisson(base * rng.uniform(0.9, 1.1)) for lib in FILE_COLUMN_ORDER}
    ).astype(float)
    # MYCL: up with overexpression, down with knockdown.
    for lib, value in {"C4-2B_Ctrl_1": 60, "C4-2B_Ctrl_2": 70, "C4-2B_MYCL_1": 9000, "C4-2B_MYCL_2": 9500,
                       "shCtrl_H660_1": 900, "shCtrl_H660_2": 800, "shMYCL_H660_1": 90, "shMYCL_H660_2": 100}.items():
        counts.loc[10, lib] = value
    tpm = counts / counts.sum() * 1e6
    return ids, names, counts, tpm


def _write_matrix(path, ids, names, values):
    frame = pd.DataFrame({"gene_id": [f"{i}.3" for i in ids], "gene_name": names})
    frame = pd.concat([frame, values.reset_index(drop=True)], axis=1)
    frame.to_csv(path, sep="\t", index=False, compression="gzip")


@pytest.fixture()
def write_dataset(tmp_path, monkeypatch):
    monkeypatch.setattr(gse324984, "CACHE_DIR", tmp_path)

    def write(mutate_tpm=None, tpm_names=None, count_names=None, mutate_counts=None):
        ids, names, counts, tpm = _synthetic()
        if mutate_counts:
            mutate_counts(counts)
            tpm = counts / counts.sum() * 1e6
        _write_series_matrix(tmp_path / "GSE324984_series_matrix.txt.gz")
        tpm = tpm.iloc[:N_TPM].copy()
        if mutate_tpm:
            mutate_tpm(tpm)
        _write_matrix(tmp_path / "GSE324984_RNAseq_gene_counts.txt.gz", ids, count_names or names, counts)
        _write_matrix(tmp_path / "GSE324984_RNAseq_TPM_matrix.txt.gz", ids[:N_TPM], (tpm_names or names)[:N_TPM], tpm)
        return ids, names, counts, tpm

    return write


def test_columns_are_mapped_by_library_name_not_position(write_dataset):
    ids, _, _, tpm = write_dataset()
    ds = gse324984.load()
    assert list(ds.matrix.columns) == list(LIBS)  # GSM ids, series order
    # File column 0 is shCtrl_H660_1 == GSM9589992, not GSM9589988.
    assert ds.matrix.loc[ids[0], "GSM9589992"] == pytest.approx(tpm.iloc[0]["shCtrl_H660_1"])
    arms = ds.samples["group_columns"].apply(lambda d: d["arm"]).value_counts()
    assert set(arms) == {2}
    assert ds.samples.loc["GSM9589994", "group_columns"]["arm"] == "NCI-H660 MYCL knockdown (shMYCL)"


def test_swapped_columns_fail_the_mycl_label_check(write_dataset):
    def swap(tpm):
        tpm[["C4-2B_Ctrl_1", "C4-2B_MYCL_1"]] = tpm[["C4-2B_MYCL_1", "C4-2B_Ctrl_1"]].to_numpy()
        tpm[["C4-2B_Ctrl_2", "C4-2B_MYCL_2"]] = tpm[["C4-2B_MYCL_2", "C4-2B_Ctrl_2"]].to_numpy()
    write_dataset(mutate_tpm=swap)
    with pytest.raises(ValueError, match="contradicts the arm labels"):
        gse324984.load()


def test_truncated_tpm_file_keeps_full_gene_set_for_statistics(write_dataset):
    ids, *_ = write_dataset()
    ds = gse324984.load()
    assert ds.matrix.shape == (N_TPM, 8)
    for contrast in ds.contrasts.values():
        assert len(contrast.table) == N_GENES
        assert not contrast.table.loc[ids[N_TPM:], "in_expression_matrix"].any()
        assert contrast.table.loc[ids[:N_TPM], "in_expression_matrix"].all()
        # Genes beyond the TPM truncation are genuinely tested, not padded.
        tail = contrast.table.loc[ids[N_TPM:]]
        assert tail["tested"].sum() > 0.5 * len(tail)
        assert tail.loc[tail["tested"].astype(bool), "p_value"].notna().all()
    assert f"{N_GENES - N_TPM:,} of {N_GENES:,} genes have no TPM" in ds.source.notes
    assert "tested genes that have TPM" in ds.source.notes


def test_contrasts_use_only_their_own_cell_line(write_dataset):
    write_dataset()
    ds = gse324984.load()
    oe, kd = ds.contrasts["c42b_mycl_oe"], ds.contrasts["h660_mycl_kd"]
    assert oe.reference_samples == ["GSM9589988", "GSM9589989"]
    assert oe.treatment_samples == ["GSM9589990", "GSM9589991"]
    assert kd.reference_samples == ["GSM9589992", "GSM9589993"]
    assert kd.treatment_samples == ["GSM9589994", "GSM9589995"]
    assert oe.table.loc[MYCL_ID, "log2_fc"] > 0 > kd.table.loc[MYCL_ID, "log2_fc"]


def test_gene_name_conflict_resolves_to_tpm_name_and_is_disclosed(write_dataset):
    ids, names, *_ = _synthetic()
    count_names = list(names)
    count_names[3] = "KDM5A"  # the real file's TP63 row
    tpm_names = list(names)
    tpm_names[3] = "TP63"
    write_dataset(tpm_names=tpm_names, count_names=count_names)
    ds = gse324984.load()
    assert ds.features.loc[ids[3], "symbol"] == "TP63"
    assert ds.contrasts["c42b_mycl_oe"].table.loc[ids[3], "symbol"] == "TP63"
    assert f"{ids[3]} TPM 'TP63' vs counts 'KDM5A'" in ds.source.notes


def test_many_name_conflicts_mean_misaligned_files(write_dataset):
    _, names, *_ = _synthetic()
    count_names = [f"OTHER{i}" if 20 <= i < 30 else n for i, n in enumerate(names)]
    write_dataset(count_names=count_names)
    with pytest.raises(ValueError, match="misaligned"):
        gse324984.load()


def test_duplicate_symbol_resolves_to_highest_mean_tpm_row(write_dataset):
    """Within the TPM rows: the higher-count row keeps the name, the other
    shows its Ensembl id and carries NO alias for the shared name (an alias
    there made the legacy-name fallback resolve to the secondary locus)."""
    _, names, counts, _ = _synthetic()
    dup = list(names)
    dup[40] = dup[41] = "DUPGENE"
    ids, *_ = write_dataset(tpm_names=dup, count_names=dup)
    ds = gse324984.load()
    # Both loci are in the TPM file, so TPM (not raw count) decides.
    hi, lo = (ids[40], ids[41]) if ds.matrix.loc[ids[40]].mean() > ds.matrix.loc[ids[41]].mean() else (ids[41], ids[40])
    assert ds.features.loc[hi, "symbol"] == "DUPGENE"
    assert ds.features.loc[lo, "symbol"] == lo
    assert ds.features.loc[lo, "aliases"] == []
    assert (ds.features["symbol"] == "DUPGENE").sum() == 1
    table = ds.contrasts["c42b_mycl_oe"].table
    assert table.loc[hi, "symbol"] == "DUPGENE" and table.loc[lo, "symbol"] == lo
    assert table.loc[lo, "gene_name"] == "DUPGENE" and table.loc[lo, "n_loci_with_name"] == 2


def test_expressed_copy_beyond_the_tpm_truncation_wins_the_name(write_dataset):
    """The real CAST case: the expressed locus exists only in the counts
    (past the TPM file's cut); a near-silent TPM row shares its name. The
    name must go to the expressed locus, in every view."""
    tail = N_TPM + 20

    def shape(counts):
        counts.loc[30, :] = 1.0  # near-silent TPM-file row
        counts.loc[tail, :] = 5000.0  # expressed, count-only

    _, names, *_ = _synthetic()
    dup = list(names)
    dup[30] = dup[tail] = "CASTLIKE"
    ids, *_ = write_dataset(tpm_names=dup, count_names=dup, mutate_counts=shape)
    ds = gse324984.load()
    assert ds.features.loc[ids[30], "symbol"] == ids[30]  # TPM row gives up the name
    assert "CASTLIKE" not in set(ds.features["symbol"])  # no matrix row claims it
    table = ds.contrasts["c42b_mycl_oe"].table
    assert table.loc[ids[tail], "symbol"] == "CASTLIKE"
    assert not table.loc[ids[tail], "in_expression_matrix"]


def test_label_check_also_guards_the_counts_file(write_dataset):
    """A column swap in the counts file alone would reach the statistics."""
    def swap(counts):
        counts[["C4-2B_Ctrl_1", "C4-2B_MYCL_1"]] = counts[["C4-2B_MYCL_1", "C4-2B_Ctrl_1"]].to_numpy()
        counts[["C4-2B_Ctrl_2", "C4-2B_MYCL_2"]] = counts[["C4-2B_MYCL_2", "C4-2B_Ctrl_2"]].to_numpy()

    ids, names, counts, tpm = _synthetic()
    write_dataset()  # TPM consistent ...
    # ... then rewrite only the counts file with swapped columns.
    swap(counts)
    _write_matrix(gse324984.CACHE_DIR / "GSE324984_RNAseq_gene_counts.txt.gz", ids, names, counts)
    with pytest.raises(ValueError, match=r"CPM \(counts file\)"):
        gse324984.load()


def test_legacy_symbols_are_recorded_as_aliases(write_dataset):
    _, names, *_ = _synthetic()
    named = list(names)
    named[50], named[51] = "CIP2A", "PABIR1"
    ids, *_ = write_dataset(tpm_names=named, count_names=named)
    ds = gse324984.load()
    assert ds.features.loc[ids[50], "aliases"] == ["KIAA1524"]
    assert ds.features.loc[ids[51], "aliases"] == ["FAM122A"]


def test_tpm_not_raw_count_decides_between_two_tpm_loci(write_dataset):
    """The C4orf36 case: near-equal counts, clearly different TPM (a long
    locus collects reads without being more expressed per base)."""
    _, names, *_ = _synthetic()
    named = list(names)
    named[40] = named[41] = "LENGTHY"

    def same_counts(counts):
        counts.loc[40, :] = 1000.0
        counts.loc[41, :] = 1010.0  # marginally more reads ...

    def tpm_says_otherwise(tpm):
        tpm.loc[40, :] = 5.0
        tpm.loc[41, :] = 1.0  # ... but far lower TPM

    ids, *_ = write_dataset(tpm_names=named, count_names=named, mutate_counts=same_counts, mutate_tpm=tpm_says_otherwise)
    ds = gse324984.load()
    assert ds.features.loc[ids[40], "symbol"] == "LENGTHY"
    assert ds.features.loc[ids[41], "symbol"] == ids[41]
