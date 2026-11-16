"""
Build the biomarker database for Biomarker Explorer.

For each drug, this script asks one question across cell lines:
does higher expression of a gene go with a lower AUC (more sensitive)
or a higher AUC (more resistant)?

Run it from the project folder, with the virtual environment turned on:

    python scripts/build_database.py

It writes three parquet files in processed/:

    biomarkers.parquet      top 50 sensitivity and top 50 resistance genes per drug
    drug_response.parquet   each drug's AUC, with the cell line name and cancer type
    expression.parquet      expression of the genes that appear in the biomarker table
"""

from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
PROCESSED_DIR = ROOT / "processed"

# A drug needs measurements in at least this many shared cell lines.
MIN_CELL_LINES = 100
# How many genes to keep at each end of the correlation list.
TOP_N = 50
# Parquet files should stay under this size.
MAX_FILE_BYTES = 50 * 1024 * 1024

EXPRESSION_FILE = "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
BOOKKEEPING_COLUMNS = [
    "Unnamed: 0",
    "SequencingID",
    "ModelConditionID",
    "IsDefaultEntryForMC",
    "IsDefaultEntryForModel",
]


def check_correlation_direction():
    """Fail early if the sign of the correlation is backwards.

    Five imaginary cell lines: expression rises while AUC falls.
    That pattern is sensitivity, so the correlation must be negative.
    """
    expression = np.array([[1.0], [2.0], [3.0], [4.0], [5.0]])
    auc = np.array([0.90, 0.70, 0.50, 0.20, 0.10])
    correlation = spearman_vs_genes(expression, auc)[0]
    reference, _ = stats.spearmanr(expression[:, 0], auc)
    if not np.isclose(correlation, reference) or correlation >= 0:
        raise RuntimeError(
            "The correlation sign check failed. High expression with a low "
            "AUC should produce a negative Spearman correlation."
        )


def spearman_vs_genes(expression_values, auc_values):
    """Spearman correlation of one AUC vector with every gene column.

    expression_values is shaped (cell lines, genes). Ranking each gene
    across cell lines is one array operation, then one correlation is
    computed for every gene at the same time.
    """
    gene_ranks = stats.rankdata(expression_values, axis=0, method="average")
    return correlations_from_ranks(gene_ranks, auc_values)


def correlations_from_ranks(gene_ranks, auc_values):
    """Pearson correlation between ranked AUC and each ranked gene."""
    auc_ranks = stats.rankdata(np.asarray(auc_values, dtype=float), method="average")
    gene_centered = gene_ranks - gene_ranks.mean(axis=0)
    auc_centered = auc_ranks - auc_ranks.mean()
    numerator = gene_centered.T @ auc_centered
    gene_scale = np.sqrt(np.sum(gene_centered ** 2, axis=0))
    auc_scale = np.sqrt(np.sum(auc_centered ** 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = numerator / (gene_scale * auc_scale)
    correlation[~np.isfinite(correlation)] = np.nan
    return correlation


def p_values_for_correlations(correlation, n_cell_lines):
    """Two-sided p-values from the usual Spearman t approximation.

    This is the same approximation scipy.stats.spearmanr uses.
    A gene with no variation has an undefined correlation and a missing p-value.
    """
    p_values = np.full(np.shape(correlation), np.nan, dtype=float)
    if n_cell_lines < 3:
        return p_values

    finite = np.isfinite(correlation)
    perfect = finite & (np.abs(correlation) >= 1.0 - 1e-12)
    testable = finite & ~perfect
    p_values[perfect] = 0.0

    degrees_of_freedom = n_cell_lines - 2
    r = np.clip(correlation[testable], -0.999999999999, 0.999999999999)
    t_statistic = r * np.sqrt(degrees_of_freedom / (1.0 - r ** 2))
    p_values[testable] = 2.0 * stats.t.sf(np.abs(t_statistic), degrees_of_freedom)
    return p_values


def adjust_fdr(p_values):
    """Benjamini-Hochberg FDR across every gene tested for one drug.

    The adjustment uses the full gene list. The top-50 cutoff happens later,
    so the saved FDR values are not recomputed on only 100 genes.
    """
    adjusted = np.full(np.shape(p_values), np.nan, dtype=float)
    valid = np.isfinite(p_values)
    if np.any(valid):
        adjusted[valid] = stats.false_discovery_control(p_values[valid], method="bh")
    return adjusted


def top_genes_for_drug(gene_names, correlation, p_values, fdr, drug_name, n_cell_lines):
    """Keep the strongest negative and strongest positive correlations."""
    table = pd.DataFrame(
        {
            "gene": gene_names,
            "correlation": correlation,
            "p_value": p_values,
            "fdr": fdr,
        }
    )
    table = table.loc[np.isfinite(table["correlation"])].copy()

    sensitive = table.loc[table["correlation"] < 0].sort_values(
        ["correlation", "gene"], ascending=[True, True]
    ).head(TOP_N)
    resistant = table.loc[table["correlation"] > 0].sort_values(
        ["correlation", "gene"], ascending=[False, True]
    ).head(TOP_N)

    sensitive = sensitive.assign(direction="sensitivity")
    resistant = resistant.assign(direction="resistance")
    chosen = pd.concat([sensitive, resistant], ignore_index=True)
    chosen.insert(0, "drug_name", drug_name)
    chosen["n_cell_lines"] = int(n_cell_lines)
    return chosen


def require_data_files():
    required = [
        DATA_DIR / "Model.csv",
        DATA_DIR / EXPRESSION_FILE,
        DATA_DIR / "GDSC2AUCMatrix.csv",
        DATA_DIR / "GDSC2Log2ViabilityConditions.csv",
    ]
    missing = [path.name for path in required if not path.exists()]
    if missing:
        names = ", ".join(missing)
        raise FileNotFoundError(f"Missing files in {DATA_DIR}: {names}")


def load_catalog():
    """Cell-line names and cancer types. One row per ModelID."""
    catalog = pd.read_csv(
        DATA_DIR / "Model.csv",
        usecols=["ModelID", "CellLineName", "OncotreeLineage"],
    )
    if catalog["ModelID"].duplicated().any():
        raise ValueError("Model.csv has more than one row for a ModelID.")
    return catalog.set_index("ModelID")


def load_expression():
    """Default log2(TPM + 1) profile for each cell line, genes as columns."""
    print("  Reading gene expression. This file is large and can take a minute.", flush=True)
    expression = pd.read_csv(DATA_DIR / EXPRESSION_FILE)
    default = expression["IsDefaultEntryForModel"].eq("Yes")
    expression = expression.loc[default].copy()
    if expression["ModelID"].duplicated().any():
        raise ValueError("More than one default expression row exists for a cell line.")

    gene_columns = [column for column in expression.columns if column not in BOOKKEEPING_COLUMNS]
    # ModelID is one of the remaining columns; keep it as the index, not as a gene.
    gene_columns = [column for column in gene_columns if column != "ModelID"]
    expression = expression.set_index("ModelID")[gene_columns]
    expression.index = expression.index.astype(str)
    return expression


def load_drug_names():
    """One row per compound: the AUC column code and the drug name."""
    conditions = pd.read_csv(
        DATA_DIR / "GDSC2Log2ViabilityConditions.csv",
        usecols=["CompoundID", "CompoundName"],
    )
    lookup = conditions.drop_duplicates()
    names_per_id = lookup.groupby("CompoundID")["CompoundName"].nunique()
    ids_per_name = lookup.groupby("CompoundName")["CompoundID"].nunique()
    if names_per_id.max() != 1 or ids_per_name.max() != 1:
        raise ValueError("Drug codes and drug names are not in a one-to-one match.")
    return lookup.set_index("CompoundID")["CompoundName"].astype(str)


def load_auc(name_by_id):
    """AUC matrix with ModelID as the index and compound codes as columns."""
    auc = pd.read_csv(DATA_DIR / "GDSC2AUCMatrix.csv")
    auc = auc.rename(columns={"Unnamed: 0": "ModelID"}).set_index("ModelID")
    auc.index = auc.index.astype(str)
    missing_names = [column for column in auc.columns if column not in name_by_id.index]
    if missing_names:
        raise ValueError(f"{len(missing_names)} AUC columns have no drug name.")
    return auc


def shared_cell_lines(expression, auc, catalog):
    """Cell lines that have expression, a catalog entry, and a row in the AUC file."""
    shared_ids = expression.index.intersection(auc.index).intersection(catalog.index)
    expression = expression.loc[shared_ids]
    auc = auc.loc[shared_ids]
    catalog = catalog.loc[shared_ids]
    if not (expression.index.equals(auc.index) and expression.index.equals(catalog.index)):
        raise RuntimeError("Cell-line tables ended up in different orders.")
    return expression, auc, catalog


def varying_genes(expression):
    """Drop genes that are constant. Their correlation is undefined."""
    if expression.isna().any().any():
        raise ValueError("Expression values are missing for some genes. This script expects a complete matrix.")
    varies = expression.std(axis=0) > 0
    n_constant = int((~varies).sum())
    return expression.loc[:, varies], n_constant


def group_drugs_by_tested_cell_lines(auc, eligible_ids):
    """Drugs measured in exactly the same cell lines can share one ranking.

    Ranking every gene is the slow step. The correlation for each drug is quick
    once those ranks exist.
    """
    groups = {}
    for compound_id in eligible_ids:
        tested = auc[compound_id].notna().to_numpy()
        group = groups.get(tested.tobytes())
        if group is None:
            groups[tested.tobytes()] = {"tested": tested, "compound_ids": [compound_id]}
        else:
            group["compound_ids"].append(compound_id)
    return list(groups.values())


def score_all_drugs(expression, auc, name_by_id):
    """Return the top sensitivity and resistance genes for every eligible drug."""
    n_measured = auc.notna().sum(axis=0)
    eligible_ids = [compound_id for compound_id in auc.columns if int(n_measured[compound_id]) >= MIN_CELL_LINES]
    skipped = [name_by_id[compound_id] for compound_id in auc.columns if compound_id not in set(eligible_ids)]

    print()
    print(
        f"Scoring {expression.shape[1]:,} genes for {len(eligible_ids):,} drugs "
        f"with at least {MIN_CELL_LINES} cell lines."
    )
    if skipped:
        print(f"Skipping {len(skipped)} drugs with fewer cell lines: {', '.join(sorted(skipped))}")
    else:
        print(f"Every drug in the AUC file clears the {MIN_CELL_LINES}-cell-line minimum.")
    print("A line that names a drug means that drug has been scored.", flush=True)

    gene_names = expression.columns.to_numpy()
    expression_values = expression.to_numpy(dtype=np.float64)
    groups = group_drugs_by_tested_cell_lines(auc, eligible_ids)
    for group in groups:
        group["compound_ids"].sort(key=lambda compound_id: name_by_id[compound_id])
    groups.sort(key=lambda group: name_by_id[group["compound_ids"][0]])

    pieces = []
    done = 0
    n_drugs = len(eligible_ids)
    started = time.perf_counter()
    for group in groups:
        tested = group["tested"]
        n_cell_lines = int(tested.sum())
        compound_ids = group["compound_ids"]
        print(
            f"  Ranking {len(gene_names):,} genes across {n_cell_lines:,} cell lines "
            f"for {len(compound_ids)} drug{'s' if len(compound_ids) != 1 else ''}...",
            flush=True,
        )
        gene_ranks = stats.rankdata(expression_values[tested], axis=0, method="average")
        for compound_id in compound_ids:
            done += 1
            drug_name = name_by_id[compound_id]
            auc_values = auc[compound_id].to_numpy(dtype=float)[tested]
            correlation = correlations_from_ranks(gene_ranks, auc_values)
            p_values = p_values_for_correlations(correlation, n_cell_lines)
            fdr = adjust_fdr(p_values)
            pieces.append(
                top_genes_for_drug(
                    gene_names, correlation, p_values, fdr, drug_name, n_cell_lines
                )
            )
            elapsed = time.perf_counter() - started
            print(
                f"    [{done:3d}/{n_drugs}] {drug_name}  —  {n_cell_lines} cell lines"
                f"  —  {elapsed / 60:.1f} min",
                flush=True,
            )

    biomarkers = pd.concat(pieces, ignore_index=True)
    column_order = [
        "drug_name",
        "gene",
        "correlation",
        "p_value",
        "fdr",
        "n_cell_lines",
        "direction",
    ]
    return biomarkers[column_order], skipped


def build_drug_response(auc, catalog, name_by_id, drug_names):
    """Long table: one row per cell line and drug that has an AUC."""
    kept_ids = [compound_id for compound_id in auc.columns if name_by_id[compound_id] in set(drug_names)]
    # copy() turns the selected columns into one tidy table before reshaping.
    kept = auc.loc[:, kept_ids].copy()
    long = kept.reset_index().melt(id_vars="ModelID", var_name="compound_id", value_name="auc")
    long = long.dropna(subset=["auc"])
    long["drug_name"] = long["compound_id"].map(name_by_id)
    long = long.merge(
        catalog.reset_index()[["ModelID", "CellLineName", "OncotreeLineage"]],
        on="ModelID",
        how="left",
    )
    long = long[["ModelID", "CellLineName", "OncotreeLineage", "drug_name", "auc"]]
    return long.sort_values(["drug_name", "ModelID"]).reset_index(drop=True)


def build_expression(expression, genes):
    """Long table of expression for biomarker genes only."""
    subset = expression.loc[:, list(genes)].copy()
    long = subset.reset_index().melt(id_vars="ModelID", var_name="gene", value_name="expression")
    return long.sort_values(["gene", "ModelID"]).reset_index(drop=True)


def save_parquet(table, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(path, index=False, compression="zstd")
    size_bytes = path.stat().st_size
    size_mb = size_bytes / (1024 * 1024)
    print(f"  {path.name}: {len(table):,} rows, {size_mb:.2f} MB")
    if size_bytes >= MAX_FILE_BYTES:
        raise RuntimeError(f"{path.name} is {size_mb:.2f} MB, above the 50 MB limit.")
    return size_mb


def check_venetoclax(biomarkers, expression, auc, name_by_id):
    """Reprint one real result beside scipy, for the drug in the case study."""
    venetoclax_ids = [compound_id for compound_id, name in name_by_id.items() if name.upper() == "VENETOCLAX"]
    if len(venetoclax_ids) != 1:
        raise RuntimeError("Expected exactly one venetoclax column in the AUC file.")
    compound_id = venetoclax_ids[0]
    drug_name = name_by_id[compound_id]
    top = biomarkers.loc[
        (biomarkers["drug_name"] == drug_name) & (biomarkers["direction"] == "sensitivity")
    ].sort_values("correlation").iloc[0]
    tested = auc[compound_id].notna()
    reference_r, reference_p = stats.spearmanr(
        expression.loc[tested, top["gene"]],
        auc.loc[tested, compound_id],
    )
    if not np.isclose(top["correlation"], reference_r, atol=1e-8):
        raise RuntimeError("Venetoclax correlation does not match scipy.stats.spearmanr.")
    if not np.isclose(top["p_value"], reference_p, atol=1e-8, rtol=1e-6):
        raise RuntimeError("Venetoclax p-value does not match scipy.stats.spearmanr.")
    print()
    print(f"Check passed for {drug_name} and {top['gene']}.")
    print(f"  Spearman r = {top['correlation']:.4f}, p = {top['p_value']:.3e}, FDR = {top['fdr']:.3e}")
    print("  The same r and p come back from scipy.stats.spearmanr.")


def main():
    started = time.perf_counter()
    print("Biomarker Explorer: building the database")
    print("Spearman correlation of gene expression with drug AUC, across cell lines.")
    print("Negative correlation: higher expression, lower AUC, more sensitive.")
    print("Positive correlation: higher expression, higher AUC, more resistant.")
    print()
    check_correlation_direction()
    require_data_files()

    print("Step 1 of 4. Loading the tables.")
    catalog = load_catalog()
    print(f"  Catalog: {len(catalog):,} cell lines")
    name_by_id = load_drug_names()
    print(f"  Drug dictionary: {len(name_by_id):,} compounds")
    auc = load_auc(name_by_id)
    print(f"  AUC matrix: {auc.shape[0]:,} cell lines by {auc.shape[1]:,} drugs")
    expression = load_expression()
    print(f"  Default expression profiles: {expression.shape[0]:,} cell lines, {expression.shape[1]:,} genes")

    expression, auc, catalog = shared_cell_lines(expression, auc, catalog)
    expression, n_constant = varying_genes(expression)
    missing_lineage = int(catalog["OncotreeLineage"].isna().sum())
    print(f"  Shared cell lines, with expression and an AUC row: {len(expression):,}")
    print(f"  Genes that vary across those cell lines: {expression.shape[1]:,}")
    print(f"  Constant genes set aside: {n_constant:,}")
    print(f"  Shared cell lines with no cancer type listed: {missing_lineage:,}")

    print()
    print("Step 2 of 4. Scoring genes. This is the long step, usually several minutes.")
    biomarkers, _skipped = score_all_drugs(expression, auc, name_by_id)
    sensitivity = biomarkers.loc[biomarkers["direction"].eq("sensitivity")]
    resistance = biomarkers.loc[biomarkers["direction"].eq("resistance")]
    print()
    print(
        f"  Biomarker rows: {len(biomarkers):,} "
        f"({sensitivity['drug_name'].nunique():,} drugs, "
        f"up to {TOP_N} sensitivity and {TOP_N} resistance genes each)"
    )
    print(f"  Unique genes in that table: {biomarkers['gene'].nunique():,}")

    print()
    print("Step 3 of 4. Building the AUC table and the expression table.")
    drug_names = set(biomarkers["drug_name"])
    drug_response = build_drug_response(auc, catalog, name_by_id, drug_names)
    selected_genes = biomarkers["gene"].drop_duplicates().tolist()
    expression_long = build_expression(expression, selected_genes)
    print(f"  AUC rows: {len(drug_response):,}")
    print(f"  Expression rows: {len(expression_long):,}")

    print()
    print("Step 4 of 4. Saving parquet files.")
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    save_parquet(biomarkers, PROCESSED_DIR / "biomarkers.parquet")
    save_parquet(drug_response, PROCESSED_DIR / "drug_response.parquet")
    save_parquet(expression_long, PROCESSED_DIR / "expression.parquet")

    check_venetoclax(biomarkers, expression, auc, name_by_id)
    elapsed = time.perf_counter() - started
    print()
    print(f"Done in {elapsed / 60:.1f} minutes. Files are in {PROCESSED_DIR}.")


if __name__ == "__main__":
    main()
