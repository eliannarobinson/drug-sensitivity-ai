"""Explore: any saved gene against any drug, with a live correlation."""

import streamlit as st

from shared import (
    drug_names,
    ensure_data,
    gene_symbols,
    index_of,
    load_biomarkers,
    page_header,
    render_scatter,
)


def render():
    ensure_data()
    page_header(
        "Explore",
        "Draw one gene against one drug, and see whether the pattern holds in a single cancer type.",
    )

    drugs = drug_names()
    symbols = gene_symbols()
    left, right = st.columns(2)
    with left:
        drug = st.selectbox("Drug", drugs, index=index_of(drugs, "VENETOCLAX"), key="explore_drug")
    with right:
        symbol = st.selectbox(
            "Gene",
            symbols,
            index=index_of(symbols, "BCL2"),
            key="explore_gene",
            help="Type a symbol to search the list.",
        )
    st.caption("The gene menu lists genes that were a top sensitivity or resistance gene for at least one drug.")

    biomarkers = load_biomarkers()
    gene_rows = biomarkers.loc[biomarkers["gene_symbol"].eq(symbol)]
    gene = gene_rows["gene"].iloc[0]
    stored = gene_rows.loc[gene_rows["drug_name"].eq(drug)]
    stored_row = None if stored.empty else stored.iloc[0]
    if stored_row is None:
        st.info(
            f"{symbol} is outside the top 50 sensitivity and resistance genes for {drug}. "
            "The plot still includes every cell line measured for both."
        )

    render_scatter(
        drug,
        gene,
        symbol,
        key_prefix=f"explore_{drug}_{gene}",
        stored_row=stored_row,
    )
