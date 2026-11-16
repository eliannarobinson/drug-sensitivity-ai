"""Search by gene: which drugs include this symbol in a top-50 list."""

import streamlit as st

from shared import (
    color_key,
    ensure_data,
    how_to_read,
    link_scale,
    load_biomarkers,
    metric_row,
    page_header,
    render_link_table,
    render_scatter,
)


def render():
    ensure_data()
    page_header(
        "Search by gene",
        "Type a symbol such as BCL2 to see which drugs that gene is linked to.",
    )
    st.caption(
        "A drug is listed when this gene was one of its 50 strongest sensitivity genes "
        "or 50 strongest resistance genes. Capital letters do not matter."
    )

    query = st.text_input(
        "Gene",
        value="BCL2",
        key="gene_query",
        placeholder="Type a gene, such as BCL2",
    ).strip()
    if not query:
        st.info("Type a gene, such as BCL2 or SLFN11.")
        return

    biomarkers = load_biomarkers()
    symbols = biomarkers["gene_symbol"].str.upper()
    matched = biomarkers.loc[symbols.eq(query.upper())].copy()
    if matched.empty:
        st.warning(f"No saved result uses the symbol {query}.")
        suggestions = sorted(
            biomarkers.loc[symbols.str.startswith(query.upper()), "gene_symbol"].drop_duplicates(),
            key=str.upper,
        )[:25]
        if suggestions:
            st.write("Symbols that start with those letters: " + ", ".join(suggestions))
        else:
            st.caption(
                "The search covers genes that reached a top-50 list for at least one drug, not every gene in the genome."
            )
        return

    symbol = matched["gene_symbol"].iloc[0]
    gene = matched["gene"].iloc[0]
    sensitivity = (
        matched.loc[matched["direction"].eq("sensitivity")]
        .sort_values(["correlation", "drug_name"])
        .reset_index(drop=True)
    )
    resistance = (
        matched.loc[matched["direction"].eq("resistance")]
        .sort_values(["correlation", "drug_name"], ascending=[False, True])
        .reset_index(drop=True)
    )

    metric_row(
        [
            {
                "label": "Linked to sensitivity",
                "value": f"{len(sensitivity)} drugs",
                "tone": "sens",
            },
            {
                "label": "Linked to resistance",
                "value": f"{len(resistance)} drugs",
                "tone": "res",
            },
        ]
    )
    how_to_read("Drug")
    color_key()

    scale = link_scale([sensitivity, resistance])
    left, right = st.columns(2)
    with left:
        render_link_table(
            sensitivity,
            "drug_name",
            "sensitivity",
            scale,
            "This gene is not on a sensitivity list.",
        )
    with right:
        render_link_table(
            resistance,
            "drug_name",
            "resistance",
            scale,
            "This gene is not on a resistance list.",
        )

    drug_options = sorted(matched["drug_name"].unique(), key=str.upper)
    default_drug = "VENETOCLAX" if "VENETOCLAX" in drug_options else drug_options[0]
    st.divider()
    plotted_drug = st.selectbox(
        "Plot this gene against a drug from the list",
        drug_options,
        index=drug_options.index(default_drug),
        key=f"gene_plot_{symbol}",
    )
    stored = matched.loc[matched["drug_name"].eq(plotted_drug)].iloc[0]
    render_scatter(
        plotted_drug,
        gene,
        symbol,
        key_prefix=f"gene_{symbol}_{plotted_drug}",
        stored_row=stored,
    )
