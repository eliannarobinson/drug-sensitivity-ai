"""Search by drug: two gene lists, a volcano plot, and a scatter plot."""

import pandas as pd
import streamlit as st

from shared import (
    color_key,
    drug_names,
    ensure_data,
    format_r,
    how_to_read,
    index_of,
    link_scale,
    load_biomarkers,
    metric_row,
    page_header,
    render_link_table,
    render_scatter,
    volcano_figure,
)


def _plot_label(symbol, direction):
    word = "Sensitivity" if direction == "sensitivity" else "Resistance"
    return f"{symbol} · {word}"


def _top_card(label, frame, tone):
    """The big text is the top gene, which changes with the drug."""
    if frame.empty:
        return {"label": label, "value": "None", "tone": tone}
    top = frame.iloc[0]
    return {
        "label": label,
        "value": str(top["gene_symbol"]),
        "detail": f"r = {format_r(float(top['correlation']))} · {int(top['n_cell_lines']):,} cell lines",
        "tone": tone,
    }


def render():
    ensure_data()
    page_header(
        "Search by drug",
        "Pick a drug and see which genes went with a better or worse response.",
    )
    st.caption(
        "A gene is listed when it was one of this drug's 50 strongest sensitivity genes "
        "or 50 strongest resistance genes."
    )

    drugs = drug_names()
    drug = st.selectbox("Drug", drugs, index=index_of(drugs, "VENETOCLAX"), key="drug_choice")

    biomarkers = load_biomarkers()
    drug_rows = biomarkers.loc[biomarkers["drug_name"].eq(drug)].copy()
    sensitivity = (
        drug_rows.loc[drug_rows["direction"].eq("sensitivity")]
        .sort_values(["correlation", "gene_symbol"])
        .reset_index(drop=True)
    )
    resistance = (
        drug_rows.loc[drug_rows["direction"].eq("resistance")]
        .sort_values(["correlation", "gene_symbol"], ascending=[False, True])
        .reset_index(drop=True)
    )
    if sensitivity.empty and resistance.empty:
        st.warning("This drug has no saved gene list.")
        return

    metric_row(
        [
            _top_card("Strongest sensitivity gene", sensitivity, "sens"),
            _top_card("Strongest resistance gene", resistance, "res"),
        ]
    )
    how_to_read("Gene")
    color_key()

    scale = link_scale([sensitivity, resistance])
    left, right = st.columns(2)
    with left:
        render_link_table(
            sensitivity,
            "gene_symbol",
            "sensitivity",
            scale,
            "No gene in the saved list points toward sensitivity.",
        )
    with right:
        render_link_table(
            resistance,
            "gene_symbol",
            "resistance",
            scale,
            "No gene in the saved list points toward resistance.",
        )

    st.plotly_chart(
        volcano_figure(drug_rows, drug),
        width="stretch",
        theme=None,
        alt=f"Volcano plot of gene correlations versus confidence for {drug}.",
    )
    st.caption(
        "Each point is one of the 100 genes saved for this drug: the 50 strongest sensitivity links "
        "and the 50 strongest resistance links. Blue is sensitivity and pink is resistance. "
        "Higher on the plot means a smaller false discovery rate, so more confidence. "
        "The strongest names that have room are labeled. Hover any point for the rest. "
        "The other genes tested for this drug are not on this plot."
    )

    ranked = pd.concat([sensitivity, resistance], ignore_index=True)
    labels = [_plot_label(row.gene_symbol, row.direction) for row in ranked.itertuples(index=False)]
    gene_for_label = dict(zip(labels, ranked["gene"]))
    if st.session_state.get("active_drug") != drug or st.session_state.get("drug_plot_label") not in gene_for_label:
        st.session_state["active_drug"] = drug
        st.session_state["drug_plot_label"] = labels[0]

    st.divider()
    choice = st.selectbox(
        "Show the scatter plot for",
        labels,
        key="drug_plot_label",
    )
    stored = ranked.loc[ranked["gene"].eq(gene_for_label[choice])].iloc[0]
    render_scatter(
        drug,
        stored["gene"],
        stored["gene_symbol"],
        key_prefix=f"drug_{drug}_{stored['gene']}",
        stored_row=stored,
    )
