"""Home page: a short introduction for someone opening the app."""

import streamlit as st

from shared import color_key, ensure_data, load_biomarkers, load_drug_response, metric_row, page_header


def render():
    ensure_data()
    biomarkers = load_biomarkers()
    n_drugs = biomarkers["drug_name"].nunique()
    n_genes = biomarkers["gene_symbol"].nunique()
    n_lines = load_drug_response()["ModelID"].nunique()

    page_header(
        "Biomarker Explorer",
        "Find genes that track how cancer cell lines respond to a drug.",
    )
    st.html(
        '<p class="bx-lead">Search a drug or a gene to see the strongest links in the DepMap and GDSC2 '
        "cell-line screens. Sky blue means higher expression went with a better response. "
        "Pink means higher expression went with a worse one. A clear first example is Venetoclax and BCL2.</p>"
    )
    metric_row(
        [
            {"label": "Drugs", "value": f"{n_drugs}"},
            {"label": "Genes in the lists", "value": f"{n_genes:,}"},
            {"label": "Cell lines", "value": f"{n_lines:,}"},
        ]
    )

    st.html(
        """
<div class="bx-cards">
  <div class="bx-card">
    <h3>Lower AUC</h3>
    <p>The drug reduced survival more. Those cells were more sensitive.</p>
  </div>
  <div class="bx-card sens">
    <h3>Sky blue: sensitivity</h3>
    <p>Higher gene expression went with a better response. The correlation is negative.</p>
  </div>
  <div class="bx-card res">
    <h3>Pink: resistance</h3>
    <p>Higher gene expression went with a worse response. The correlation is positive.</p>
  </div>
</div>
"""
    )
    color_key()

    st.html(
        """
<h2 class="bx-block">Where to start</h2>
<ol class="bx-start">
  <li><b>Search by drug.</b> <span>Pick a drug and see the genes most linked to sensitivity or resistance.</span></li>
  <li><b>Search by gene.</b> <span>Type a symbol such as BCL2 and see which drugs that gene is linked to.</span></li>
  <li><b>Explore.</b> <span>Draw expression against drug response, color the points by cancer type, and filter.</span></li>
  <li><b>Case study: Venetoclax.</b> <span>A written read of the BCL2 result and what to test next.</span></li>
</ol>
"""
    )
    st.caption("A useful first look is Venetoclax and BCL2. Choose a page on the left.")
