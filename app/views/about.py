"""Where the numbers come from, and what they do not mean."""

import streamlit as st

from shared import ensure_data, page_header


def render():
    ensure_data()
    page_header(
        "About and limitations",
        "Where the numbers come from, and what a correlation does not prove.",
    )

    st.subheader("Data sources")
    st.markdown(
        """
**DepMap models.** `Model.csv` supplies the cell-line name and the cancer type (`OncotreeLineage`).

**DepMap expression.** Expression is log2(TPM + 1) for human protein-coding genes.
TPM is transcripts per million, a measure of how much RNA a gene produced.
When a cell line was sequenced more than once, the app keeps the row marked
`IsDefaultEntryForModel = Yes`, so each cell line contributes one profile.

**GDSC2 drug response.** GDSC2 is the Genomics of Drug Sensitivity in Cancer screen, version 2.
The sensitivity measure is AUC, the area under the dose-response curve.
A lower AUC means the drug reduced survival more.
The matrix columns are compound codes such as `DPC-006825`.
The GDSC2 conditions file supplies the drug name for each code. Venetoclax is `DPC-006825`.
        """
    )

    st.subheader("Methods")
    st.markdown(
        """
The comparison uses the 708 cell lines that have a default expression profile, a GDSC2 row, and a catalog entry.
Every drug in this release clears 100 measured cell lines in that overlap, so all 286 drugs are included.

For each drug, the script ranks each gene's expression across those cell lines and ranks the AUC values the same way.
The Spearman correlation is the correlation of the two ranks.
Ranks keep a few extreme values from dominating the result.
A **negative** correlation means higher expression tends to sit with a lower AUC, which is greater sensitivity.
A **positive** correlation means higher expression tends to sit with a higher AUC, which is resistance.

Each correlation has a p-value. Thousands of genes are tested for every drug, so some small p-values would appear by chance.
The Benjamini-Hochberg procedure turns those p-values into a false discovery rate (FDR), separately for each drug, using every gene tested.
The saved lists are then the 50 most negative genes and the 50 most positive genes.
The FDR on a row was computed before that cutoff.

Gene labels in the files look like `BCL2 (596)`. The number is an Entrez gene ID. The app shows the symbol, BCL2.
        """
    )

    st.subheader("Limitations")
    st.markdown(
        """
**Cell lines are models, and patients are not in this dataset.**
A cell line is cancer cells grown in a dish. It has no immune system, no neighboring tissue, and often a long history in culture.
A pattern here is a hypothesis to take into a better model. It is not evidence about how a drug will behave in a person.

**A correlation records that two measurements move together.**
The gene can track drug response without being the reason the drug works. It may travel with some other difference between the cell lines.

**Cancer type can create the correlation.**
Cell lines from different tissues differ in thousands of genes, and they can differ in drug response too.
A correlation across the whole panel can mean that the gene marks a cancer type which responds, rather than that the gene tracks response inside one cancer type.
The filter on Explore and on the drug page recalculates the Spearman correlation for the cancer types you leave on the plot.

**A few practical limits of this app.**
The gene search only knows genes that reached a top-50 list for at least one drug.
A blank AUC means that cell line was not tested with that drug, so it is left out of the correlation. It is not entered as zero.
The expression value is one default profile per cell line, not a measurement from a tumor sample.
        """
    )
