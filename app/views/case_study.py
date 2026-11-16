"""Case study page: the written venetoclax interpretation."""

import base64
from pathlib import Path

import streamlit as st

from shared import ensure_data, page_header

# Path is based on this file, not the folder the app was started from.
MECHANISM_SVG = Path(__file__).resolve().parent.parent / "assets" / "bcl2_mechanism.svg"


def render():
    ensure_data()
    page_header(
        "Case study: venetoclax",
        "Does the app find the biomarker we already know about?",
    )

    st.subheader("Why I picked this drug")
    st.markdown(
        """
Venetoclax is a good test for this app because we already know the answer. It blocks BCL-2, a protein that some cancer cells rely on to avoid apoptosis, the cell's built-in self-destruct program. Take BCL-2 away and those cells finally die. That's why venetoclax is approved for chronic lymphocytic leukemia (CLL), and for acute myeloid leukemia (AML) in combination with other drugs.

So if the app is doing its job, BCL2 should come out near the top when I ask which genes predict venetoclax response, without me telling it where to look.
        """
    )

    diagram = base64.b64encode(MECHANISM_SVG.read_bytes()).decode("ascii")
    st.markdown(
        (
            '<div style="margin:1rem 0">'
            f'<img src="data:image/svg+xml;base64,{diagram}" '
            'alt="How Venetoclax works: blocking BCL-2 triggers apoptosis" '
            'style="width:100%;max-width:800px;display:block;margin:0 auto">'
            "</div>"
        ),
        unsafe_allow_html=True,
    )
    st.caption(
        "How Venetoclax works: blocking BCL-2 releases the cell's death signal and triggers apoptosis."
    )

    st.subheader("How I measured it")
    st.markdown(
        """
For each cell line I used two numbers: how much BCL2 it expresses (from DepMap RNA sequencing, as log2(TPM + 1)) and how well venetoclax killed it (from the GDSC2 drug screen, as AUC). A lower AUC means the drug worked better. Then I used a Spearman correlation to check whether those two numbers move together across cell lines.
        """
    )

    st.subheader("What I found")
    st.markdown(
        """
**BCL2 came out on top.** Across the 698 cell lines with both measurements, the correlation was r = −0.30 (FDR = 7.9 × 10⁻¹²). The negative sign is the good news: cell lines with more BCL2 tended to have a lower AUC, so venetoclax worked better on them.

An r of −0.30 is a moderate correlation. The evidence is still very strong. The false discovery rate is tiny because the same pattern holds across 698 lines after thousands of genes were tested. On the drug page the confidence label says “Very strong evidence.” That label is about the evidence, not about how large r is.

The next sensitivity genes, CRYBB1, HPS4, and MAZ, are a bit weaker (about −0.25 to −0.27).

**The link is tighter in lymphoid lines.** Keeping only lymphoid cell lines leaves 85 lines and moves the correlation to r = −0.35 (p = 0.001). The pattern is still there inside one cancer type. That fits the clinical use of venetoclax in a lymphoid leukemia. These are cell lines, not patients.

**The strongest resistance gene is less obvious.** CDC42EP2 has r = 0.25 across the same 698 lines. YAP1 is also on that list (r = 0.24). Higher expression of these genes went with a higher AUC, so a worse response. Neither is as strong as BCL2, and neither is the drug's intended target.
        """
    )

    st.subheader("Strategy memo")
    st.markdown(
        """
1. Carry **BCL2 expression** forward. It is the drug target, it is the strongest link in the screen, and the association remains when the panel is limited to lymphoid lines.
2. Start the next experiment in lymphoid models. Measure whether changing BCL2 levels changes the venetoclax response in those lines, then ask whether the same expression pattern shows up in patient samples from a lymphoid leukemia.
3. Keep CDC42EP2 and YAP1 on a watch list. They mark cell lines that responded worse, so they are reasonable resistance hypotheses. They are not yet a reason to pick a patient population. A correlation records that expression and response move together. It does not show that the gene is why the drug works.
        """
    )
