"""Load the parquet tables and draw the pieces used by more than one page."""

import html
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

from colors import BACKGROUND, MUTED, PANEL, RESISTANT, SENSITIVE, TEXT, apply_dark_layout

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "processed"
BIOMARKERS_PATH = PROCESSED / "biomarkers.parquet"
DRUG_RESPONSE_PATH = PROCESSED / "drug_response.parquet"
EXPRESSION_PATH = PROCESSED / "expression.parquet"

# Sky blue always means sensitivity. Pink always means resistance.
SENSITIVITY = SENSITIVE
RESISTANCE = RESISTANT

# Cancer-type colors stay off that blue and that orange, so a point's color
# never looks like a direction.
CANCER_COLORS = [
    "#009E73",
    "#CC79A7",
    "#882255",
    "#44AA99",
    "#332288",
    "#C4B5FD",
    "#117733",
    "#AA4499",
    "#999933",
    "#661100",
    "#6699CC",
    "#888888",
    "#CC6677",
    "#DDCC77",
    "#225555",
    "#994455",
    "#66AA77",
    "#445588",
    "#AA7744",
    "#775577",
    "#337777",
    "#774433",
    "#556677",
    "#99AA66",
    "#663355",
    "#228866",
    "#886688",
    "#446655",
    "#777744",
    "#555577",
]

DIRECTION_TITLE = {
    "sensitivity": "Linked to sensitivity (cells responded better)",
    "resistance": "Linked to resistance (cells responded worse)",
}

DESIGN_CSS = f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Source+Sans+3:wght@400;600;700&display=swap');
html, body, [class*="css"], .stApp {{
  font-family: "Source Sans 3", "Segoe UI", sans-serif;
}}
.stApp, [data-testid="stAppViewContainer"], [data-testid="stMain"] {{
  background: {BACKGROUND};
  color: {TEXT};
}}
header[data-testid="stHeader"] {{ background: rgba(11, 18, 32, 0.94); }}
div[data-testid="stDecoration"] {{ display: none; }}
footer {{ visibility: hidden; height: 0; }}
section[data-testid="stSidebar"] {{
  background: {PANEL};
  border-right: 1px solid #24344F;
}}
[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p,
[data-testid="stSidebar"] span,
[data-testid="stSidebar"] a {{
  color: {TEXT};
}}
label[data-testid="stWidgetLabel"],
label[data-testid="stWidgetLabel"] p,
label[data-testid="stWidgetLabel"] span {{
  color: {TEXT} !important;
}}
div[data-testid="stCaptionContainer"],
div[data-testid="stCaptionContainer"] p {{
  color: #C5D2E8 !important;
}}
div[data-testid="stMarkdownContainer"] p,
.stApp li {{
  color: {TEXT};
}}
div[data-baseweb="select"] > div {{
  background-color: {PANEL} !important;
  color: {TEXT} !important;
  border-color: #2A3A56 !important;
}}
div[data-baseweb="select"] span,
div[data-baseweb="select"] input {{
  color: {TEXT} !important;
}}
.st-key-gene_query label p,
.st-key-drug_choice label p {{
  font-size: 1.25rem !important;
  font-weight: 700 !important;
}}
.st-key-gene_query input,
.st-key-drug_choice div[data-baseweb="select"] > div {{
  font-size: 1.15rem !important;
  min-height: 3rem;
  border: 1px solid {SENSITIVE} !important;
  border-radius: 12px !important;
}}
.st-key-gene_query input {{
  background: {PANEL} !important;
  color: {TEXT} !important;
  padding: 0.4rem 0.85rem !important;
}}
div[data-testid="stMetric"] {{
  background: {PANEL};
  border: 1px solid #2A3A56;
  border-top: 3px solid {SENSITIVE};
  border-radius: 12px;
  padding: 0.7rem 0.95rem 0.45rem;
}}
div[data-testid="stMetricLabel"],
div[data-testid="stMetricLabel"] p,
div[data-testid="stMetricLabel"] label {{
  color: #C5D2E8 !important;
  font-weight: 600;
}}
div[data-testid="stMetricValue"] {{ color: {TEXT} !important; }}
div[data-testid="stAlert"] {{ border-radius: 12px; }}
hr {{ border-color: #24344F; }}
.bx-head {{ margin: 0.15rem 0 1rem; }}
.bx-title {{
  font-size: 2rem;
  line-height: 1.15;
  font-weight: 700;
  letter-spacing: -0.03em;
  color: {TEXT};
  margin: 0;
  padding: 0;
}}
.bx-head p {{ margin: 0.35rem 0 0; color: #C5D2E8; font-size: 1.05rem; }}
.bx-lead {{
  font-size: 1.08rem;
  line-height: 1.55;
  color: {TEXT};
  background: {PANEL};
  border-left: 4px solid {SENSITIVE};
  border-radius: 0 12px 12px 0;
  padding: 0.9rem 1.1rem;
  margin: 0.2rem 0 0.9rem;
}}
.bx-metrics {{ display: flex; gap: 12px; margin: 0.2rem 0 1rem; }}
.bx-metric {{
  flex: 1;
  background: {PANEL};
  border: 1px solid #2A3A56;
  border-top: 3px solid {SENSITIVE};
  border-radius: 12px;
  padding: 0.75rem 0.95rem 0.7rem;
}}
.bx-metric.sens {{ border-top-color: {SENSITIVE}; }}
.bx-metric.res {{ border-top-color: {RESISTANT}; }}
.bx-metric .k {{
  font-size: 0.75rem;
  letter-spacing: 0.04em;
  text-transform: uppercase;
  color: #C5D2E8;
  font-weight: 600;
}}
.bx-metric .v {{
  font-size: 1.45rem;
  font-weight: 700;
  color: {TEXT};
  margin-top: 0.15rem;
  letter-spacing: -0.02em;
}}
.bx-metric .d {{ color: #C5D2E8; font-size: 0.92rem; margin-top: 0.1rem; }}
.bx-howto {{
  border: 1px solid #2A3A56;
  border-radius: 12px;
  padding: 0.85rem 1rem 0.95rem;
  margin: 0.2rem 0 0.9rem;
  background: {PANEL};
}}
.bx-howto-title {{
  margin: 0 0 0.55rem;
  font-size: 1rem;
  font-weight: 700;
  letter-spacing: -0.01em;
  color: {TEXT};
}}
.bx-howto-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 0.7rem 1.2rem; }}
.bx-howto-grid b {{ display: block; color: {TEXT}; margin-bottom: 0.12rem; }}
.bx-howto-grid span, .bx-note {{ color: #C5D2E8; font-size: 0.92rem; line-height: 1.4; }}
.bx-note {{ display: block; margin-top: 0.7rem; }}
.bx-key {{ display: flex; gap: 1.2rem; margin: 0.1rem 0 0.8rem; color: {TEXT}; font-size: 0.95rem; }}
.bx-key i {{
  display: inline-block;
  width: 0.7rem;
  height: 0.7rem;
  border-radius: 99px;
  margin-right: 0.4rem;
}}
.bx-key .sens i {{ background: {SENSITIVE}; }}
.bx-key .res i {{ background: {RESISTANT}; }}
.bx-table-card {{
  border: 1px solid #2A3A56;
  border-radius: 12px;
  padding: 0.75rem 0.85rem 0.4rem;
  margin-bottom: 0.6rem;
  background: {PANEL};
}}
.bx-block {{
  display: flex;
  align-items: center;
  gap: 0.45rem;
  margin: 0 0 0.55rem;
  font-size: 1.02rem;
  line-height: 1.3;
  color: {TEXT};
}}
.bx-block.sens {{ color: {SENSITIVE}; }}
.bx-block.res {{ color: {RESISTANT}; }}
.bx-block .swatch {{ width: 0.7rem; height: 0.7rem; border-radius: 99px; flex: none; }}
.bx-block.sens .swatch {{ background: {SENSITIVE}; }}
.bx-block.res .swatch {{ background: {RESISTANT}; }}
.bx-scroll {{ max-height: 460px; overflow: auto; }}
table.bx {{ width: 100%; border-collapse: collapse; font-size: 0.95rem; }}
table.bx th {{
  position: sticky;
  top: 0;
  background: {PANEL};
  text-align: left;
  font-size: 0.72rem;
  letter-spacing: 0.04em;
  text-transform: uppercase;
  color: #C5D2E8;
  font-weight: 600;
  padding: 0.35rem 0.35rem 0.45rem;
  border-bottom: 1px solid #2A3A56;
}}
table.bx td {{
  padding: 0.42rem 0.35rem;
  border-bottom: 1px solid #24344F;
  color: {TEXT};
  vertical-align: middle;
}}
table.bx td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
.bx-strength {{ display: flex; align-items: center; gap: 8px; min-width: 150px; }}
.bx-track {{ flex: 1; height: 8px; background: #24344F; border-radius: 99px; overflow: hidden; }}
.bx-fill {{ height: 8px; border-radius: 99px; }}
.bx-fill.sens {{ background: {SENSITIVE}; }}
.bx-fill.res {{ background: {RESISTANT}; }}
.bx-rval {{ font-variant-numeric: tabular-nums; font-weight: 600; min-width: 3.1rem; color: {TEXT}; }}
.bx-badge {{
  display: inline-block;
  border-radius: 99px;
  padding: 0.12rem 0.55rem;
  font-size: 0.8rem;
  font-weight: 600;
  cursor: help;
  white-space: nowrap;
}}
.bx-badge.vs {{ background: {SENSITIVE}; color: {BACKGROUND}; }}
.bx-badge.st {{ background: #1B3A52; color: {TEXT}; }}
.bx-badge.md {{ background: #24344F; color: {TEXT}; }}
.bx-badge.sg {{ background: transparent; color: #C5D2E8; border: 1px solid #2A3A56; }}
.bx-empty {{ color: #C5D2E8; margin: 0.2rem 0 0.6rem; }}
.bx-pill {{
  display: inline-flex;
  align-items: center;
  gap: 0.4rem;
  border-radius: 99px;
  padding: 0.28rem 0.7rem;
  font-weight: 600;
  font-size: 0.92rem;
  margin: 0.3rem 0 0.7rem;
}}
.bx-pill .dot {{ width: 0.65rem; height: 0.65rem; border-radius: 99px; }}
.bx-pill.sens {{ background: #123044; color: {SENSITIVE}; }}
.bx-pill.sens .dot {{ background: {SENSITIVE}; }}
.bx-pill.res {{ background: #3A2030; color: {RESISTANT}; }}
.bx-pill.res .dot {{ background: {RESISTANT}; }}
.bx-cards {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; margin: 0.4rem 0 1.1rem; }}
.bx-card {{
  border: 1px solid #2A3A56;
  border-top: 3px solid {SENSITIVE};
  border-radius: 12px;
  padding: 0.85rem 1rem 0.95rem;
  background: {PANEL};
}}
.bx-card.sens {{ border-top-color: {SENSITIVE}; }}
.bx-card.res {{ border-top-color: {RESISTANT}; }}
.bx-card h3 {{ margin: 0 0 0.3rem; font-size: 1rem; color: {TEXT}; }}
.bx-card.sens h3 {{ color: {SENSITIVE}; }}
.bx-card.res h3 {{ color: {RESISTANT}; }}
.bx-card p {{ margin: 0; color: #C5D2E8; line-height: 1.45; }}
.bx-start {{ margin: 0.2rem 0 0.8rem; padding-left: 1.2rem; color: {TEXT}; }}
.bx-start li {{ margin: 0.35rem 0; line-height: 1.45; }}
.bx-start span {{ color: #C5D2E8; }}
</style>
"""


def apply_design():
    """Load the shared type, spacing, and color styles."""
    st.html(DESIGN_CSS)


def page_header(title: str, subtitle: str):
    apply_design()
    st.html(
        f'<div class="bx-head"><div class="bx-title">{html.escape(title)}</div>'
        f"<p>{html.escape(subtitle)}</p></div>"
    )


def ensure_data():
    """Stop the page with a plain message if the database has not been built."""
    missing = [
        path.name
        for path in (BIOMARKERS_PATH, DRUG_RESPONSE_PATH, EXPRESSION_PATH)
        if not path.exists()
    ]
    if missing:
        apply_design()
        st.error(
            "Missing data file"
            + ("s" if len(missing) > 1 else "")
            + ": "
            + ", ".join(missing)
            + ". From the project folder, run python scripts/build_database.py, then start the app again."
        )
        st.stop()


def _file_time(path: Path) -> float:
    return path.stat().st_mtime


@st.cache_data(show_spinner="Loading biomarker results...")
def _load_biomarkers(_file_time: float) -> pd.DataFrame:
    # gene looks like "BCL2 (596)": symbol, then the Entrez number.
    biomarkers = pd.read_parquet(BIOMARKERS_PATH)
    parts = biomarkers["gene"].str.extract(r"^(?P<gene_symbol>.+?) \((?P<entrez>\d+)\)$")
    biomarkers["gene_symbol"] = parts["gene_symbol"]
    biomarkers["entrez"] = parts["entrez"]
    return biomarkers


@st.cache_data(show_spinner="Loading drug response...")
def _load_drug_response(_file_time: float) -> pd.DataFrame:
    response = pd.read_parquet(DRUG_RESPONSE_PATH)
    response["OncotreeLineage"] = response["OncotreeLineage"].fillna("Not listed")
    return response


@st.cache_data(show_spinner="Loading gene expression. This happens once.")
def _load_expression(_file_time: float) -> pd.DataFrame:
    return pd.read_parquet(EXPRESSION_PATH)


def load_biomarkers() -> pd.DataFrame:
    return _load_biomarkers(_file_time(BIOMARKERS_PATH))


def load_drug_response() -> pd.DataFrame:
    return _load_drug_response(_file_time(DRUG_RESPONSE_PATH))


def load_expression() -> pd.DataFrame:
    return _load_expression(_file_time(EXPRESSION_PATH))


def drug_names() -> list[str]:
    names = load_biomarkers()["drug_name"].unique().tolist()
    return sorted(names, key=str.upper)


def gene_symbols() -> list[str]:
    symbols = load_biomarkers()["gene_symbol"].unique().tolist()
    return sorted(symbols, key=str.upper)


def index_of(options: list[str], value: str) -> int:
    try:
        return options.index(value)
    except ValueError:
        return 0


def prose_name(name: str) -> str:
    """VENETOCLAX becomes venetoclax. Mixed codes stay as written."""
    if re.fullmatch(r"[A-Z0-9][A-Z0-9 .+\-()/]*", name) and re.search(r"[A-Z]", name):
        return name.lower()
    return name


def format_r(value: float, digits: int = 2) -> str:
    """A correlation with a true minus sign, so it does not look like a hyphen."""
    sign = "−" if value < 0 else ""
    return f"{sign}{abs(value):.{digits}f}"


def confidence_label(fdr: float) -> str:
    """Label how strong the evidence is, from the false discovery rate. This is not the size of the correlation."""
    if pd.isna(fdr):
        return "Not available"
    if fdr < 0.001:
        return "Very strong evidence"
    if fdr < 0.01:
        return "Strong evidence"
    if fdr < 0.05:
        return "Moderate evidence"
    return "Suggestive evidence"


def link_scale(frames) -> float:
    """The strongest |correlation| on the page, so bar lengths can be compared."""
    peaks = []
    for frame in frames:
        if frame is not None and not frame.empty:
            peaks.append(float(frame["correlation"].abs().max()))
    if not peaks or max(peaks) == 0:
        return 1.0
    return max(peaks)


def metric_row(cards: list[dict]):
    pieces = []
    for card in cards:
        tone = card.get("tone", "")
        detail = card.get("detail")
        detail_html = f'<div class="d">{html.escape(detail)}</div>' if detail else ""
        pieces.append(
            f'<div class="bx-metric {tone}">'
            f'<div class="k">{html.escape(card["label"])}</div>'
            f'<div class="v">{html.escape(str(card["value"]))}</div>'
            f"{detail_html}</div>"
        )
    st.html(f'<div class="bx-metrics">{"".join(pieces)}</div>')


def color_key():
    st.html(
        '<div class="bx-key">'
        '<span class="sens"><i></i>Sky blue: sensitivity, cells responded better</span>'
        '<span class="res"><i></i>Pink: resistance, cells responded worse</span>'
        "</div>"
    )


def how_to_read(first_column: str):
    """Plain-language key for the four table columns."""
    if first_column == "Gene":
        first = "The short name of the gene, such as BCL2."
    else:
        first = "The drug this gene is linked to."
    st.html(
        f"""
<div class="bx-howto">
  <div class="bx-howto-title">How to read this</div>
  <div class="bx-howto-grid">
    <div><b>{html.escape(first_column)}</b><span>{first}</span></div>
    <div><b>Strength of link</b><span>How tightly expression and drug response move together. A longer bar is a stronger link. The number is the correlation. A minus sign means higher expression went with a better response.</span></div>
    <div><b>Confidence</b><span>How strong the evidence is that the link is not a fluke after testing thousands of genes. This is separate from how large the correlation is. Hover a label to see the false discovery rate. Very strong evidence is below 0.001, Strong evidence is below 0.01, Moderate evidence is below 0.05, and Suggestive evidence is 0.05 or higher.</span></div>
    <div><b>Cell lines tested</b><span>How many cell lines had both a gene measurement and a drug measurement.</span></div>
  </div>
  <span class="bx-note">In the resistance table the correlation is positive: higher expression went with a worse response. Bar length is compared with the strongest link in the lists on this page.</span>
</div>
"""
    )


def _badge(fdr: float) -> str:
    label = confidence_label(fdr)
    slug = {
        "Very strong evidence": "vs",
        "Strong evidence": "st",
        "Moderate evidence": "md",
        "Suggestive evidence": "sg",
    }.get(label, "sg")
    if pd.isna(fdr):
        hover = "No false discovery rate was saved for this row."
    else:
        hover = (
            f"False discovery rate = {fdr:.2e}. "
            "A smaller number means more confidence that the link is not a chance result."
        )
    return f'<span class="bx-badge {slug}" title="{html.escape(hover)}">{label}</span>'


def render_link_table(frame: pd.DataFrame, name_column: str, direction: str, scale: float, empty_message: str):
    """A results table with a colored strength bar and a confidence label."""
    tone = "sens" if direction == "sensitivity" else "res"
    title = DIRECTION_TITLE[direction]
    name_header = "Gene" if name_column == "gene_symbol" else "Drug"
    if frame.empty:
        st.html(
            f'<div class="bx-table-card"><div class="bx-block {tone}"><span class="swatch"></span>'
            f"{html.escape(title)}</div><p class=\"bx-empty\">{html.escape(empty_message)}</p></div>"
        )
        return

    rows = []
    for record in frame.to_dict("records"):
        name = html.escape(str(record[name_column]))
        correlation = float(record["correlation"])
        width = 0.0 if scale <= 0 else min(100.0, abs(correlation) / scale * 100.0)
        exact = f"Correlation = {correlation:.3f}"
        fdr = record["fdr"]
        fdr_value = float(fdr) if pd.notna(fdr) else np.nan
        rows.append(
            "<tr>"
            f"<td>{name}</td>"
            '<td><div class="bx-strength">'
            f'<div class="bx-track"><div class="bx-fill {tone}" style="width:{width:.1f}%"></div></div>'
            f'<span class="bx-rval" title="{html.escape(exact)}">{format_r(correlation)}</span>'
            "</div></td>"
            f"<td>{_badge(fdr_value)}</td>"
            f'<td class="num">{int(record["n_cell_lines"]):,}</td>'
            "</tr>"
        )
    st.html(
        f'<div class="bx-table-card"><div class="bx-block {tone}"><span class="swatch"></span>'
        f"{html.escape(title)}</div>"
        '<div class="bx-scroll"><table class="bx"><thead><tr>'
        f"<th>{name_header}</th><th>Strength of link</th><th>Confidence</th>"
        '<th style="text-align:right">Cell lines tested</th>'
        "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></div></div>"
    )


def volcano_figure(frame: pd.DataFrame, drug_name: str):
    """Correlation against confidence for the genes saved for one drug."""
    plotted = frame.loc[frame["direction"].isin(["sensitivity", "resistance"])].copy()
    plotted["fdr_plot"] = plotted["fdr"].clip(lower=1e-300)
    plotted["confidence_y"] = -np.log10(plotted["fdr_plot"])
    plotted["confidence_word"] = plotted["fdr"].map(confidence_label)

    figure = go.Figure()
    label_rows = []
    for direction, color, ascending in (
        ("sensitivity", SENSITIVITY, True),
        ("resistance", RESISTANCE, False),
    ):
        part = plotted.loc[plotted["direction"].eq(direction)]
        if part.empty:
            continue
        figure.add_trace(
            go.Scatter(
                x=part["correlation"],
                y=part["confidence_y"],
                mode="markers",
                marker={"size": 11, "color": color, "line": {"width": 1, "color": "white"}},
                name="Sensitivity" if direction == "sensitivity" else "Resistance",
                customdata=list(
                    zip(
                        part["gene_symbol"].tolist(),
                        part["fdr"].astype(float).tolist(),
                        part["confidence_word"].tolist(),
                        part["n_cell_lines"].astype(int).tolist(),
                    )
                ),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Correlation: %{x:.3f}<br>"
                    "FDR: %{customdata[1]:.2e}<br>"
                    "Confidence: %{customdata[2]}<br>"
                    "Cell lines: %{customdata[3]}<extra></extra>"
                ),
            )
        )
        # Keep the strongest names that have room to be read. Nearby genes stay on hover.
        chosen = part.sort_values(["correlation", "gene_symbol"], ascending=[ascending, True]).head(5)
        chosen = chosen.sort_values("confidence_y", ascending=False)
        last_y = None
        for row in chosen.itertuples(index=False):
            label_y = float(row.confidence_y)
            if last_y is not None and last_y - label_y < 0.9:
                continue
            label_rows.append((direction, row, label_y))
            last_y = label_y

    span = float(plotted["correlation"].max() - plotted["correlation"].min())
    pad = max(0.1, span * 0.28)
    ymax = float(plotted["confidence_y"].max())
    figure.add_hline(
        y=-np.log10(0.05),
        line_dash="dot",
        line_color=MUTED,
        annotation_text="FDR 0.05",
        annotation_position="bottom right",
        annotation_font_color=MUTED,
    )
    figure.add_vline(x=0, line_color=MUTED, line_width=1)
    for direction, row, label_y in label_rows:
        figure.add_annotation(
            x=float(row.correlation),
            y=label_y,
            text=str(row.gene_symbol),
            showarrow=False,
            xanchor="right" if direction == "sensitivity" else "left",
            xshift=-8 if direction == "sensitivity" else 8,
            font={"size": 12, "color": TEXT, "family": "Source Sans 3, sans-serif"},
        )
    figure.update_layout(
        title={"text": f"Links for {drug_name}", "font": {"size": 20, "color": TEXT}},
        height=560,
        margin={"l": 20, "r": 24, "t": 70, "b": 20},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "x": 0, "title": ""},
    )
    figure.update_xaxes(
        title="Correlation (left: sensitivity, right: resistance)",
        range=[float(plotted["correlation"].min()) - pad, float(plotted["correlation"].max()) + pad],
        automargin=True,
    )
    figure.update_yaxes(
        title="Confidence (−log10 FDR)",
        range=[-0.3, ymax * 1.12 + 0.5],
        automargin=True,
    )
    return apply_dark_layout(figure)


def points_for(drug_name: str, gene: str) -> pd.DataFrame:
    """One row per cell line that has both this gene's expression and this drug's AUC."""
    expression = load_expression()
    response = load_drug_response()
    gene_values = expression.loc[expression["gene"].eq(gene), ["ModelID", "expression"]]
    drug_values = response.loc[
        response["drug_name"].eq(drug_name),
        ["ModelID", "CellLineName", "OncotreeLineage", "auc"],
    ]
    merged = gene_values.merge(drug_values, on="ModelID", how="inner")
    return merged.drop_duplicates(subset="ModelID")


def correlation_of(points: pd.DataFrame):
    """Spearman correlation for the points on the plot. None when it is undefined."""
    if len(points) < 3:
        return None
    if points["expression"].nunique(dropna=True) < 2 or points["auc"].nunique(dropna=True) < 2:
        return None
    coefficient, p_value = stats.spearmanr(points["expression"], points["auc"])
    if pd.isna(coefficient) or pd.isna(p_value):
        return None
    return float(coefficient), float(p_value)


def format_p(p_value: float) -> str:
    if p_value < 0.001:
        return f"{p_value:.2e}"
    return f"{p_value:.3f}"


def direction_pill(correlation: float):
    """Color cue for the scatter: sky blue for a better response, pink for a worse one."""
    if correlation < 0:
        tone = "sens"
        word = "Sensitivity"
        detail = "higher expression went with a better response"
    elif correlation > 0:
        tone = "res"
        word = "Resistance"
        detail = "higher expression went with a worse response"
    else:
        return
    st.html(
        f'<div class="bx-pill {tone}"><span class="dot"></span>'
        f"{word}: {detail} (r = {format_r(correlation, digits=3)}).</div>"
    )


def scatter_figure(points: pd.DataFrame, drug_name: str, gene_symbol: str):
    cancer_types = sorted(points["OncotreeLineage"].unique())
    figure = px.scatter(
        points,
        x="expression",
        y="auc",
        color="OncotreeLineage",
        hover_name="CellLineName",
        category_orders={"OncotreeLineage": cancer_types},
        color_discrete_sequence=CANCER_COLORS,
        labels={
            "expression": f"{gene_symbol} expression, log2(TPM + 1)",
            "auc": "AUC (lower = more sensitive)",
            "OncotreeLineage": "Cancer type",
        },
        template="plotly_white",
    )
    figure.update_traces(marker={"size": 9, "opacity": 0.85, "line": {"width": 0}})
    figure.update_layout(
        title={"text": f"{gene_symbol} and {drug_name}", "font": {"size": 20}},
        height=640,
        legend_title_text="Cancer type",
        margin={"l": 10, "r": 10, "t": 50, "b": 10},
    )
    return apply_dark_layout(figure)


def render_scatter(drug_name, gene, gene_symbol, key_prefix, stored_row=None):
    """Plot expression against AUC, with a cancer-type filter and a live correlation."""
    points = points_for(drug_name, gene)
    if points.empty:
        st.warning("No cell line has both expression for this gene and an AUC for this drug.")
        return

    cancer_types = sorted(points["OncotreeLineage"].unique())
    selected = st.multiselect(
        "Filter by cancer type",
        cancer_types,
        key=f"{key_prefix}_cancer",
        placeholder="All cancer types",
        help="Leave this empty to include every cancer type. Choose one or more types to recalculate the correlation on that subset.",
    )
    shown = points.loc[points["OncotreeLineage"].isin(selected)] if selected else points
    if shown.empty:
        st.info("None of the selected cancer types have a point for this pair. Choose a different type.")
        return

    result = correlation_of(shown)
    filtered = bool(selected)
    first, second, third = st.columns(3)
    first.metric("Cell lines shown", f"{len(shown):,}")
    if result is None:
        second.metric("Strength of link", "—")
        third.metric("Confidence", "—")
        st.caption("The correlation needs at least three cell lines, and the values cannot all be identical.")
    else:
        coefficient, p_value = result
        second.metric("Strength of link", format_r(coefficient, digits=3))
        if stored_row is not None and not filtered:
            third.metric("Confidence", confidence_label(float(stored_row["fdr"])))
            st.caption(
                f"FDR = {float(stored_row['fdr']):.2e}, from the original test of every gene, across all cancer types. "
                "A negative link means higher expression goes with a lower AUC (more sensitive)."
            )
        else:
            third.metric("p-value", format_p(p_value))
            if filtered:
                st.caption(
                    "This correlation uses only the cancer types selected above. "
                    "A negative value means higher expression goes with a lower AUC."
                )
            else:
                st.caption(
                    "This gene was outside the top 50 for this drug, so there is no saved FDR. "
                    "The p-value here is for this one comparison."
                )
        direction_pill(coefficient)

    st.plotly_chart(
        scatter_figure(shown, drug_name, gene_symbol),
        width="stretch",
        theme=None,
        alt=f"Scatter plot of {gene_symbol} expression versus {drug_name} AUC, colored by cancer type.",
    )
    st.caption(
        "Each point is one cell line. Hover a point for its name. "
        "Point color is the cancer type. Sky blue and pink are reserved for sensitivity and resistance, so they are not used as cancer-type colors."
    )
