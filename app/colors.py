"""One set of colours for every page and chart."""

SENSITIVE = "#5CC8FF"  # sky blue: genes linked to drug sensitivity
RESISTANT = "#FF7AB6"  # pink: genes linked to resistance
MUTED = "#7F8BA6"  # axis labels and gridlines
BACKGROUND = "#0B1220"
PANEL = "#111A2E"
TEXT = "#E4ECFA"


def apply_dark_layout(figure):
    """Dark Plotly layout: transparent chart, muted axes, light titles."""
    figure.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font={"color": TEXT, "family": "Source Sans 3, sans-serif"},
        title_font_color=TEXT,
        legend_font_color=TEXT,
        hoverlabel={"bgcolor": PANEL, "bordercolor": PANEL, "font": {"color": TEXT, "family": "Source Sans 3, sans-serif"}},
    )
    axis = {
        "gridcolor": MUTED,
        "linecolor": MUTED,
        "zerolinecolor": MUTED,
        "tickfont": {"color": MUTED},
        "title": {"font": {"color": MUTED}},
    }
    figure.update_xaxes(**axis)
    figure.update_yaxes(**axis)
    return figure
