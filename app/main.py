"""Biomarker Explorer.

From the project folder, with the virtual environment on:

    streamlit run app/main.py

Then open http://localhost:8501
"""

import streamlit as st

from views.about import render as about
from views.case_study import render as case_study
from views.explore import render as explore
from views.home import render as home
from views.search_drug import render as search_drug
from views.search_gene import render as search_gene

st.set_page_config(
    page_title="Biomarker Explorer",
    layout="wide",
    initial_sidebar_state="expanded",
)

navigation = st.navigation(
    [
        st.Page(home, title="Home", icon=":material/home:", url_path="home", default=True),
        st.Page(search_drug, title="Search by drug", icon=":material/medication:", url_path="search-by-drug"),
        st.Page(search_gene, title="Search by gene", icon=":material/genetics:", url_path="search-by-gene"),
        st.Page(explore, title="Explore", icon=":material/scatter_plot:", url_path="explore"),
        st.Page(case_study, title="Case study: Venetoclax", icon=":material/menu_book:", url_path="case-study"),
        st.Page(about, title="About & limitations", icon=":material/info:", url_path="about"),
    ]
)
navigation.run()
