import streamlit as st
from pathlib import Path
import pandas as pd
from graphs.show_pca import show_pca
from graphs.show_gene_expression import show_expression
from graphs.show_dge import show_volcano, link_to_string, download_filtered_hits
from processing.load_data import load_data, load_from_bucket
from check_password import check_password

datadir = "/Users/ansintsova/polybox/Vorholt/akeppler/22AK07"
page_name = "AK07"
sample_id = 'Novogene_sample'
annotation_cols = ()


def app():
    st.title(page_name)
    results, tpms, vsd, sd = load_from_bucket('vorholt', experiment_name='AK07', gene_name='')
    #results, tpms, vsd, sd = load_data(datadir, gene_name='')

    st.header("PCA")
    with st.expander('Show PCA'):
        show_pca(tpms, sd)
    st.header("Expression")
    with st.expander('Show Gene Expression'):
        genes = show_expression(tpms, sd, annotation_cols=annotation_cols, sample_col=sample_id)

if check_password():
    app()


