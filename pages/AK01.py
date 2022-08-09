import streamlit as st
from pathlib import Path
import pandas as pd
from processing.load_data import load_data, load_from_bucket
from graphs.show_pca import show_pca
from graphs.show_gene_expression import show_expression
from graphs.show_dge import show_volcano, link_to_string, download_filtered_hits
from check_password import check_password


datadir = "/Users/ansintsova/polybox/Vorholt/akeppler/21AK01"
page_name = "AK01"
sample_id = 'Sample'
annotation_cols = ()

#st.set_page_config(page_title=page_name, layout='wide')


def app():
    st.title(page_name)
    #results, tpms, vsd, sd = load_data(datadir, gene_name='')
    results, tpms, vsd, sd = load_from_bucket('vorholt', experiment_name='AK01', gene_name='')
    st.header("PCA")
    with st.expander('Show PCA'):
        show_pca(vsd, sd)
    st.header("Expression")
    with st.expander('Show Gene Expression'):
        tpms = tpms.drop(["SYMBOL", "KEGG_Pathway"], axis=1)
        genes = show_expression(tpms, sd, sample_col=sample_id)
    st.header("DGE Results")
    with st.expander('Show DGE Results'):
        c1, c2, c3 = st.columns([1, 3, 1])
        hits_df = show_volcano(results, c1, c2, genes_to_highlight=genes)
        c3.markdown("### STRING PPI Network")
        link_to_string(hits_df, c3)
        c3.markdown("### Download results")
        download_filtered_hits(hits_df, c3)

if check_password():
    app()


