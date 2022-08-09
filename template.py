import streamlit as st
from pathlib import Path
import pandas as pd
from graphs.show_pca import show_pca
from graphs.show_gene_expression import show_expression
from graphs.show_dge import show_volcano, link_to_string, download_filtered_hits

datadir = "/Users/ansintsova/polybox/Vorholt/akeppler/22AK06"
page_name = "Template"
sample_id = 'sampleID'
annotation_cols = ()
gene_name = ''




def app():
    st.title(page_name)
    results, tpms, vsd, sd = load_data(datadir)
    st.header("PCA")
    with st.expander('Show PCA'):
        show_pca(vsd, sd)
    st.header("Expression")
    with st.expander('Show Gene Expression'):
        genes = show_expression(tpms, sd,annotation_cols=annotation_cols, sample_col=sample_id)
    st.header("DGE Results")
    with st.expander('Show DGE Results'):
        c1, c2, c3 = st.columns([1, 3, 1])
        hits_df = show_volcano(results, c1, c2, genes_to_highlight=genes)
        c3.markdown("### STRING PPI Network")
        link_to_string(hits_df, c3)
        c3.markdown("### Download results")
        download_filtered_hits(hits_df, c3)
app()


