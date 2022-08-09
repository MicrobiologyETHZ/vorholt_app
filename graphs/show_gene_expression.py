import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px


def show_expression(countData, sampleData, annotation_cols=(), sample_col="sampleID"):
    df = countData.copy()
    if annotation_cols:
        gene_name = st.radio('Choose gene annotation', annotation_cols)
        df = df.set_index(gene_name)
    else:
        df.index.name = 'Gene'
        gene_name = 'Gene'
    df = df.apply(lambda x: np.log2(x + 0.5) if np.issubdtype(x.dtype, np.number) else x)
    sampleDataAb = sampleData.reset_index()
    df = df.reset_index()
    c1, c2 = st.columns(2)
    compare_by = c1.selectbox('Compare by', sampleDataAb.columns)
    categories = c1.multiselect(f'Categories of {compare_by} to display',
                                ['All'] + list(sampleDataAb[compare_by].unique()), key='cats')
    filter_by = c2.selectbox("Filter by",  sampleDataAb.columns)
    filter_out = c2.selectbox(f'Which category of {filter_by} to keep?',
                               [None] + list(sampleDataAb[filter_by].unique()))

    genes = st.multiselect("Choose gene(s) of interest", df[gene_name].unique(), key='gois')

    if 'All' in categories:
        categories = list(sampleDataAb[compare_by].unique())
    if genes:
        if len(genes) * len(categories) > 40:
            st.write('Too many genes/categories to display, consider choosing fewer genes')
            st.stop()
        c3, c4 = st.columns(2)
        tpm_label = 'log2 (TPM)'
        gene_df = df[df[gene_name].isin(genes)]

        sample_df = sampleDataAb[sampleDataAb[compare_by].isin(categories)]
        if filter_out:
            sample_df = sample_df[sample_df[filter_by] == filter_out]
        gene_df2 = (gene_df.melt(id_vars=[gene_name], value_name=tpm_label, var_name=sample_col)
                    .merge(sample_df, how='inner', on=sample_col))
        groupby = st.radio('Group by', [gene_name, compare_by])
        color_by = [c for c in [gene_name, compare_by] if c != groupby][0]
        fig = px.box(gene_df2, x=groupby, y=tpm_label, color=color_by,
                     hover_data=[gene_name] + list(sampleData.columns), points='all')
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=16))
        fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
        st.plotly_chart(fig, use_container_width=True)
    return genes

def show_expression_heatmap():
    pass

