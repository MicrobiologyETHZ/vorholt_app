import streamlit as st

st.set_page_config(layout="wide", page_title="Metatranscriptomics")

import ast
from io import StringIO
import pandas as pd
import plotly.express as px
import numpy as np
from google.oauth2 import service_account
from google.cloud import storage
from check_password import check_password

# Custom CSS to reduce padding and maximize width
st.markdown("""
    <style>
        .block-container {
            padding-top: 1rem;
            padding-bottom: 1rem;
            padding-left: 2rem;
            padding-right: 2rem;
            max-width: 100%;
        }
    </style>
""", unsafe_allow_html=True)

page_name = "Metatranscriptomics"

# Strain data — use "bucket"/"folder" for GCS, or "results"/"annotations" for local files
STRAIN_DATA = {
    "L130": {"bucket": "vorholt", "folder": "BARA23-1"},
    "L137": {"bucket": "vorholt", "folder": "BARA23-1"},
    "L257": {"bucket": "vorholt", "folder": "BARA23-1"},
    "L415": {"bucket": "vorholt", "folder": "BARA23-1"},
    "L68":  {"bucket": "vorholt", "folder": "BARA23-1"},
}


@st.cache_data
def load_results(file_path):
    return pd.read_csv(file_path, index_col=0)


@st.cache_data
def load_annotations(file_path):
    return pd.read_csv(file_path, index_col=0)


@st.cache_data(ttl=600)
def load_from_bucket_mt(bucket_name, folder, strain):
    """Load results and annotations CSVs for a strain from a GCS bucket folder."""
    credentials = service_account.Credentials.from_service_account_info(
        st.secrets["gcp_service_account"]
    )
    client = storage.Client(credentials=credentials)

    blobs = [b for b in client.list_blobs(bucket_name, prefix=folder + '/') if b.size > 0]

    results_blob = next((b for b in blobs if f'{strain}-results' in b.name), None)
    ann_blob     = next((b for b in blobs if f'{strain}-annotations' in b.name), None)

    def read_blob(blob):
        content = blob.download_as_string().decode('utf-8')
        return pd.read_csv(StringIO(content), index_col=0)

    if results_blob is None:
        raise FileNotFoundError(f"No results CSV for {strain} in {bucket_name}/{folder}/")
    if ann_blob is None:
        raise FileNotFoundError(f"No annotations CSV for {strain} in {bucket_name}/{folder}/")

    return read_blob(results_blob), read_blob(ann_blob)


def parse_list_str(s):
    """Parse a Python list string like \"['foo', 'bar']\" into an actual list."""
    if pd.isna(s):
        return []
    try:
        result = ast.literal_eval(str(s))
        return result if isinstance(result, list) else [result]
    except Exception:
        return []


@st.cache_data
def enrich_results(results, annotations):
    """Merge results with processed annotation columns."""
    ann = annotations.copy()

    def first_val(s):
        vals = parse_list_str(s)
        return vals[0] if vals else None

    def join_vals(s):
        vals = parse_list_str(s)
        return ', '.join(str(v) for v in vals) if vals else None

    ann['ann_Preferred_name'] = ann['Preferred_name'].apply(first_val)
    ko_first = ann['KEGG_ko_list'].apply(first_val)
    ann['ann_KEGG_ko_url'] = ko_first.apply(
        lambda ko: f"https://www.genome.jp/dbget-bin/www_bget?{ko}"
        if ko and str(ko).startswith('ko:') else None
    )
    ann['ann_KEGG_Pathway'] = ann['KEGG_Pathway_list'].apply(join_vals)
    ann['ann_substrate'] = ann['substrate'].apply(join_vals)
    ann['ann_antismash_product'] = ann['antismash_product'].apply(join_vals)

    ann_cols = ['locus_tag', 'protein_id', 'Description', 'ann_Preferred_name', 'ann_KEGG_ko_url',
                'ann_KEGG_Pathway', 'ann_substrate', 'ann_antismash_product']
    # Drop protein_id from results if already present (annotations is authoritative)
    results_clean = results.drop(columns=['protein_id'], errors='ignore')
    return results_clean.merge(ann[ann_cols], on='locus_tag', how='left')


@st.cache_data
def get_unique_values(annotations, col):
    """Get sorted unique non-null values from a list-string column."""
    all_vals = set()
    for val in annotations[col].dropna():
        all_vals.update(parse_list_str(val))
    return sorted(str(v) for v in all_vals if v)



def show_volcano_plot(volcano_df, contrast_name, key_prefix, gene_name='Gene', lfc_col='log2FoldChange',
                      pval_col='padj', genes_to_highlight=None, fdr=0.05, lfc_th=1.0):
    """Display volcano plot with selection support. Returns selected genes."""
    plot_df = volcano_df.copy().reset_index(drop=True)
    plot_df['log10FDR'] = -np.log10(plot_df[pval_col].clip(lower=1e-300))

    highlight_set = set(genes_to_highlight) if genes_to_highlight else set()

    if highlight_set:
        def _category(row):
            if row[gene_name] in highlight_set:
                return 'Highlighted'
            return 'Hit' if row['Hit'] else 'Background'
        plot_df['Category'] = plot_df.apply(_category, axis=1)
        plot_df['GOI'] = plot_df[gene_name].isin(highlight_set).astype(int) * 50 + 10
        size_kwargs = dict(size='GOI', size_max=25)
    else:
        plot_df['Category'] = plot_df['Hit'].map({True: 'Hit', False: 'Background'})
        size_kwargs = {}

    color_map = {'Background': '#d3d3d3', 'Hit': '#378b84', 'Highlighted': '#e07b39'}
    cat_order = {'Category': ['Background', 'Hit', 'Highlighted']}

    hover_cols = {
        lfc_col: ':.2f',
        pval_col: ':.2e',
        'log10FDR': False,
        'Hit': False,
        'Category': False,
    }
    for col in ['ann_Preferred_name', 'Description']:
        if col in plot_df.columns:
            hover_cols[col] = True

    fig = px.scatter(
        plot_df,
        x=lfc_col,
        y='log10FDR',
        color='Category',
        height=600,
        title=f"Volcano: {contrast_name}",
        color_discrete_map=color_map,
        category_orders=cat_order,
        hover_name=plot_df[gene_name],
        hover_data=hover_cols,
        custom_data=[gene_name],
        **size_kwargs
    )

    fig.add_vline(x=lfc_th, line_width=2, line_dash="dash", line_color="grey")
    fig.add_vline(x=-lfc_th, line_width=2, line_dash="dash", line_color="grey")
    fig.add_hline(y=-np.log10(fdr), line_width=2, line_dash="dash", line_color="grey")

    fig.update_layout(
        {'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
        autosize=True,
        font=dict(size=14),
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(padj)",
        dragmode='select'
    )
    base_marker = dict(line=dict(width=1, color='DarkSlateGrey'), opacity=0.8)
    if not highlight_set:
        base_marker['size'] = 8
    fig.update_traces(marker=base_marker, selector=dict(mode='markers'))

    event = st.plotly_chart(fig, use_container_width=True, on_select="rerun", key=f"{key_prefix}_volcano")

    selected_genes = []
    if event and event.selection and event.selection.points:
        for point in event.selection.points:
            customdata = point.get('customdata')
            if customdata and len(customdata) > 0:
                selected_genes.append(customdata[0])

    return selected_genes


def show_comparison_plot(results, contrast1, contrast2, key_prefix, gene_name='Gene',
                         lfc_col='log2FoldChange', pval_col='padj', fdr=0.05, lfc_th=1.0,
                         genes_to_highlight=None):
    """Display log2FC comparison between two contrasts. Returns (merged_df, selected_genes)."""
    df1 = results[results['contrast'] == contrast1][[gene_name, lfc_col, pval_col]].copy()
    df1.columns = [gene_name, 'lfc1', 'padj1']

    df2 = results[results['contrast'] == contrast2][[gene_name, lfc_col, pval_col]].copy()
    df2.columns = [gene_name, 'lfc2', 'padj2']

    merged = df1.merge(df2, on=gene_name, how='outer')
    merged['hit1'] = (abs(merged['lfc1']) > lfc_th) & (merged['padj1'] < fdr)
    merged['hit2'] = (abs(merged['lfc2']) > lfc_th) & (merged['padj2'] < fdr)
    merged = merged[merged['hit1'] | merged['hit2']].copy().reset_index(drop=True)

    if merged.empty:
        st.warning("No significant hits found in either contrast with current thresholds.")
        return pd.DataFrame(), []

    def categorize(row):
        if row['hit1'] and row['hit2']:
            return 'Both'
        elif row['hit1']:
            return f'{contrast1} only'
        else:
            return f'{contrast2} only'

    merged['Category'] = merged.apply(categorize, axis=1)

    # Size encodes highlight status (color is already used for contrast category)
    highlight_set = set(genes_to_highlight) if genes_to_highlight else set()
    merged['GOI'] = merged[gene_name].isin(highlight_set).astype(int) * 50 + 10
    # Only use size encoding when highlighted genes are actually present in the hits data;
    # if none are hits, all GOI values are equal and Plotly would scale all dots to size_max.
    highlight_in_data = bool(highlight_set & set(merged[gene_name]))
    size_max = 25 if highlight_in_data else None

    # Pull annotation columns from results for hover
    for col in ['ann_Preferred_name', 'Description']:
        if col in results.columns:
            name_map = results.drop_duplicates(gene_name).set_index(gene_name)[col].to_dict()
            merged[col] = merged[gene_name].map(name_map)

    color_map = {
        f'{contrast1} only': '#e41a1c',
        f'{contrast2} only': '#377eb8',
        'Both': '#984ea3'
    }

    hover_data = {
        'lfc1': ':.2f', 'lfc2': ':.2f', 'padj1': ':.2e', 'padj2': ':.2e',
        'GOI': False, 'hit1': False, 'hit2': False,
    }
    for col in ['ann_Preferred_name', 'Description']:
        if col in merged.columns:
            hover_data[col] = True

    scatter_size_kwargs = dict(size='GOI', size_max=size_max) if highlight_in_data else {}

    fig = px.scatter(
        merged,
        x='lfc1',
        y='lfc2',
        color='Category',
        height=600,
        title=f"log2FC Comparison: {contrast1} vs {contrast2}",
        color_discrete_map=color_map,
        category_orders={'Category': [f'{contrast1} only', f'{contrast2} only', 'Both']},
        hover_name=gene_name,
        hover_data=hover_data,
        custom_data=[gene_name],
        **scatter_size_kwargs
    )

    max_val = max(abs(merged['lfc1'].max()), abs(merged['lfc2'].max()),
                  abs(merged['lfc1'].min()), abs(merged['lfc2'].min()))
    fig.add_shape(
        type='line',
        x0=-max_val, y0=-max_val,
        x1=max_val, y1=max_val,
        line=dict(color='grey', dash='dash', width=1)
    )

    fig.update_layout(
        {'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
        autosize=True,
        font=dict(size=14),
        xaxis_title=f"log2FC ({contrast1})",
        yaxis_title=f"log2FC ({contrast2})",
        dragmode='select'
    )
    base_marker = dict(line=dict(width=1, color='DarkSlateGrey'), opacity=0.8)
    if not highlight_in_data:
        base_marker['size'] = 8
    fig.update_traces(marker=base_marker, selector=dict(mode='markers'))

    event = st.plotly_chart(fig, use_container_width=True, on_select="rerun", key=f"{key_prefix}_compare")

    selected_genes = []
    if event and event.selection and event.selection.points:
        for point in event.selection.points:
            customdata = point.get('customdata')
            if customdata and len(customdata) > 0:
                selected_genes.append(customdata[0])

    counts = merged['Category'].value_counts()
    c1, c2, c3 = st.columns(3)
    c1.metric(f"{contrast1} only", counts.get(f'{contrast1} only', 0))
    c2.metric(f"{contrast2} only", counts.get(f'{contrast2} only', 0))
    c3.metric("Both", counts.get('Both', 0))

    return merged, selected_genes


def show_hits_table(volcano_df, display_genes, contrast_name):
    """Render hits table for one contrast."""
    if display_genes:
        table_df = volcano_df[volcano_df['Gene'].isin(display_genes)].copy()
        label = f"{len(table_df)} selected"
    else:
        table_df = volcano_df[volcano_df['Hit'] == True].copy()
        label = f"{len(table_df)} significant hits"

    if table_df.empty:
        st.warning("No genes to display.")
        return

    table_df = table_df.sort_values('padj', ascending=True)

    # Build clickable protein_id column
    if 'protein_id' in table_df.columns:
        table_df['protein_id_url'] = table_df['protein_id'].apply(
            lambda x: f"https://www.ncbi.nlm.nih.gov/protein/{x}" if pd.notna(x) else None
        )
        id_col = ['protein_id_url']
    else:
        id_col = ['Gene']

    ann_cols = ['ann_Preferred_name', 'ann_KEGG_ko_url', 'ann_KEGG_Pathway',
                'ann_substrate', 'ann_antismash_product', 'Description']
    locus_col = ['locus_tag'] if 'locus_tag' in table_df.columns else []
    display_cols = id_col + locus_col + ['log2FoldChange', 'padj'] + [c for c in ann_cols if c in table_df.columns]

    column_config = {
        "protein_id_url": st.column_config.LinkColumn("Protein ID", display_text=r"/protein/(.*)"),
        "locus_tag": st.column_config.TextColumn("Locus Tag"),
        "log2FoldChange": st.column_config.NumberColumn(format="%.2f"),
        "padj": st.column_config.NumberColumn(format="%.2e"),
        "ann_Preferred_name": st.column_config.TextColumn("Preferred Name"),
        "ann_KEGG_ko_url": st.column_config.LinkColumn("KEGG ko", display_text=r"\?(.*)"),
        "ann_KEGG_Pathway": st.column_config.TextColumn("KEGG Pathway"),
        "ann_substrate": st.column_config.TextColumn("Substrate"),
        "ann_antismash_product": st.column_config.TextColumn("antiSMASH Product"),
        "Description": st.column_config.TextColumn("Description"),
    }

    st.dataframe(
        table_df[display_cols],
        height=400,
        use_container_width=True,
        column_config=column_config
    )
    st.caption(label)


def app():
    st.title(page_name)
    st.markdown("Visualize differential expression results from metatranscriptomic data.")

    strain = st.selectbox("Strain", list(STRAIN_DATA.keys()), key="strain_select")
    strain_info = STRAIN_DATA[strain]

    # Load results + annotations (bucket or local)
    try:
        if 'bucket' in strain_info:
            results, annotations = load_from_bucket_mt(
                strain_info['bucket'], strain_info['folder'], strain
            )
        else:
            results = load_results(strain_info['results'])
            ann_path = strain_info.get('annotations')
            annotations = load_annotations(ann_path) if ann_path else None
    except Exception as e:
        st.error(f"Error loading data: {e}")
        st.stop()

    # Merge annotations into results
    try:
        if annotations is not None and not annotations.empty:
            results = enrich_results(results, annotations)
    except Exception as e:
        st.warning(f"Could not merge annotations: {e}")
        annotations = None

    if results.empty:
        st.warning("No results found in the data file.")
        st.stop()

    st.info(f"Loaded {len(results):,} entries across {results['contrast'].nunique()} contrast(s)")

    # Unique values for annotation filters
    all_pathways = get_unique_values(annotations, 'KEGG_Pathway_list') if annotations is not None else []
    all_substrates = get_unique_values(annotations, 'substrate') if annotations is not None else []
    all_antismash = get_unique_values(annotations, 'antismash_product') if annotations is not None else []

    # === CONTROLS ===
    c1, c2 = st.columns(2)
    fdr = c1.number_input('FDR cutoff', value=0.05, min_value=0.0, max_value=1.0,
                          step=0.01, key=f'{strain}_fdr')
    lfc_th = c2.number_input('Log FC cutoff', value=1.0, min_value=0.0,
                             step=0.5, key=f'{strain}_lfc')

    # Annotation-based highlight filters
    c4, c5, c6 = st.columns(3)
    selected_pathways = c4.multiselect("Highlight by KEGG Pathway", all_pathways,
                                       default=[], key=f'{strain}_pathways')
    selected_substrates = c5.multiselect("Highlight by Substrate", all_substrates,
                                         default=[], key=f'{strain}_substrates')
    selected_antismash = c6.multiselect("Highlight by antiSMASH Product", all_antismash,
                                        default=[], key=f'{strain}_antismash')

    # Contrast selectors
    contrasts = list(results['contrast'].unique())
    c1, c2 = st.columns(2)
    contrast1 = c1.selectbox('Contrast 1', contrasts, index=0, key=f'{strain}_c1')
    contrast2 = c2.selectbox('Contrast 2', contrasts, index=min(1, len(contrasts) - 1),
                             key=f'{strain}_c2')

    # Compute filter_genes from annotation selections
    filter_genes = set()
    if annotations is not None and (selected_pathways or selected_substrates or selected_antismash):
        locus_to_gene = (results.drop_duplicates('locus_tag')
                         .set_index('locus_tag')['Gene'].to_dict())

        if selected_pathways:
            selected_set = set(selected_pathways)
            mask = annotations['KEGG_Pathway_list'].apply(
                lambda s: bool(set(parse_list_str(s)) & selected_set)
            )
            for lt in annotations.loc[mask, 'locus_tag']:
                if lt in locus_to_gene:
                    filter_genes.add(locus_to_gene[lt])

        if selected_substrates:
            selected_set = set(selected_substrates)
            mask = annotations['substrate'].apply(
                lambda s: bool(set(parse_list_str(s)) & selected_set)
            )
            for lt in annotations.loc[mask, 'locus_tag']:
                if lt in locus_to_gene:
                    filter_genes.add(locus_to_gene[lt])

        if selected_antismash:
            selected_set = set(selected_antismash)
            mask = annotations['antismash_product'].apply(
                lambda s: bool(set(parse_list_str(s)) & selected_set)
            )
            for lt in annotations.loc[mask, 'locus_tag']:
                if lt in locus_to_gene:
                    filter_genes.add(locus_to_gene[lt])

    filter_genes = list(filter_genes)
    genes_to_highlight = filter_genes

    st.markdown("---")

    # Compute hits for each contrast
    volcano_df = results[results['contrast'] == contrast1].copy()
    volcano_df['Hit'] = (abs(volcano_df['log2FoldChange']) > lfc_th) & (volcano_df['padj'] < fdr)

    volcano_df2 = results[results['contrast'] == contrast2].copy()
    volcano_df2['Hit'] = (abs(volcano_df2['log2FoldChange']) > lfc_th) & (volcano_df2['padj'] < fdr)

    # === GRAPHS ===
    selected_genes1, selected_genes2, selected_genes3 = [], [], []

    if contrast1 != contrast2:
        graph_col1, graph_col2, graph_col3 = st.columns(3)

        with graph_col1:
            selected_genes1 = show_volcano_plot(
                volcano_df, contrast1, key_prefix=f"{strain}_v1",
                genes_to_highlight=genes_to_highlight, fdr=fdr, lfc_th=lfc_th
            )
        with graph_col2:
            _, selected_genes2 = show_comparison_plot(
                results, contrast1, contrast2, key_prefix=f"{strain}_comp",
                fdr=fdr, lfc_th=lfc_th, genes_to_highlight=genes_to_highlight
            )
        with graph_col3:
            selected_genes3 = show_volcano_plot(
                volcano_df2, contrast2, key_prefix=f"{strain}_v2",
                genes_to_highlight=genes_to_highlight, fdr=fdr, lfc_th=lfc_th
            )
    else:
        st.warning("Select two different contrasts to see plots.")
        selected_genes1 = show_volcano_plot(
            volcano_df, contrast1, key_prefix=f"{strain}_v1",
            genes_to_highlight=genes_to_highlight, fdr=fdr, lfc_th=lfc_th
        )

    # Combine graph-click selections with annotation filter genes for table display
    graph_selected = list(set(selected_genes1 + selected_genes2 + selected_genes3))
    display_genes = list(set(graph_selected + filter_genes))

    st.markdown("---")

    if display_genes:
        gene_to_protein = (results.drop_duplicates('Gene')
                           .set_index('Gene')['protein_id'].to_dict())
        links = []
        for g in sorted(display_genes)[:10]:
            pid = gene_to_protein.get(g)
            if pid:
                links.append(f"[{pid}](https://www.ncbi.nlm.nih.gov/protein/{pid})")
            else:
                links.append(g)
        suffix = ' ...' if len(display_genes) > 10 else ''
        st.info(f"Showing {len(display_genes)} gene(s): {', '.join(links)}{suffix}")

    # === TABLES ===
    table_col1, table_col2 = st.columns(2)

    with table_col1:
        st.markdown(f"**{contrast1}**")
        show_hits_table(volcano_df, display_genes, contrast1)

    with table_col2:
        st.markdown(f"**{contrast2}**")
        show_hits_table(volcano_df2, display_genes, contrast2)



if check_password():
    app()
