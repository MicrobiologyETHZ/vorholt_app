# Vorholt App Development Notes

## Last Updated: 2026-02-27

## Overview
This is a Streamlit app for visualizing RNA-seq and metatranscriptomic data. The app includes multiple experiment pages (AK01, AK06, AK07) and a new Metatranscriptomics page.

---

## What Was Done (2026-02-27)

### 1. Annotation Integration

Added a separate annotations CSV per strain (e.g. `L137-annotations.csv`) stored alongside results in GCS. Annotations are merged into results on `locus_tag` via `enrich_results()`.

**Annotation file format** (columns stored as Python list-string literals, e.g. `"['ko:K02963']"`):
- `locus_tag` — join key
- `protein_id` — NCBI protein accession
- `Preferred_name` — human-readable gene name
- `KEGG_ko_list` — KEGG orthology codes
- `KEGG_Pathway_list` — KEGG pathway names
- `substrate` — substrate annotations
- `antismash_product` — antiSMASH biosynthetic cluster product

**`parse_list_str()`** helper parses list-string literals with `ast.literal_eval`, returning `[]` on NaN or error.

**`enrich_results()`** (cached): processes annotation columns (`ann_*` prefix to avoid conflicts), adds clickable KEGG URL, merges on `locus_tag`. Drops `protein_id` from results before merge so annotations file is authoritative.

### 2. Annotation-Based Highlight Controls

Added two multiselects (`Highlight by KEGG Pathway`, `Highlight by Substrate`) that collect matching `locus_tag`s from annotations and pass the corresponding `Gene` IDs as `genes_to_highlight` to all three plots.

Highlighted genes appear in **orange (#e07b39)** with larger markers. Color/size encoding:
- No highlights → fixed `marker.size=8`, no bubble mode
- Highlights present → `size='GOI'` (10 for background, 60 for highlighted), `size_max=25`

Scatter plot additionally checks `highlight_in_data = bool(highlight_set & set(merged[gene_name]))` before using bubble mode to avoid all-dots-oversized bug when highlighted genes are not hits.

### 3. Table Improvements

- **Clickable protein_id**: replaces `Gene` column; links to `https://www.ncbi.nlm.nih.gov/protein/{id}` using `LinkColumn(display_text=r"/protein/(.*)")`
- **locus_tag column**: shown for investigation of missing annotations
- **Annotation columns shown**: Preferred Name, KEGG ko (clickable), KEGG Pathway, Substrate, antiSMASH Product
- Hit banner shows clickable protein_id links (up to 10) instead of internal gene IDs

### 4. GCS Bucket Loading

Implemented `load_from_bucket_mt(bucket_name, folder, strain)` (TTL=600s):
- Lists blobs under `{folder}/` in bucket
- Matches `{strain}-results` and `{strain}-annotations` in blob name
- Returns two DataFrames: results, annotations

**Bucket structure** (`vorholt/BARA23-1/`):
```
BARA23-1/
├── L130-results.csv
├── L130-annotations.csv
├── L137-results.csv
├── L137-annotations.csv
├── L257-results.csv
├── L257-annotations.csv
├── L415-results.csv
├── L415-annotations.csv
├── L68-results.csv
└── L68-annotations.csv
```

**Results file columns** (new strains): `Gene`, `baseMean`, `log2FoldChange`, `padj`, `contrast`, `locus_tag`
`protein_id` is only in the annotations file.

### 5. Multi-Strain Support

`STRAIN_DATA` now covers 5 strains:
```python
STRAIN_DATA = {
    "L130": {"bucket": "vorholt", "folder": "BARA23-1"},
    "L137": {"bucket": "vorholt", "folder": "BARA23-1"},
    "L257": {"bucket": "vorholt", "folder": "BARA23-1"},
    "L415": {"bucket": "vorholt", "folder": "BARA23-1"},
    "L68":  {"bucket": "vorholt", "folder": "BARA23-1"},
}
```

### 6. Strain Selector: Tabs → Selectbox

Originally used `st.tabs()` for strain selection. This caused **Plotly frowny face (crash)** for L257 because `st.tabs()` renders all tab content simultaneously — 5 strains × 3 charts = 15 Plotly charts loaded at once, causing browser memory pressure.

**Fix**: replaced `st.tabs()` with `st.selectbox()`. Only the selected strain's 3 charts render at a time.

Also fixed: y-axis label changed from `-10*log10(padj)` (non-standard) to `-log10(padj)` (standard convention).

### 7. Known Data Issue: L257 Missing Annotations

458 `locus_tag`s in L257 results are absent from L257 annotations file → `protein_id` is NA for those genes. This is a **data problem**, not a code issue. The `locus_tag` column in the table is visible to help investigate.

---

## What Was Done (2026-02-24)

### 1. Dependency Updates
- Updated `requirements.txt` to modern versions:
  - streamlit >= 1.54.0 (was 1.10.0)
  - numpy >= 1.24.0
  - pandas >= 2.0.0
  - plotly >= 5.18.0
  - scipy >= 1.11.0
  - scikit-learn >= 1.3.0
  - google-cloud-storage >= 2.14.0

### 2. Fixed Deprecated Streamlit APIs
All deprecated caching decorators were replaced:
- `@st.cache` → `@st.cache_data`
- `@st.experimental_memo` → `@st.cache_data`

Files updated:
- `processing/load_data.py`
- `processing/eda.py`
- `processing/connect_to_gsheet.py`
- `graphs/show_dge.py`

### 3. Re-enabled Experiment Pages
- Uncommented and restored `pages/AK01.py`, `pages/AK06.py`, `pages/AK07.py`
- These pages load data from GCS bucket 'vorholt'

### 4. Cleanup
- Removed unused `app.py` (was a Dash app, not Streamlit)
- Cleaned up `Home.py` - removed unnecessary Google Cloud imports that were causing slow loading

### 5. Created Metatranscriptomics Page (`pages/Metatranscriptomics.py`)

#### Features Implemented:
- **Wide layout** with custom CSS for maximum width
- **Tab-based strain selection** (currently only L137 configured)
- **Controls row**: FDR cutoff, Log FC cutoff, gene highlighter, contrast selectors
- **3 graphs in a row**:
  1. Contrast 1 Volcano plot
  2. log2FC Comparison scatter plot (Contrast 1 vs Contrast 2)
  3. Contrast 2 Volcano plot
- **Interactive selection**: Click/drag to select points on any graph, selection filters both tables
- **2 tables side-by-side**: Hits for Contrast 1 and Contrast 2
- **Table columns**: Gene, log2FoldChange, padj, KEGG_ko, KEGG_Pathway, Preferred_name, Description, EC, PFAMs
- **Clickable KEGG_ko links**: Values starting with "ko:" link to KEGG database
- **Download button** for comparison data CSV

#### Data Format Expected:
CSV file with columns:
- `Gene` - gene identifier
- `log2FoldChange` - fold change value
- `padj` - adjusted p-value
- `contrast` - contrast name
- Optional annotation columns: `KEGG_ko`, `KEGG_Pathway`, `Preferred_name`, `Description`, `EC`, `PFAMs`

#### Current Data Configuration:
```python
STRAIN_DATA = {
    "L137": "test_data/2026-02-24_L137-results.csv",
}
```

---

## What Still Needs To Be Done

### High Priority
1. ~~**Add more strains**~~ — Done: L130, L137, L257, L415, L68 all loading from GCS
2. ~~**Test GCS bucket loading**~~ — Done: `load_from_bucket_mt()` implemented and tested
3. **Verify AK01/AK06/AK07 pages work** - These pages use GCS bucket and haven't been tested after the updates
4. **Fix L257 annotations** — 458 `locus_tag`s in L257 results have no match in the annotations file; protein_id is NA for those genes. The annotations CSV needs to be regenerated/extended to cover all locus_tags in the results.

### Medium Priority
5. **Data pivot script** - User has wide-format data that needs pivoting to long format. Code provided:
```python
# Extract contrast names from column headers
contrasts = [col.replace('log2FoldChange_', '')
             for col in df.columns if col.startswith('log2FoldChange_')]

# Build long format
long_dfs = []
for contrast in contrasts:
    temp = df[['Geneid', 'baseMean']].copy()
    temp['log2FoldChange'] = df[f'log2FoldChange_{contrast}']
    temp['padj'] = df[f'padj_{contrast}']
    temp['contrast'] = contrast
    long_dfs.append(temp)

df_long = pd.concat(long_dfs, ignore_index=True)
df_long = df_long.rename(columns={'Geneid': 'Gene'})
```

6. **Update projects.csv** - Add metatranscriptomics projects to the home page table

### Low Priority / Nice to Have
6. **Add PCA visualization** to Metatranscriptomics page (like other experiment pages)
7. **Add gene expression plots** to Metatranscriptomics page
8. **STRING PPI network integration** for metatranscriptomics (exists in other pages)

---

## File Structure

```
vorholt_app/
├── Home.py                     # Main entry point
├── check_password.py           # Authentication
├── requirements.txt            # Dependencies
├── projects.csv                # Project listing
├── DEVELOPMENT_NOTES.md        # This file
├── .streamlit/
│   └── secrets.toml            # GCP credentials (gitignored)
├── test_data/
│   └── 2026-02-24_L137-results.csv  # Test metatranscriptomic data
├── pages/
│   ├── AK01.py                 # RNA-seq experiment page
│   ├── AK06.py                 # RNA-seq experiment page
│   ├── AK07.py                 # RNA-seq experiment page
│   └── Metatranscriptomics.py  # New metatranscriptomics page
├── processing/
│   ├── load_data.py            # Data loading from local/GCS
│   ├── eda.py                  # PCA computation
│   └── connect_to_gsheet.py    # Google Sheets connector
└── graphs/
    ├── show_pca.py             # PCA visualization
    ├── show_gene_expression.py # Expression plots
    └── show_dge.py             # Volcano plots, STRING integration
```

---

## Running the App

```bash
conda activate vorholt_app  # the streamlit env had numpy/pandas incompatibility
cd /Users/ansintsova/git_repos/vorholt_app
streamlit run Home.py
```

---

## Key Technical Notes

### Graph Selection Implementation
- Uses Plotly's `on_select="rerun"` with `custom_data` to capture selected genes
- `custom_data=[gene_name]` stores gene names for reliable selection (point_index is per-trace, not per-dataframe)
- Selections from all 3 graphs are combined to filter both tables

### KEGG Links
- Values like "ko:K9890,ko:K0282" are split, only first value used for link
- URL format: `https://www.genome.jp/dbget-bin/www_bget?ko:K9890`
- Uses `st.column_config.LinkColumn(display_text=r"\?(.*)")` to show original value

### Caching
- `@st.cache_data` used for data loading functions
- Data is cached with TTL where appropriate (600s for GCS data)
