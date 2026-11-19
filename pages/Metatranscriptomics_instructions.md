# Metatranscriptomics Page — User Guide

## Strain & Thresholds

- Select a **strain** from the dropdown at the top of the page.
- Set the **FDR cutoff** and **Log FC cutoff** to define what counts as a significant hit. Genes meeting both thresholds are shown as teal points on the volcano plots and appear in the tables by default.

---

## Contrasts

- Choose two contrasts using the **Contrast 1** and **Contrast 2** selectors.
- The **left volcano** shows Contrast 1; the **right volcano** shows Contrast 2.
- The **middle scatter plot** plots log2FC of Contrast 1 (x-axis) against log2FC of Contrast 2 (y-axis), showing only genes that are significant hits in at least one contrast. Point color indicates whether a gene is a hit in Contrast 1 only, Contrast 2 only, or both.

---

## Highlighting Genes by Annotation

Use any combination of the three highlight dropdowns to mark genes of interest across all plots:

| Dropdown | Annotation used |
|---|---|
| Highlight by KEGG Pathway | KEGG metabolic pathway membership |
| Highlight by Substrate | Predicted substrate |
| Highlight by antiSMASH Product | Biosynthetic gene cluster product type |

- Highlighted genes appear in **orange** on all three plots, with larger markers.
- Selecting multiple entries highlights genes matching **any** of them.
- Highlighted genes are listed first in the tables below the plots.

---

## Selecting Genes from Plots

- **Click** a single point to select it.
- **Click and drag** to draw a box around multiple points.
- To **deselect**, double-click anywhere on the plot.
- Selections from all three plots are **combined** — you can select genes from the volcano and scatter simultaneously.
- Selected genes are shown in the tables below.

---

## Tables

Two tables are shown side-by-side, one per contrast.

- **Default view**: all significant hits (passing FDR and log FC thresholds), sorted by p-value.
- **When genes are selected or highlighted**: only those genes are shown.

### Columns

| Column | Description |
|---|---|
| Protein ID | NCBI protein accession — click to open the NCBI protein page |
| Locus Tag | Genome locus tag identifier |
| log2FoldChange | Log2 fold change for this contrast |
| padj | Benjamini-Hochberg adjusted p-value |
| Preferred Name | Short gene name (e.g. *rpsB*) |
| KEGG ko | KEGG orthology entry — click to open the KEGG database page |
| KEGG Pathway | KEGG metabolic pathway(s) |
| Substrate | Predicted substrate |
| antiSMASH Product | Biosynthetic gene cluster product type |
| Description | Functional description of the gene |
