# DESeq2 Heatmap Generator (Shiny + pheatmap)

Web app built with R Shiny and `pheatmap` to generate row-scaled clustered heatmaps from DESeq2 output.

## Inputs

- **DESeq2 result file** (`csv/tsv/txt`, required) with:
  - `gene_name` (used as heatmap row label)
  - `padj`
  - `log2FoldChange`
- **Expression columns priority**:
  - first tries expression columns inside DESeq2 (columns between `gene_name` and `baseMean`)
  - if none exist, uses optional expression matrix file fallback
- **Optional selected genes file**:
  - one gene per line, or a table where the first column is gene ID

## Filters

- `padj.cutoff`: must be within `0` to `1`
- `log2FC.cutoff`: must be within `0` to max `|log2FoldChange|` from the DESeq2 file

## Plot behavior

- genes are filtered by:
  - `padj <= padj.cutoff`
  - `|log2FoldChange| >= log2FC.cutoff`
- if selected genes file is provided, genes must also be in that list
- each gene row is z-score scaled across samples (`scale = "row"`)
- row and column clustering are enabled (`cluster_rows = TRUE`, `cluster_cols = TRUE`)
- color palette uses a Qlik Sense-style diverging scheme

## Usage

Install packages (once):

`install.packages(c("shiny", "pheatmap"))`

Run the app:

`shiny::runApp("app.R")`

Then upload files, set thresholds, click **Generate Heatmap**, and use **Export PDF** to save `deseq2_heatmap.pdf`.
