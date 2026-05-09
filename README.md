# DESeq2 Heatmap Generator (Standalone HTML/JS)

Client-side web app that parses DESeq2 output and renders a clustered, row-scaled heatmap with PDF export.

## Files

- `index.html`: app layout
- `styles.css`: styling
- `script.js`: parsing, filtering, clustering, plotting, and export logic
- `app.R`: older R Shiny implementation (kept unchanged)

## Input Requirements

Upload a DESeq2 result file (`.csv`, `.tsv`, or `.txt`) containing at least:

- `gene_name`
- `padj`
- `log2FoldChange`
- `baseMean`

Expression values are read from columns strictly between `gene_name` and `baseMean`.

Optional:

- Gene list file (`.csv`, `.tsv`, or plain text)
  - one gene per line, or
  - table where the first column is gene name

## Filtering Logic

1. If duplicate `gene_name` values exist, keep only the row with higher `log2FoldChange`.
2. If a gene list is uploaded, overlap it with DESeq2 `gene_name` first.
3. Apply cutoffs on the overlapped set:
   - `padj <= padj.cutoff`
   - `abs(log2FoldChange) >= log2FC.cutoff`

Input validation:

- `padj.cutoff` must be between `0` and `1`
- `log2FC.cutoff` must be between `0` and max `|log2FoldChange|`

## Heatmap Behavior

- Row-scaled (z-score per gene)
- Both rows and columns clustered
- Qlik-style diverging palette
- Gene labels:
  - shown when plotted genes are `< 50`
  - hidden when plotted genes are `>= 50`

## PDF Export

- Click **Export as PDF** after generating a heatmap
- Output filename: `deseq2_heatmap.pdf`

## Run

No build step required.

Open `index.html` in a browser, upload files, set cutoffs, click **Generate Heatmap**, then export PDF if needed.
