# DESeq2 Heatmap Generator

Client-side web app that parses DESeq2 output and renders a clustered, row-scaled heatmap with PDF export.

## How To Use

Download the app to your computer. Unzip if necessary.

Open `index.html` in a browser, upload files, set cutoffs, click **Generate Heatmap**, then export PDF if needed.

## Input Requirements

Upload a DESeq2 result file (`.csv`, `.tsv`, or `.txt`) containing at least:

- `gene_name`
- `padj`
- `log2FoldChange`
- `baseMean`

** Expression values are read from columns strictly between `gene_name` and `baseMean`.**

Optional:

- Gene list file (`.csv`, `.tsv`, or plain text)
  - one gene per line, or
  - table where the first column is gene name


## Heatmap Behavior

- Row-scaled (z-score per gene)
- Both rows and columns clustered

## PDF Export

- Click **Export as PDF** after generating a heatmap
- Output filename: `deseq2_heatmap.pdf`
