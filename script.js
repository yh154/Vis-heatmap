const state = {
  deseqRows: [],
  selectedGenesSet: null,
  lastFigureReady: false,
  maxAbsLog2Fc: 0
};

const qlikSenseDiverging = [
  [0.0, "#2f4b7c"],
  [0.2, "#4f81bd"],
  [0.4, "#9bbcdf"],
  [0.5, "#f7f7f7"],
  [0.6, "#f7c6b8"],
  [0.8, "#e26d5a"],
  [1.0, "#a73030"]
];

const deseqFileInput = document.getElementById("deseqFile");
const geneFileInput = document.getElementById("geneFile");
const deseqFileName = document.getElementById("deseqFileName");
const geneFileName = document.getElementById("geneFileName");
const padjCutoffInput = document.getElementById("padjCutoff");
const log2fcCutoffInput = document.getElementById("log2FCCutoff");
const log2fcRangeText = document.getElementById("log2FCRange");
const updatePlotButton = document.getElementById("updatePlot");
const exportPdfButton = document.getElementById("exportPDF");
const plotInfo = document.getElementById("plotInfo");

function normalizeKey(key) {
  return String(key || "").trim().toLowerCase().replace(/[^a-z0-9]/g, "");
}

function normalizeGeneName(value) {
  return String(value || "")
    .trim()
    .replace(/^['"]+|['"]+$/g, "")
    .toUpperCase();
}

function toNumber(value) {
  if (value === null || value === undefined || value === "") {
    return NaN;
  }
  const n = Number(value);
  return Number.isFinite(n) ? n : NaN;
}

function setInfo(message, isError = false) {
  plotInfo.textContent = message;
  plotInfo.className = `plot-info ${isError ? "error" : "ok"}`;
}

function detectDelimiter(text) {
  const firstNonEmptyLine = text.split(/\r?\n/).find((line) => line.trim().length > 0) || "";
  const tabs = (firstNonEmptyLine.match(/\t/g) || []).length;
  const commas = (firstNonEmptyLine.match(/,/g) || []).length;
  return tabs > commas ? "\t" : ",";
}

function parseTable(text) {
  const delimiter = detectDelimiter(text);
  const parsed = Papa.parse(text, {
    header: true,
    delimiter,
    skipEmptyLines: true,
    dynamicTyping: false
  });
  if (parsed.errors.length > 0) {
    throw new Error(`File parse error: ${parsed.errors[0].message}`);
  }
  return parsed.data;
}

function getColumnNameByAliases(row, aliases) {
  const keys = Object.keys(row);
  const normalizedMap = {};
  keys.forEach((k) => {
    normalizedMap[normalizeKey(k)] = k;
  });
  for (const alias of aliases) {
    const hit = normalizedMap[normalizeKey(alias)];
    if (hit) {
      return hit;
    }
  }
  return null;
}

function getExpressionColumnsBetweenGeneAndBaseMean(row) {
  const columns = Object.keys(row);
  const normalizedColumns = columns.map(normalizeKey);
  const geneIdx = normalizedColumns.findIndex((col) => col === "genename");
  const baseMeanIdx = normalizedColumns.findIndex((col) => col === "basemean");
  if (geneIdx < 0 || baseMeanIdx < 0 || baseMeanIdx <= geneIdx + 1) {
    return [];
  }
  return columns.slice(geneIdx + 1, baseMeanIdx);
}

function rowZScore(values) {
  const mean = values.reduce((sum, val) => sum + val, 0) / values.length;
  const variance = values.reduce((sum, val) => sum + (val - mean) ** 2, 0) / values.length;
  const sd = Math.sqrt(variance);
  if (sd === 0) {
    return values.map(() => 0);
  }
  return values.map((val) => (val - mean) / sd);
}

function euclideanDistance(a, b) {
  let sum = 0;
  for (let i = 0; i < a.length; i += 1) {
    const d = a[i] - b[i];
    sum += d * d;
  }
  return Math.sqrt(sum);
}

function nearestNeighborOrder(vectors) {
  const n = vectors.length;
  if (n <= 1) {
    return [...Array(n).keys()];
  }
  let start = 0;
  let maxNorm = -Infinity;
  for (let i = 0; i < n; i += 1) {
    const norm = vectors[i].reduce((acc, val) => acc + val * val, 0);
    if (norm > maxNorm) {
      maxNorm = norm;
      start = i;
    }
  }
  const used = new Array(n).fill(false);
  const order = [start];
  used[start] = true;

  while (order.length < n) {
    const last = order[order.length - 1];
    let bestIdx = -1;
    let bestDist = Infinity;
    for (let i = 0; i < n; i += 1) {
      if (used[i]) {
        continue;
      }
      const dist = euclideanDistance(vectors[last], vectors[i]);
      if (dist < bestDist) {
        bestDist = dist;
        bestIdx = i;
      }
    }
    if (bestIdx < 0) {
      break;
    }
    used[bestIdx] = true;
    order.push(bestIdx);
  }
  return order;
}

function transpose(matrix) {
  if (matrix.length === 0) {
    return [];
  }
  return matrix[0].map((_, colIdx) => matrix.map((row) => row[colIdx]));
}

async function readGeneListSet(file) {
  if (!file) {
    return null;
  }
  const text = await file.text();
  const lines = text.split(/\r?\n/).map((x) => x.trim()).filter(Boolean);
  if (lines.length === 0) {
    return null;
  }

  if (lines[0].includes(",") || lines[0].includes("\t")) {
    const table = parseTable(text);
    if (table.length === 0) {
      return null;
    }
    const firstCol = Object.keys(table[0])[0];
    return new Set(
      table
        .map((row) => normalizeGeneName(row[firstCol]))
        .filter(Boolean)
    );
  }
  return new Set(lines.map((line) => normalizeGeneName(line)).filter(Boolean));
}

function updateLog2FcRangeDisplay(maxAbsLog2Fc) {
  state.maxAbsLog2Fc = maxAbsLog2Fc;
  log2fcCutoffInput.max = String(maxAbsLog2Fc);
  log2fcRangeText.textContent = `(0 - ${maxAbsLog2Fc.toFixed(4)})`;
  if (toNumber(log2fcCutoffInput.value) > maxAbsLog2Fc) {
    log2fcCutoffInput.value = String(maxAbsLog2Fc);
  }
}

function validateCutoffs() {
  const padjCutoff = toNumber(padjCutoffInput.value);
  if (!Number.isFinite(padjCutoff) || padjCutoff < 0 || padjCutoff > 1) {
    throw new Error("padj.cutoff must be between 0 and 1.");
  }

  const log2fcCutoff = toNumber(log2fcCutoffInput.value);
  if (!Number.isFinite(log2fcCutoff) || log2fcCutoff < 0 || log2fcCutoff > state.maxAbsLog2Fc) {
    throw new Error(`log2FC.cutoff must be between 0 and ${state.maxAbsLog2Fc.toFixed(4)}.`);
  }

  return { padjCutoff, log2fcCutoff };
}

async function parseDESeq2File() {
  const file = deseqFileInput.files[0];
  if (!file) {
    throw new Error("Please upload a DESeq2 result file.");
  }

  const text = await file.text();
  const rows = parseTable(text);
  if (rows.length === 0) {
    throw new Error("DESeq2 file has no usable rows.");
  }

  const geneCol = getColumnNameByAliases(rows[0], ["gene_name"]);
  const padjCol = getColumnNameByAliases(rows[0], ["padj", "adjustedpvalue", "fdr"]);
  const lfcCol = getColumnNameByAliases(rows[0], ["log2FoldChange", "log2FC", "lfc"]);
  const baseMeanCol = getColumnNameByAliases(rows[0], ["baseMean"]);
  if (!geneCol || !padjCol || !lfcCol || !baseMeanCol) {
    throw new Error("DESeq2 file must contain gene_name, padj, log2FoldChange, and baseMean columns.");
  }

  const expressionCols = getExpressionColumnsBetweenGeneAndBaseMean(rows[0]);
  if (expressionCols.length === 0) {
    throw new Error("No expression columns found between gene_name and baseMean.");
  }

  let maxAbsLog2Fc = 0;
  rows.forEach((row) => {
    const v = Math.abs(toNumber(row[lfcCol]));
    if (Number.isFinite(v) && v > maxAbsLog2Fc) {
      maxAbsLog2Fc = v;
    }
  });
  if (!Number.isFinite(maxAbsLog2Fc) || maxAbsLog2Fc <= 0) {
    throw new Error("Could not determine maximum |log2FoldChange| from DESeq2 file.");
  }

  updateLog2FcRangeDisplay(maxAbsLog2Fc);

  state.deseqRows = rows;
  return { rows, geneCol, padjCol, lfcCol, expressionCols };
}

function buildHeatmapMatrix(parsed, cutoffs, selectedGenesSet) {
  const { rows, geneCol, padjCol, lfcCol, expressionCols } = parsed;
  const bestRowByGene = new Map();

  rows.forEach((row) => {
    const geneRaw = String(row[geneCol] || "").trim();
    const geneKey = normalizeGeneName(geneRaw);
    if (!geneRaw || !geneKey) {
      return;
    }
    const lfc = toNumber(row[lfcCol]);
    if (!Number.isFinite(lfc)) {
      return;
    }
    const prev = bestRowByGene.get(geneKey);
    if (!prev || lfc > prev.lfc) {
      bestRowByGene.set(geneKey, { row, lfc, geneRaw });
    }
  });

  const dedupedRows = [...bestRowByGene.values()];
  const overlappedRows = selectedGenesSet
    ? dedupedRows.filter((item) => selectedGenesSet.has(normalizeGeneName(item.geneRaw)))
    : dedupedRows;

  if (selectedGenesSet && overlappedRows.length === 0) {
    throw new Error("No overlap between uploaded gene list and DESeq2 gene_name values.");
  }

  const eligibleGenes = new Set();
  overlappedRows.forEach((item) => {
    const row = item.row;
    const gene = normalizeGeneName(item.geneRaw);
    if (!gene) {
      return;
    }
    const padj = toNumber(row[padjCol]);
    const lfc = toNumber(row[lfcCol]);
    if (!Number.isFinite(padj) || !Number.isFinite(lfc)) {
      return;
    }
    if (padj <= cutoffs.padjCutoff && Math.abs(lfc) >= cutoffs.log2fcCutoff) {
      eligibleGenes.add(gene);
    }
  });

  const finalGenes = eligibleGenes;

  if (finalGenes.size === 0) {
    throw new Error("No genes passed filtering criteria.");
  }

  const geneNames = [];
  const matrix = [];
  overlappedRows.forEach((item) => {
    const row = item.row;
    const geneRaw = String(item.geneRaw || "").trim();
    const gene = normalizeGeneName(geneRaw);
    if (!finalGenes.has(gene)) {
      return;
    }
    const values = expressionCols.map((col) => toNumber(row[col]));
    if (values.some((v) => !Number.isFinite(v))) {
      return;
    }
    geneNames.push(geneRaw);
    matrix.push(rowZScore(values));
  });

  if (matrix.length === 0) {
    throw new Error("Filtered genes were found but none had fully numeric expression values.");
  }

  return { geneNames, sampleNames: expressionCols, matrix, eligibleCount: eligibleGenes.size, finalCount: finalGenes.size };
}

function clusterAndReorder(data) {
  const rowOrder = nearestNeighborOrder(data.matrix);
  const rowClusteredMatrix = rowOrder.map((i) => data.matrix[i]);
  const rowClusteredGenes = rowOrder.map((i) => data.geneNames[i]);

  const colVectors = transpose(rowClusteredMatrix);
  const colOrder = nearestNeighborOrder(colVectors);
  const colClusteredSamples = colOrder.map((i) => data.sampleNames[i]);
  const clusteredMatrix = rowClusteredMatrix.map((row) => colOrder.map((i) => row[i]));

  return {
    matrix: clusteredMatrix,
    geneNames: rowClusteredGenes,
    sampleNames: colClusteredSamples
  };
}

async function renderHeatmap() {
  try {
    exportPdfButton.disabled = true;
    setInfo("Processing input files...");

    const parsed = await parseDESeq2File();
    const cutoffs = validateCutoffs();
    const selectedGenesSet = await readGeneListSet(geneFileInput.files[0]);

    const matrixData = buildHeatmapMatrix(parsed, cutoffs, selectedGenesSet);
    const clusteredData = clusterAndReorder(matrixData);
    const showAllGeneLabels = clusteredData.geneNames.length < 50;

    const figureData = [{
      type: "heatmap",
      z: clusteredData.matrix,
      y: clusteredData.geneNames,
      x: clusteredData.sampleNames,
      colorscale: qlikSenseDiverging,
      zmid: 0,
      colorbar: { title: "Row z-score" },
      hovertemplate: "Gene: %{y}<br>Sample: %{x}<br>z-score: %{z:.3f}<extra></extra>"
    }];

    const figureLayout = {
      title: `DESeq2 Heatmap (${clusteredData.geneNames.length} genes x ${clusteredData.sampleNames.length} samples)`,
      xaxis: { side: "top", tickangle: -45 },
      yaxis: { automargin: true, showticklabels: showAllGeneLabels },
      margin: { l: 130, r: 40, t: 90, b: 80 },
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#ffffff",
      height: Math.max(620, 260 + clusteredData.geneNames.length * 11)
    };

    await Plotly.newPlot("heatmapPlot", figureData, figureLayout, { responsive: true, displaylogo: false });
    state.lastFigureReady = true;
    exportPdfButton.disabled = false;
    setInfo(
      //`Heatmap ready. Eligible genes: ${matrixData.eligibleCount}. Plotted genes: ${clusteredData.geneNames.length}.`
    );
  } catch (err) {
    state.lastFigureReady = false;
    exportPdfButton.disabled = true;
    setInfo(err.message || "Failed to render heatmap.", true);
  }
}

async function exportHeatmapPdf() {
  if (!state.lastFigureReady) {
    setInfo("Generate a heatmap before exporting.", true);
    return;
  }
  try {
    setInfo("Exporting heatmap to PDF...");
    const imageData = await Plotly.toImage("heatmapPlot", {
      format: "png",
      width: 1800,
      height: 1200,
      scale: 2
    });
    const { jsPDF } = window.jspdf;
    const pdf = new jsPDF({ orientation: "landscape", unit: "pt", format: "a4" });
    const pageWidth = pdf.internal.pageSize.getWidth();
    const pageHeight = pdf.internal.pageSize.getHeight();
    const targetWidth = pageWidth - 40;
    const targetHeight = pageHeight - 40;
    pdf.addImage(imageData, "PNG", 20, 20, targetWidth, targetHeight, undefined, "FAST");
    pdf.save("deseq2_heatmap.pdf");
    setInfo("Export completed: deseq2_heatmap.pdf");
  } catch (err) {
    setInfo(`PDF export failed: ${err.message || "Unknown error"}`, true);
  }
}

async function refreshMaxRangeFromFile() {
  if (!deseqFileInput.files[0]) {
    return;
  }
  try {
    const parsed = await parseDESeq2File();
    setInfo(`Loaded ${parsed.rows.length} DESeq2 rows. Ready to generate heatmap.`);
  } catch (err) {
    state.lastFigureReady = false;
    exportPdfButton.disabled = true;
    setInfo(err.message || "Could not parse DESeq2 file.", true);
  }
}

function bindFileNamePreview(inputEl, outputEl) {
  inputEl.addEventListener("change", () => {
    const file = inputEl.files[0];
    outputEl.textContent = file ? file.name : "";
  });
}

bindFileNamePreview(deseqFileInput, deseqFileName);
bindFileNamePreview(geneFileInput, geneFileName);
deseqFileInput.addEventListener("change", refreshMaxRangeFromFile);
updatePlotButton.addEventListener("click", renderHeatmap);
exportPdfButton.addEventListener("click", exportHeatmapPdf);
