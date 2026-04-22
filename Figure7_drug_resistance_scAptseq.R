#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Figure 7 and related supplementary figures: drug-resistance reversal scApt-seq
# Manuscript: Molecular Pattern Recognition by Sequencing Multiplex Aptamers
#             for Cell Phenotyping at Single-cell Resolution
#
# Purpose
#   Reproduce the analyses corresponding to Figure 7 and the related
#   supplementary figures for the colorectal cancer drug-resistance reversal
#   model. The manuscript compares parental HCT-8 cells, 5-FU-resistant cells,
#   and evodiamine-treated 5-FU-resistant cells using paired transcriptome and
#   aptamer-derived single-cell readouts.
#
# Main outputs
#   Figure 7a: transcriptome-based UMAP by RNA clusters and experimental group
#   Figure 7b: aptamer-based UMAP by APT clusters and experimental group
#   Figure 7c: heatmap of mean aptamer signals across HCT-8, 5-FU, 5-FU-Evo
#   Figure 7d: CD49C/ZAJ2c RNA, aptamer, and optional bulk-protein comparison
#   Figure 7e: MET/CLN0003 RNA, aptamer, and optional bulk-protein comparison
#   Figure 7f: CD71/HG19 RNA, aptamer, and optional bulk-protein comparison
#
# Supplementary outputs
#   Figure S10: QC violin plots for the three experimental groups
#   Figure S11: UMAP feature plots for target mRNAs and aptamer signals
#   Figure S12: cluster-level RNA-aptamer correlation plots
#   Figure S13: target transcript box plots across groups
#   Figure S14: selected aptamer-signal box plots across groups
#
# Example
#   DRUG_RNA_DIR=/path/to/10Aptamer_naiyao_auto_analysis_result \
#   DRUG_APT_DIR=/path/to/10Aptamer_naiyao_auto_analysis_result \
#   DRUG_PROTEIN_TABLE=/path/to/bulk_TMT_target_protein_values.csv \
#   FIG7_OUT_DIR=Figure7_drug_resistance_outputs \
#   SAVE_PNG=true \
#   Rscript Figure7_drug_resistance_scAptseq_GitHub.R
#
# Notes
#   1. The script avoids hard-coded absolute paths and writes all outputs to
#      FIG7_OUT_DIR.
#   2. Raw aptamer counts are CLR-normalized and interpreted comparatively across
#      cells/conditions for the same aptamer, consistent with the revised text.
#   3. If package versions or input data differ from the submitted analysis,
#      clustering resolutions may be adjusted using environment variables.
# -----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)
set.seed(20260422)

# ----------------------------- package management ----------------------------
required_pkgs <- c(
  "Seurat", "dplyr", "tidyr", "ggplot2", "patchwork", "data.table",
  "Matrix", "pheatmap", "RColorBrewer"
)
optional_pkgs <- c("harmony", "DoubletFinder", "ggsignif")

missing_required <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_required, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
  library(pheatmap)
  library(RColorBrewer)
})

HAS_HARMONY <- requireNamespace("harmony", quietly = TRUE)
HAS_DF <- requireNamespace("DoubletFinder", quietly = TRUE)
HAS_GGSIGNIF <- requireNamespace("ggsignif", quietly = TRUE)

# ------------------------------- configuration -------------------------------
get_env <- function(name, default = NULL) {
  x <- Sys.getenv(name)
  if (identical(x, "")) default else x
}

RNA_DIR <- get_env("DRUG_RNA_DIR", get_env("NAIYAO_RNA_DIR", "./10Aptamer_naiyao_auto_analysis_result"))
APT_DIR <- get_env("DRUG_APT_DIR", get_env("NAIYAO_APT_DIR", RNA_DIR))
PROTEIN_TABLE <- get_env("DRUG_PROTEIN_TABLE", get_env("NAIYAO_PROTEIN_TABLE", ""))
OUT_DIR <- get_env("FIG7_OUT_DIR", get_env("DRUG_OUT_DIR", "./Figure7_drug_resistance_outputs"))
SAVE_PNG <- tolower(get_env("SAVE_PNG", "false")) %in% c("true", "t", "1", "yes", "y")
RUN_DOUBLET <- tolower(get_env("RUN_DOUBLET", "false")) %in% c("true", "t", "1", "yes", "y")

MIN_FEATURE_RNA <- as.numeric(get_env("MIN_FEATURE_RNA", "400"))
MAX_FEATURE_RNA <- as.numeric(get_env("MAX_FEATURE_RNA", "4500"))
MAX_PERCENT_MT <- as.numeric(get_env("MAX_PERCENT_MT", "10"))
N_PCS <- as.integer(get_env("N_PCS", "30"))
RNA_CLUSTER_RES <- as.numeric(get_env("RNA_CLUSTER_RES", "0.5"))
APT_CLUSTER_RES <- as.numeric(get_env("APT_CLUSTER_RES", "0.3"))
APT_PCS <- as.integer(get_env("APT_PCS", "8"))
CORRELATION_GROUP_BY <- get_env("CORRELATION_GROUP_BY", "rna_cluster")

if (RUN_DOUBLET && !HAS_DF) {
  warning("RUN_DOUBLET=true, but DoubletFinder is not installed. Doublet filtering will be skipped.")
  RUN_DOUBLET <- FALSE
}

# ------------------------------- output paths --------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIG_DIR <- file.path(OUT_DIR, "Figure7")
SUPP_DIR <- file.path(OUT_DIR, "Supplementary")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDS_DIR <- file.path(OUT_DIR, "rds")
for (d in c(FIG_DIR, SUPP_DIR, TAB_DIR, RDS_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot, name, width = 8, height = 6, dir = FIG_DIR) {
  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, units = "in", bg = "white")
  if (SAVE_PNG) {
    png_path <- file.path(dir, paste0(name, ".png"))
    ggplot2::ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  }
}

save_pheatmap <- function(ph, name, width = 6, height = 7, dir = FIG_DIR) {
  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  grDevices::pdf(pdf_path, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  grDevices::dev.off()
  if (SAVE_PNG) {
    png_path <- file.path(dir, paste0(name, ".png"))
    grDevices::png(png_path, width = width, height = height, units = "in", res = 300)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    grDevices::dev.off()
  }
}

# ------------------------------ metadata maps --------------------------------
group_levels <- c("HCT-8", "5-FU", "5-FU-Evo")
group_colors <- c("HCT-8" = "#4B727E", "5-FU" = "#70B199", "5-FU-Evo" = "#CDE6CC")

# The figure legend uses common protein names, whereas RNA data use gene symbols.
target_pairs <- data.frame(
  display_name = c("PTK7", "MET", "CD71", "PD-L1", "CD49C", "PTPRF", "CD318", "EpCAM", "ZNRF3"),
  gene = c("PTK7", "MET", "TFRC", "CD274", "ITGA3", "PTPRF", "CDCP1", "EPCAM", "ZNRF3"),
  aptamer_feature = c("PTK7-APT", "MET-APT", "CD71-APT", "PD-L1-APT", "CD49C-APT", "PTPRF-APT", "CD318-APT", "EPCAM-APT", "ZNRF3-APT"),
  aptamer_name = c("sgc8c", "CLN0003", "HG19", "PD-L1", "ZAJ2c", "XQ20a", "XQ3b", "SYL3C", "ZNRF3"),
  stringsAsFactors = FALSE
)

# Figure 7c order follows the biological interpretation described in the text.
fig7c_aptamer_order <- c(
  "CD49C-APT", "PTK7-APT",          # enriched in HCT-8
  "PD-L1-APT", "EPCAM-APT", "PTPRF-APT", "MET-APT", "ZNRF3-APT", # enriched in 5-FU
  "CD318-APT", "CD71-APT"           # enriched in 5-FU-Evo
)
fig7c_row_labels <- target_pairs$aptamer_name[match(fig7c_aptamer_order, target_pairs$aptamer_feature)]

# Panels Figure 7d-f.
fig7_detail_pairs <- data.frame(
  panel = c("Figure7d_CD49C_ZAJ2c", "Figure7e_MET_CLN0003", "Figure7f_CD71_HG19"),
  display_name = c("CD49C", "MET", "CD71"),
  gene = c("ITGA3", "MET", "TFRC"),
  aptamer_feature = c("CD49C-APT", "MET-APT", "CD71-APT"),
  aptamer_name = c("ZAJ2c", "CLN0003", "HG19"),
  stringsAsFactors = FALSE
)

# ------------------------------- utilities -----------------------------------
message_header <- function(txt) {
  message("\n", paste(rep("=", 78), collapse = ""))
  message(txt)
  message(paste(rep("=", 78), collapse = ""))
}

clean_sample_id <- function(x) {
  x <- basename(x)
  x <- sub("/+$", "", x)
  x <- sub("_P_cDNA$", "", x)
  x <- sub("_cDNA$", "", x)
  x
}

infer_group <- function(sample_id) {
  sid <- clean_sample_id(sample_id)
  out <- rep(NA_character_, length(sid))
  out[grepl("^HCT[-_]?8|^HCT8|parent", sid, ignore.case = TRUE)] <- "HCT-8"
  out[grepl("5[-_]?FU|FU", sid, ignore.case = TRUE) & !grepl("EVO", sid, ignore.case = TRUE)] <- "5-FU"
  out[grepl("EVO|Evo|evodiamine", sid, ignore.case = TRUE)] <- "5-FU-Evo"
  if (any(is.na(out))) {
    warning("Could not infer experimental group for sample(s): ", paste(sid[is.na(out)], collapse = ", "))
  }
  factor(out, levels = group_levels)
}

find_matrix_dir <- function(sample_dir) {
  candidates <- c(
    file.path(sample_dir, "filter_matrix"),
    file.path(sample_dir, "filtered_feature_bc_matrix"),
    file.path(sample_dir, "filtered_gene_bc_matrices"),
    sample_dir
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) stop("No matrix directory found for ", sample_dir)
  # Prefer directories containing matrix.mtx or matrix.mtx.gz.
  for (d in candidates) {
    if (file.exists(file.path(d, "matrix.mtx")) || file.exists(file.path(d, "matrix.mtx.gz"))) return(d)
  }
  candidates[1]
}

list_sample_dirs <- function(base_dir) {
  if (!dir.exists(base_dir)) stop("Input directory does not exist: ", base_dir)
  dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl("HCT|5[-_]?FU|EVO|_P_cDNA", basename(dirs), ignore.case = TRUE)]
  if (length(dirs) == 0) {
    dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  }
  if (length(dirs) == 0) stop("No sample directories found in ", base_dir)
  dirs
}

standardize_cell_barcode <- function(x) {
  x <- as.character(x)
  x <- sub("-1$", "", x)
  x <- gsub("\\s+", "", x)
  x
}

standardize_apt_feature <- function(x) {
  y <- trimws(as.character(x))
  y <- gsub("_", "-", y)
  y <- gsub("\\s+", "", y)
  y_upper <- toupper(y)
  out <- y
  out[grepl("SGC8|PTK7", y_upper)] <- "PTK7-APT"
  out[grepl("CLN0003|^MET", y_upper)] <- "MET-APT"
  out[grepl("HG19|CD71|TFRC", y_upper)] <- "CD71-APT"
  out[grepl("PD-L1|PDL1|CD274", y_upper)] <- "PD-L1-APT"
  out[grepl("ZAJ2C|CD49C|ITGA3", y_upper)] <- "CD49C-APT"
  out[grepl("XQ20A|PTPRF", y_upper)] <- "PTPRF-APT"
  out[grepl("XQ3B|CD318|CDCP1", y_upper)] <- "CD318-APT"
  out[grepl("SYL3C|EPCAM", y_upper)] <- "EPCAM-APT"
  out[grepl("ZNRF3", y_upper)] <- "ZNRF3-APT"
  out[grepl("CTRL|CONTROL|SCRAMBLE|SCRAMBLED", y_upper)] <- "Scramble-APT"
  out
}

safe_fetch <- function(obj, vars, assay = NULL, slot = "data") {
  if (!is.null(assay)) DefaultAssay(obj) <- assay
  keep <- vars[vars %in% c(rownames(obj), colnames(obj@meta.data))]
  missing <- setdiff(vars, keep)
  if (length(missing) > 0) warning("Feature(s) not found and skipped: ", paste(missing, collapse = ", "))
  if (length(keep) == 0) return(NULL)
  FetchData(obj, vars = keep, slot = slot)
}

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

# ----------------------------- read RNA matrices -----------------------------
message_header("Reading RNA count matrices")
rna_sample_dirs <- list_sample_dirs(RNA_DIR)
rna_sample_ids <- clean_sample_id(rna_sample_dirs)
names(rna_sample_dirs) <- rna_sample_ids

rna_objects <- list()
for (sid in names(rna_sample_dirs)) {
  matrix_dir <- find_matrix_dir(rna_sample_dirs[[sid]])
  message("Reading RNA sample: ", sid, " from ", matrix_dir)
  counts <- Read10X(data.dir = matrix_dir, gene.column = 1)
  if (is.list(counts)) {
    # Prefer Gene Expression if a multi-modal 10X-like matrix is supplied.
    if ("Gene Expression" %in% names(counts)) counts <- counts[["Gene Expression"]] else counts <- counts[[1]]
  }
  obj <- CreateSeuratObject(counts = counts, project = sid, min.cells = 1)
  obj <- RenameCells(obj, add.cell.id = sid)
  obj$sample_id <- sid
  obj$group <- infer_group(sid)
  rna_objects[[sid]] <- obj
}

if (length(rna_objects) == 1) {
  sce <- rna_objects[[1]]
} else {
  sce <- merge(x = rna_objects[[1]], y = rna_objects[-1], merge.data = TRUE)
}

sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce$sample_id <- as.character(sce$sample_id)
sce$group <- factor(as.character(sce$group), levels = group_levels)

# ------------------------ optional doublet filtering -------------------------
if (RUN_DOUBLET) {
  message_header("Running DoubletFinder")
  # A simple per-merged-object implementation. For strict sample-wise filtering,
  # run DoubletFinder before merging or adapt this block by sample.
  sce <- NormalizeData(sce, verbose = FALSE)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sce <- ScaleData(sce, verbose = FALSE)
  sce <- RunPCA(sce, npcs = 30, verbose = FALSE)
  sweep_res <- DoubletFinder::paramSweep_v3(sce, PCs = 1:20, sct = FALSE)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep_stats)
  pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  n_exp <- round(0.05 * ncol(sce))
  sce <- DoubletFinder::doubletFinder_v3(sce, PCs = 1:20, pN = 0.25, pK = pk, nExp = n_exp, reuse.pANN = FALSE, sct = FALSE)
  df_col <- tail(grep("DF.classifications", colnames(sce@meta.data), value = TRUE), 1)
  sce$doublet_info <- sce@meta.data[[df_col]]
  sce <- subset(sce, subset = doublet_info == "Singlet")
}

# ------------------------------ RNA QC filter --------------------------------
message_header("Applying RNA QC filters")
qc_before <- sce@meta.data %>%
  group_by(group, sample_id) %>%
  summarise(n_cells_before_qc = dplyr::n(), .groups = "drop")

sce <- subset(
  sce,
  subset = nFeature_RNA >= MIN_FEATURE_RNA &
    nFeature_RNA <= MAX_FEATURE_RNA &
    percent.mt <= MAX_PERCENT_MT
)

qc_after <- sce@meta.data %>%
  group_by(group, sample_id) %>%
  summarise(n_cells_after_qc = dplyr::n(), .groups = "drop")
qc_summary <- full_join(qc_before, qc_after, by = c("group", "sample_id"))
write.csv(qc_summary, file.path(TAB_DIR, "QC_cells_before_after_filtering.csv"), row.names = FALSE)

# -------------------------- transcriptome analysis ---------------------------
message_header("Running transcriptome analysis")
DefaultAssay(sce) <- "RNA"
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, features = VariableFeatures(sce), npcs = max(N_PCS, 10), verbose = FALSE)

if (HAS_HARMONY && length(unique(sce$sample_id)) > 1) {
  message("Using Harmony integration with group.by.vars = sample_id")
  sce <- harmony::RunHarmony(sce, group.by.vars = "sample_id", reduction = "pca", assay.use = "RNA", verbose = FALSE)
  rna_reduction <- "harmony"
} else {
  message("Harmony not available or single sample only; using PCA reduction.")
  rna_reduction <- "pca"
}

sce <- RunUMAP(sce, reduction = rna_reduction, dims = 1:N_PCS, reduction.name = "rna_umap", verbose = FALSE)
sce <- FindNeighbors(sce, reduction = rna_reduction, dims = 1:N_PCS, verbose = FALSE)
sce <- FindClusters(sce, algorithm = 2, resolution = RNA_CLUSTER_RES, verbose = FALSE)
sce$rna_cluster <- factor(sce$seurat_clusters)

# ---------------------------- read aptamer counts ----------------------------
message_header("Reading aptamer count matrices")
apt_sample_dirs <- list_sample_dirs(APT_DIR)
apt_sample_ids <- clean_sample_id(apt_sample_dirs)
names(apt_sample_dirs) <- apt_sample_ids

read_apt_cell_table <- function(sample_dir, sample_id) {
  candidates <- c(
    file.path(sample_dir, "ALL_APT_cell.csv"),
    file.path(sample_dir, "all_APT_cell.csv"),
    file.path(sample_dir, "aptamer_cell_counts.csv"),
    file.path(sample_dir, "APT_cell.csv")
  )
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f)) {
    warning("No ALL_APT_cell.csv-like file found for sample ", sample_id, " in ", sample_dir)
    return(NULL)
  }
  dt <- data.table::fread(f, data.table = FALSE)
  cell_col <- intersect(c("cell", "barcode", "Cell", "Barcode", "cell_id", "CellID"), colnames(dt))[1]
  if (is.na(cell_col)) {
    stop("No cell/barcode column found in aptamer table: ", f)
  }
  cell <- standardize_cell_barcode(dt[[cell_col]])
  mat_df <- dt[, setdiff(colnames(dt), cell_col), drop = FALSE]
  # Remove non-aptamer metadata columns if present.
  non_apt_cols <- c("sample", "group", "patient", "cellnumber", "sum", "total", "mean")
  mat_df <- mat_df[, setdiff(colnames(mat_df), non_apt_cols), drop = FALSE]
  colnames(mat_df) <- standardize_apt_feature(colnames(mat_df))
  # Collapse duplicate standardized columns by summing.
  unique_features <- unique(colnames(mat_df))
  mat_df2 <- data.frame(row.names = paste(sample_id, cell, sep = "_"))
  for (feat in unique_features) {
    cols <- which(colnames(mat_df) == feat)
    mat_df2[[feat]] <- rowSums(as.data.frame(mat_df[, cols, drop = FALSE]), na.rm = TRUE)
  }
  mat_df2$sample_id <- sample_id
  mat_df2
}

apt_tables <- lapply(names(apt_sample_dirs), function(sid) read_apt_cell_table(apt_sample_dirs[[sid]], sid))
apt_tables <- Filter(Negate(is.null), apt_tables)
if (length(apt_tables) == 0) stop("No cell-level aptamer tables were loaded from DRUG_APT_DIR: ", APT_DIR)

# Bind with missing columns filled with zero.
all_cols <- Reduce(union, lapply(apt_tables, colnames))
apt_tables <- lapply(apt_tables, function(df) {
  miss <- setdiff(all_cols, colnames(df))
  for (m in miss) df[[m]] <- 0
  df[, all_cols, drop = FALSE]
})
apt_df <- dplyr::bind_rows(apt_tables)
apt_features <- intersect(target_pairs$aptamer_feature, colnames(apt_df))
if (length(apt_features) < 3) {
  warning("Fewer than three canonical aptamer features detected. Available columns: ", paste(colnames(apt_df), collapse = ", "))
}

# Harmonize cell IDs between RNA and aptamer tables.
common_cells <- intersect(colnames(sce), rownames(apt_df))
if (length(common_cells) == 0) {
  stop(
    "No overlapping cell barcodes between RNA and aptamer data.\n",
    "Check whether sample prefixes in ALL_APT_cell.csv match RNA sample folder names."
  )
}
message("Matched ", length(common_cells), " cells between RNA and aptamer data.")
sce <- sce[, common_cells]
apt_df <- apt_df[common_cells, , drop = FALSE]

# Add APT assay with canonical feature order.
apt_features_ordered <- target_pairs$aptamer_feature[target_pairs$aptamer_feature %in% colnames(apt_df)]
if ("Scramble-APT" %in% colnames(apt_df)) apt_features_ordered <- c(apt_features_ordered, "Scramble-APT")
apt_mat <- as.matrix(apt_df[, apt_features_ordered, drop = FALSE])
apt_mat <- t(apt_mat)
apt_mat <- Matrix::Matrix(apt_mat, sparse = TRUE)
sce[["APT"]] <- CreateAssayObject(counts = apt_mat)
DefaultAssay(sce) <- "APT"
sce <- NormalizeData(sce, normalization.method = "CLR", margin = 2, verbose = FALSE)
sce$nCount_APT <- Matrix::colSums(GetAssayData(sce, assay = "APT", slot = "counts"))
sce$nFeature_APT <- Matrix::colSums(GetAssayData(sce, assay = "APT", slot = "counts") > 0)

# Save quality-control summary after aptamer merge.
qc_metrics <- sce@meta.data %>%
  group_by(group) %>%
  summarise(
    cell_number = dplyr::n(),
    median_gene_number = median(nFeature_RNA, na.rm = TRUE),
    median_umi = median(nCount_RNA, na.rm = TRUE),
    median_mitochondrial_percentage = median(percent.mt, na.rm = TRUE),
    median_aptamer_counts = median(nCount_APT, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(qc_metrics, file.path(TAB_DIR, "FigureS10_QC_summary_by_group.csv"), row.names = FALSE)

# ---------------------------- aptamer clustering -----------------------------
message_header("Running aptamer-based clustering")
DefaultAssay(sce) <- "APT"
sce <- ScaleData(sce, verbose = FALSE)
# Use all aptamer features for APT PCA; number of PCs cannot exceed number of features.
n_apt_pcs <- min(APT_PCS, length(rownames(sce[["APT"]])) - 1)
if (n_apt_pcs < 2) n_apt_pcs <- 2
sce <- RunPCA(sce, features = rownames(sce[["APT"]]), npcs = n_apt_pcs, reduction.name = "apt_pca", verbose = FALSE)
sce <- RunUMAP(sce, reduction = "apt_pca", dims = 1:n_apt_pcs, reduction.name = "apt_umap", verbose = FALSE)
sce <- FindNeighbors(sce, reduction = "apt_pca", dims = 1:n_apt_pcs, graph.name = "apt_snn", verbose = FALSE)
sce <- FindClusters(sce, graph.name = "apt_snn", algorithm = 2, resolution = APT_CLUSTER_RES, verbose = FALSE)
apt_cluster_col <- paste0("apt_snn_res.", APT_CLUSTER_RES)
if (!apt_cluster_col %in% colnames(sce@meta.data)) {
  apt_cluster_col <- tail(grep("apt_snn_res", colnames(sce@meta.data), value = TRUE), 1)
}
sce$apt_cluster <- factor(sce@meta.data[[apt_cluster_col]])

# ---------------------------- Figure S10: QC ---------------------------------
message_header("Generating Figure S10 outputs")
qc_plot <- VlnPlot(
  sce,
  features = c("nCount_APT", "nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "group",
  cols = group_colors,
  pt.size = 0,
  ncol = 4
) & theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(qc_plot, "FigureS10_QC_violin_plots", width = 14, height = 4.5, dir = SUPP_DIR)

# ---------------------------- Figure 7a / 7b ---------------------------------
message_header("Generating Figure 7a and 7b outputs")
cluster_palette <- colorRampPalette(brewer.pal(12, "Paired"))(max(12, length(unique(sce$rna_cluster)), length(unique(sce$apt_cluster))))

fig7a_cluster <- DimPlot(sce, reduction = "rna_umap", group.by = "rna_cluster", label = TRUE, repel = TRUE, cols = cluster_palette) +
  ggtitle("Transcriptome clusters") + theme(plot.title = element_text(hjust = 0.5))
fig7a_group <- DimPlot(sce, reduction = "rna_umap", group.by = "group", cols = group_colors) +
  ggtitle("Experimental group") + theme(plot.title = element_text(hjust = 0.5))
save_plot(fig7a_cluster + fig7a_group, "Figure7a_transcriptome_UMAP", width = 10, height = 5)

fig7b_cluster <- DimPlot(sce, reduction = "apt_umap", group.by = "apt_cluster", label = TRUE, repel = TRUE, cols = cluster_palette) +
  ggtitle("Aptamer clusters") + theme(plot.title = element_text(hjust = 0.5))
fig7b_group <- DimPlot(sce, reduction = "apt_umap", group.by = "group", cols = group_colors) +
  ggtitle("Experimental group") + theme(plot.title = element_text(hjust = 0.5))
save_plot(fig7b_cluster + fig7b_group, "Figure7b_aptamer_UMAP", width = 10, height = 5)

# ---------------------------- Figure 7c heatmap ------------------------------
message_header("Generating Figure 7c heatmap")
DefaultAssay(sce) <- "APT"
apt_data <- GetAssayData(sce, assay = "APT", slot = "data")
apt_group_means <- t(apply(apt_data, 1, function(x) tapply(x, sce$group, mean, na.rm = TRUE)))
apt_group_means <- apt_group_means[intersect(fig7c_aptamer_order, rownames(apt_group_means)), group_levels, drop = FALSE]
apt_group_scaled <- zscore_rows(apt_group_means)
rownames(apt_group_scaled) <- target_pairs$aptamer_name[match(rownames(apt_group_scaled), target_pairs$aptamer_feature)]
write.csv(apt_group_means, file.path(TAB_DIR, "Figure7c_mean_CLR_aptamer_signal_by_group.csv"))

fig7c_heatmap <- pheatmap::pheatmap(
  apt_group_scaled,
  color = colorRampPalette(rev(brewer.pal(9, "GnBu")))(100),
  border_color = "grey75",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = data.frame(Group = factor(group_levels, levels = group_levels), row.names = group_levels),
  annotation_colors = list(Group = group_colors),
  main = "Mean aptamer signal",
  fontsize_row = 10,
  fontsize_col = 10,
  silent = TRUE
)
save_pheatmap(fig7c_heatmap, "Figure7c_mean_aptamer_signal_heatmap", width = 5, height = 6)

# ---------------------------- Figure S11 feature plots -----------------------
message_header("Generating Figure S11 feature plots")
DefaultAssay(sce) <- "RNA"
rna_features <- target_pairs$gene[target_pairs$gene %in% rownames(sce)]
rna_feature_plot <- FeaturePlot(
  sce,
  reduction = "rna_umap",
  features = rna_features,
  cols = c("lightgrey", "#545493"),
  pt.size = 0.35,
  ncol = 3,
  combine = TRUE
) & theme(plot.title = element_text(hjust = 0.5))
save_plot(rna_feature_plot, "FigureS11_target_mRNA_UMAP_feature_plots", width = 10, height = 10, dir = SUPP_DIR)

DefaultAssay(sce) <- "APT"
apt_features_present <- target_pairs$aptamer_feature[target_pairs$aptamer_feature %in% rownames(sce[["APT"]])]
apt_feature_plot <- FeaturePlot(
  sce,
  reduction = "rna_umap",
  features = apt_features_present,
  cols = c("lightgrey", "#3F693B"),
  pt.size = 0.35,
  ncol = 3,
  combine = TRUE
) & theme(plot.title = element_text(hjust = 0.5))
save_plot(apt_feature_plot, "FigureS11_aptamer_UMAP_feature_plots", width = 10, height = 10, dir = SUPP_DIR)

# ---------------------------- Figure S12 correlations ------------------------
message_header("Generating Figure S12 correlation plots")
# Use RNA clusters by default, matching the manuscript's cluster-level comparison.
if (!CORRELATION_GROUP_BY %in% colnames(sce@meta.data)) {
  warning("CORRELATION_GROUP_BY=", CORRELATION_GROUP_BY, " not found. Using rna_cluster.")
  CORRELATION_GROUP_BY <- "rna_cluster"
}
DefaultAssay(sce) <- "RNA"
rna_avg <- AverageExpression(sce, assays = "RNA", group.by = CORRELATION_GROUP_BY, slot = "data", verbose = FALSE)$RNA
DefaultAssay(sce) <- "APT"
apt_avg <- AverageExpression(sce, assays = "APT", group.by = CORRELATION_GROUP_BY, slot = "data", verbose = FALSE)$APT

cor_plots <- list()
cor_summary <- data.frame()
for (i in seq_len(nrow(target_pairs))) {
  gene <- target_pairs$gene[i]
  apt <- target_pairs$aptamer_feature[i]
  label <- paste0(target_pairs$display_name[i], " / ", target_pairs$aptamer_name[i])
  if (!gene %in% rownames(rna_avg) || !apt %in% rownames(apt_avg)) next
  df <- data.frame(
    RNA = as.numeric(rna_avg[gene, ]),
    APT = as.numeric(apt_avg[apt, ]),
    Cluster = colnames(rna_avg)
  )
  if (nrow(df) >= 3 && stats::sd(df$RNA) > 0 && stats::sd(df$APT) > 0) {
    ct <- suppressWarnings(cor.test(df$RNA, df$APT, method = "pearson"))
    r_val <- as.numeric(ct$estimate)
    p_val <- ct$p.value
  } else {
    r_val <- NA_real_
    p_val <- NA_real_
  }
  cor_summary <- rbind(cor_summary, data.frame(
    display_name = target_pairs$display_name[i],
    gene = gene,
    aptamer = target_pairs$aptamer_name[i],
    aptamer_feature = apt,
    pearson_r = r_val,
    p_value = p_val,
    stringsAsFactors = FALSE
  ))
  cor_plots[[label]] <- ggplot(df, aes(x = RNA, y = APT)) +
    geom_point(size = 2.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70", linewidth = 0.6) +
    labs(
      title = label,
      subtitle = paste0("r = ", ifelse(is.na(r_val), "NA", sprintf("%.3f", r_val))),
      x = paste0(gene, " RNA"),
      y = paste0(target_pairs$aptamer_name[i], " aptamer")
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
}
write.csv(cor_summary, file.path(TAB_DIR, "FigureS12_cluster_level_RNA_aptamer_correlations.csv"), row.names = FALSE)
if (length(cor_plots) > 0) {
  cor_grid <- wrap_plots(cor_plots, ncol = 3) + plot_annotation(title = "Cluster-level RNA-aptamer correlations")
  save_plot(cor_grid, "FigureS12_cluster_level_RNA_aptamer_correlations", width = 12, height = 10, dir = SUPP_DIR)
}

# ---------------------------- boxplot helpers --------------------------------
make_boxplot_from_seurat <- function(obj, feature, assay, title, ylab = "Expression") {
  DefaultAssay(obj) <- assay
  df <- FetchData(obj, vars = c(feature, "group"), slot = "data")
  colnames(df) <- c("Expression", "Group")
  df$Group <- factor(df$Group, levels = group_levels)
  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, color = "#2C3E50", linewidth = 0.45, alpha = 0.95) +
    scale_fill_manual(values = group_colors) +
    labs(title = title, x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  if (HAS_GGSIGNIF) {
    y_max <- max(df$Expression, na.rm = TRUE)
    if (is.finite(y_max) && y_max > 0) {
      comparisons <- list(c("HCT-8", "5-FU"), c("5-FU", "5-FU-Evo"), c("HCT-8", "5-FU-Evo"))
      p <- p + ggsignif::geom_signif(
        comparisons = comparisons,
        test = "wilcox.test",
        map_signif_level = TRUE,
        y_position = y_max * c(1.05, 1.17, 1.29),
        tip_length = 0.01,
        textsize = 3
      ) + coord_cartesian(ylim = c(0, y_max * 1.35))
    }
  }
  p
}

# ---------------------------- Figure S13 / S14 -------------------------------
message_header("Generating Figure S13 and S14 box plots")
gene_boxplots <- list()
DefaultAssay(sce) <- "RNA"
for (i in seq_len(nrow(target_pairs))) {
  gene <- target_pairs$gene[i]
  if (gene %in% rownames(sce)) {
    gene_boxplots[[target_pairs$display_name[i]]] <- make_boxplot_from_seurat(
      sce, gene, "RNA", paste0(target_pairs$display_name[i], " RNA"), ylab = "RNA expression"
    )
  }
}
if (length(gene_boxplots) > 0) {
  s13 <- wrap_plots(gene_boxplots, ncol = 3) + plot_annotation(title = "Target transcript expression across groups")
  save_plot(s13, "FigureS13_target_transcript_boxplots", width = 12, height = 10, dir = SUPP_DIR)
}

apt_boxplots <- list()
DefaultAssay(sce) <- "APT"
for (i in seq_len(nrow(target_pairs))) {
  apt <- target_pairs$aptamer_feature[i]
  if (apt %in% rownames(sce[["APT"]])) {
    apt_boxplots[[target_pairs$aptamer_name[i]]] <- make_boxplot_from_seurat(
      sce, apt, "APT", target_pairs$aptamer_name[i], ylab = "CLR-transformed aptamer signal"
    )
  }
}
if (length(apt_boxplots) > 0) {
  s14 <- wrap_plots(apt_boxplots, ncol = 3) + plot_annotation(title = "Selected aptamer signals across groups")
  save_plot(s14, "FigureS14_aptamer_signal_boxplots", width = 12, height = 10, dir = SUPP_DIR)
}

# Export raw long-format expression summaries for reproducibility.
export_expression_by_group <- function(obj, features, assay, output_name) {
  DefaultAssay(obj) <- assay
  features <- features[features %in% rownames(obj[[assay]])]
  if (length(features) == 0) return(NULL)
  df <- FetchData(obj, vars = c(features, "group"), slot = "data")
  out <- tidyr::pivot_longer(df, cols = all_of(features), names_to = "feature", values_to = "value") %>%
    group_by(group, feature) %>%
    summarise(
      n_cells = dplyr::n(),
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(out, file.path(TAB_DIR, output_name), row.names = FALSE)
}
export_expression_by_group(sce, target_pairs$gene, "RNA", "FigureS13_target_transcript_summary_by_group.csv")
export_expression_by_group(sce, target_pairs$aptamer_feature, "APT", "FigureS14_aptamer_signal_summary_by_group.csv")

# ---------------------------- optional protein table -------------------------
read_protein_table <- function(path) {
  if (is.null(path) || path == "" || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  protein <- if (ext %in% c("tsv", "txt")) {
    data.table::fread(path, sep = "\t", data.table = FALSE)
  } else {
    data.table::fread(path, data.table = FALSE)
  }
  colnames(protein) <- trimws(colnames(protein))
  id_col <- intersect(tolower(colnames(protein)), c("target", "protein", "gene", "display_name", "marker"))[1]
  if (is.na(id_col)) {
    id_col <- tolower(colnames(protein))[1]
  }
  id_col_real <- colnames(protein)[match(id_col, tolower(colnames(protein)))]
  group_col <- intersect(tolower(colnames(protein)), c("group", "condition", "sample_group"))[1]
  value_col <- intersect(tolower(colnames(protein)), c("value", "abundance", "protein", "intensity", "expression"))[1]
  if (!is.na(group_col) && !is.na(value_col)) {
    group_col_real <- colnames(protein)[match(group_col, tolower(colnames(protein)))]
    value_col_real <- colnames(protein)[match(value_col, tolower(colnames(protein)))]
    out <- protein %>%
      rename(target = all_of(id_col_real), group = all_of(group_col_real), value = all_of(value_col_real)) %>%
      mutate(target = as.character(target), group = as.character(group), value = as.numeric(value))
  } else {
    # Wide table: target plus group columns.
    out <- protein %>%
      rename(target = all_of(id_col_real)) %>%
      pivot_longer(cols = -target, names_to = "group", values_to = "value") %>%
      mutate(target = as.character(target), group = as.character(group), value = as.numeric(value))
  }
  group_chr <- as.character(out$group)
  group_norm <- rep(NA_character_, length(group_chr))
  group_norm[grepl("^HCT[-_ ]?8|^HCT8|parent", group_chr, ignore.case = TRUE)] <- "HCT-8"
  group_norm[grepl("5[-_ ]?FU|FU", group_chr, ignore.case = TRUE) & !grepl("EVO|Evo|evodiamine", group_chr, ignore.case = TRUE)] <- "5-FU"
  group_norm[grepl("EVO|Evo|evodiamine", group_chr, ignore.case = TRUE)] <- "5-FU-Evo"
  group_norm[is.na(group_norm)] <- group_chr[is.na(group_norm)]
  out$group <- factor(group_norm, levels = group_levels)
  out
}

protein_df <- read_protein_table(PROTEIN_TABLE)
if (!is.null(protein_df)) {
  write.csv(protein_df, file.path(TAB_DIR, "bulk_TMT_protein_values_used_for_Figure7d_f.csv"), row.names = FALSE)
} else {
  message("No DRUG_PROTEIN_TABLE supplied. Figure 7d-f will include RNA and aptamer panels only; protein panels will be skipped.")
}

make_protein_boxplot <- function(protein_df, display_name) {
  if (is.null(protein_df)) return(NULL)
  alias <- unique(c(display_name, target_pairs$gene[target_pairs$display_name == display_name]))
  df <- protein_df %>% filter(toupper(target) %in% toupper(alias), !is.na(group), is.finite(value))
  if (nrow(df) == 0) return(NULL)
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, color = "#2C3E50", linewidth = 0.45, alpha = 0.95) +
    geom_point(position = position_jitter(width = 0.08), size = 1.3, alpha = 0.7) +
    scale_fill_manual(values = group_colors) +
    labs(title = paste0(display_name, " protein"), x = NULL, y = "Bulk protein abundance") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

# ---------------------------- Figure 7d-f panels -----------------------------
message_header("Generating Figure 7d-f panels")
for (i in seq_len(nrow(fig7_detail_pairs))) {
  row <- fig7_detail_pairs[i, ]
  plots <- list()
  if (row$gene %in% rownames(sce)) {
    plots[[paste0(row$display_name, " RNA")]] <- make_boxplot_from_seurat(
      sce, row$gene, "RNA", paste0(row$display_name, " RNA"), ylab = "RNA expression"
    )
  }
  if (row$aptamer_feature %in% rownames(sce[["APT"]])) {
    plots[[paste0(row$aptamer_name, " aptamer")]] <- make_boxplot_from_seurat(
      sce, row$aptamer_feature, "APT", paste0(row$aptamer_name, " aptamer"), ylab = "CLR-transformed aptamer signal"
    )
  }
  prot_plot <- make_protein_boxplot(protein_df, row$display_name)
  if (!is.null(prot_plot)) plots[[paste0(row$display_name, " protein")]] <- prot_plot
  if (length(plots) > 0) {
    panel <- wrap_plots(plots, nrow = 1) + plot_annotation(title = row$panel)
    save_plot(panel, row$panel, width = max(7, 3.6 * length(plots)), height = 4.5)
  }
}

# Combined Figure 7d-f if possible.
combined_detail <- list()
for (i in seq_len(nrow(fig7_detail_pairs))) {
  row <- fig7_detail_pairs[i, ]
  p_rna <- if (row$gene %in% rownames(sce)) make_boxplot_from_seurat(sce, row$gene, "RNA", paste0(row$display_name, " RNA"), ylab = "RNA") else NULL
  p_apt <- if (row$aptamer_feature %in% rownames(sce[["APT"]])) make_boxplot_from_seurat(sce, row$aptamer_feature, "APT", paste0(row$aptamer_name, " aptamer"), ylab = "APT") else NULL
  p_prot <- make_protein_boxplot(protein_df, row$display_name)
  combined_detail <- c(combined_detail, Filter(Negate(is.null), list(p_rna, p_apt, p_prot)))
}
if (length(combined_detail) > 0) {
  combined_7d_f <- wrap_plots(combined_detail, ncol = ifelse(is.null(protein_df), 2, 3)) +
    plot_annotation(title = "Figure 7d-f: RNA, aptamer, and bulk-protein comparisons")
  save_plot(combined_7d_f, "Figure7d_f_RNA_aptamer_protein_boxplots_combined", width = ifelse(is.null(protein_df), 9, 12), height = 10)
}

# ------------------------------ save objects ---------------------------------
message_header("Saving analysis objects and run manifest")
saveRDS(sce, file.path(RDS_DIR, "Figure7_drug_resistance_scAptseq_Seurat_object.rds"))
write.csv(sce@meta.data, file.path(TAB_DIR, "cell_metadata_after_QC_and_APT_merge.csv"))

manifest <- data.frame(
  item = c("RNA_DIR", "APT_DIR", "PROTEIN_TABLE", "OUT_DIR", "SAVE_PNG", "RUN_DOUBLET", "MIN_FEATURE_RNA", "MAX_FEATURE_RNA", "MAX_PERCENT_MT", "N_PCS", "RNA_CLUSTER_RES", "APT_CLUSTER_RES", "APT_PCS", "CORRELATION_GROUP_BY", "n_cells_final"),
  value = c(RNA_DIR, APT_DIR, PROTEIN_TABLE, OUT_DIR, SAVE_PNG, RUN_DOUBLET, MIN_FEATURE_RNA, MAX_FEATURE_RNA, MAX_PERCENT_MT, N_PCS, RNA_CLUSTER_RES, APT_CLUSTER_RES, APT_PCS, CORRELATION_GROUP_BY, ncol(sce)),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(TAB_DIR, "run_manifest.csv"), row.names = FALSE)

sink(file.path(OUT_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Analysis complete. Outputs written to: ", OUT_DIR)
