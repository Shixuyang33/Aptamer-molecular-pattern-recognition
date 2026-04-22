#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Figure 5 and supplementary ascites scApt-seq analysis
# Manuscript: Molecular Pattern Recognition by Sequencing Multiplex Aptamers for
# Cell Phenotyping at Single-cell Resolution
#
# Purpose
#   Reproduce the ascites-sample analyses corresponding to Figure 5 and the
#   related supplementary figures in the revised manuscript:
#     - Figure 5a: transcriptome UMAP by annotated cell type and patient/sample
#     - Figure 5b: patient-level mean aptamer-signal heatmap
#     - Figure 5c: CD49C/ZAJ2c mRNA-aptamer UMAP projection and correlation
#     - Figure 5d: CD318/XQ3b mRNA-aptamer UMAP projection and correlation
#     - Figure 5e: ZAJ2c aptamer signal by patient and cell cluster
#     - Figure 5f: XQ3b aptamer signal by patient and cell cluster
#     - Figure S7a: single-cell QC metrics
#     - Figure S7b: canonical marker-gene dot plot
#     - Figure S7c-d: aptamer-based clustering and patient/sample identity
#     - Figure S8: cluster-level mRNA-aptamer correlations for all target pairs
#
# Input directory structure
#   AS_RNA_DIR should contain one folder per RNA library, each with a 10X-style
#   matrix folder, for example:
#       AS_RNA_DIR/AS-1-1/filter_matrix/
#       AS_RNA_DIR/AS-1-2/filter_matrix/
#       AS_RNA_DIR/AS-3-1/filter_matrix/
#       ...
#
#   AS_APT_DIR should contain one folder per aptamer-count library, each with
#   ALL_APT_cell.csv and/or Count_mean.txt, for example:
#       AS_APT_DIR/AS-1-1_9/ALL_APT_cell.csv
#       AS_APT_DIR/AS-1-2_9/ALL_APT_cell.csv
#       ...
#
# Example
#   AS_RNA_DIR=/path/to/10Aptamer_AS \
#   AS_APT_DIR=/path/to/auto_analysis_new \
#   AS_OUT_DIR=Figure5_ascites_outputs \
#   AS_SAVE_PNG=true \
#   Rscript Figure5_ascites_scAptseq_Fig5_Supplementary_GitHub.R
# -----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)
set.seed(20260421)

# ---------------------------- package management -----------------------------
required_pkgs <- c(
  "Seurat", "harmony", "dplyr", "ggplot2", "cowplot", "data.table",
  "pheatmap", "RColorBrewer", "Matrix"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(Matrix)
})

# ------------------------------- configuration -------------------------------
RNA_DIR <- Sys.getenv("AS_RNA_DIR", unset = "./10Aptamer_AS")
APT_DIR <- Sys.getenv("AS_APT_DIR", unset = "./auto_analysis_new")
OUT_DIR <- Sys.getenv("AS_OUT_DIR", unset = "./Figure5_ascites_outputs")
SAVE_PNG <- tolower(Sys.getenv("AS_SAVE_PNG", unset = "false")) %in% c("true", "t", "1", "yes", "y")
RUN_DOUBLET_FINDER <- tolower(Sys.getenv("AS_RUN_DOUBLET_FINDER", unset = "false")) %in% c("true", "t", "1", "yes", "y")

MIN_FEATURE_RNA <- as.numeric(Sys.getenv("AS_MIN_FEATURE_RNA", unset = "400"))
MAX_FEATURE_RNA <- as.numeric(Sys.getenv("AS_MAX_FEATURE_RNA", unset = "4500"))
MAX_PERCENT_MT <- as.numeric(Sys.getenv("AS_MAX_PERCENT_MT", unset = "10"))
RNA_CLUSTER_RES <- as.numeric(Sys.getenv("AS_RNA_CLUSTER_RES", unset = "0.7"))
APT_CLUSTER_RES <- as.numeric(Sys.getenv("AS_APT_CLUSTER_RES", unset = "0.6"))
CORRELATION_GROUP_BY <- Sys.getenv("AS_CORRELATION_GROUP_BY", unset = "rna_seurat_clusters")

# P1 and P2 definitions used in the revised manuscript.
# P1 = CRC ascites sample; P2 = GC ascites sample.
sample_to_patient <- c(
  "AS-1-1" = "P1", "AS-1-2" = "P1",
  "AS-3-1" = "P2", "AS-3-2" = "P2", "AS-3-3" = "P2"
)
patient_labels <- c("P1" = "P1 (CRC)", "P2" = "P2 (GC)")

# Transcript-aptamer target map. The manuscript uses the protein/common marker
# names CD49C and CD318 in Figure 5, while the corresponding gene symbols are
# ITGA3 and CDCP1, respectively.
target_pairs <- data.frame(
  display_name = c("PTK7", "MET", "CD71", "PD-L1", "CD49C", "PTPRF", "CD318", "EpCAM", "ZNRF3"),
  gene = c("PTK7", "MET", "TFRC", "CD274", "ITGA3", "PTPRF", "CDCP1", "EPCAM", "ZNRF3"),
  aptamer_feature = c("PTK7-APT", "MET-APT", "CD71-APT", "PD-L1-APT", "CD49C-APT", "PTPRF-APT", "CD318-APT", "EPCAM-APT", "ZNRF3-APT"),
  aptamer_name = c("sgc8c", "CLN0003", "HG19", "PD-L1", "ZAJ2c", "XQ20a", "XQ3b", "SYL3C", "ZNRF3"),
  stringsAsFactors = FALSE
)
aptamer_order <- target_pairs$aptamer_feature

markergenes <- c("LTB","TCF7","CD3D",#Na_ve CD4+ T
                 "CD8A","GZMA","GZMB","NKG7",#NKT
                 "STMN1","MKI67",#Cycling T
                 "CD14","S100A12","S100A8","LYZ", #CD14+ Monocyte
                 "FCGR3A","MS4A7",#CD16+ Monocyte
                 "CDK1","TK1",#cycling mythoid cell
                 "CLEC10A","CD1C",#Dendritic cell
                 "CD79B","MS4A1","MZB1",#B
                 "IGHM","IGKC","AFF3","PAX5",#Plasma
                 "LAMC2","S100A2","SPARC","NNMT","COL3A1", #fibroblast
                 "EPCAM","TACSTD2","FABP1","TM4SF1" #tumor cell
                 )

# ------------------------------- output paths --------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
fig5_dir <- file.path(OUT_DIR, "Figure5")
supp_dir <- file.path(OUT_DIR, "Supplementary")
tab_dir <- file.path(OUT_DIR, "tables")
rds_dir <- file.path(OUT_DIR, "rds")
for (d in c(fig5_dir, supp_dir, tab_dir, rds_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ------------------------------- utilities -----------------------------------
message2 <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)

standardize_sample_id <- function(x) {
  # Converts folder names such as AS-1-1_9 to AS-1-1.
  x <- basename(as.character(x))
  sub("_.*$", "", x)
}

standardize_cell_barcode <- function(x) {
  # Removes run/library prefixes and formatting differences to allow matching
  # between RNA barcodes and aptamer barcode tables.
  x <- as.character(x)
  x <- basename(x)
  x <- sub("-1$", "", x)
  x <- sub("^[0-9]+[_-]", "", x)
  x <- sub("^[0-9]+", "", x)
  gsub("[^A-Za-z0-9]", "", x)
}

make_cell_key <- function(sample_id, cell_barcode) {
  paste(standardize_sample_id(sample_id), standardize_cell_barcode(cell_barcode), sep = "|")
}

rna_cell_key_from_colname <- function(cell_name) {
  # RNA cells are renamed as sample_barcode by RenameCells(add.cell.id = sample).
  sample_hit <- vapply(names(sample_to_patient), function(s) startsWith(cell_name, paste0(s, "_")), logical(1))
  if (any(sample_hit)) {
    sample_id <- names(sample_to_patient)[which(sample_hit)[1]]
    bc <- sub(paste0("^", sample_id, "_"), "", cell_name)
  } else {
    parts <- strsplit(cell_name, "_", fixed = TRUE)[[1]]
    sample_id <- standardize_sample_id(parts[1])
    bc <- paste(parts[-1], collapse = "_")
  }
  make_cell_key(sample_id, bc)
}

standardize_aptamer_feature <- function(x) {
  y <- toupper(gsub("[ .]", "_", as.character(x)))
  y <- gsub("-", "_", y)
  dplyr::recode(
    y,
    "PTK7" = "PTK7-APT", "PTK7_APT" = "PTK7-APT", "SGC8C" = "PTK7-APT", "SGC8C_APT" = "PTK7-APT",
    "MET" = "MET-APT", "MET_APT" = "MET-APT", "CLN0003" = "MET-APT", "CLN0003_APT" = "MET-APT",
    "CD71" = "CD71-APT", "CD71_APT" = "CD71-APT", "TFRC" = "CD71-APT", "TFRC_APT" = "CD71-APT", "HG19" = "CD71-APT", "HG19_APT" = "CD71-APT",
    "PD_L1" = "PD-L1-APT", "PD_L1_APT" = "PD-L1-APT", "PDL1" = "PD-L1-APT", "PDL1_APT" = "PD-L1-APT", "CD274" = "PD-L1-APT", "CD274_APT" = "PD-L1-APT",
    "CD49C" = "CD49C-APT", "CD49C_APT" = "CD49C-APT", "ITGA3" = "CD49C-APT", "ITGA3_APT" = "CD49C-APT", "ZAJ2C" = "CD49C-APT", "ZAJ2C_APT" = "CD49C-APT",
    "PTPRF" = "PTPRF-APT", "PTPRF_APT" = "PTPRF-APT", "XQ20A" = "PTPRF-APT", "XQ20A_APT" = "PTPRF-APT",
    "CD318" = "CD318-APT", "CD318_APT" = "CD318-APT", "CDCP1" = "CD318-APT", "CDCP1_APT" = "CD318-APT", "XQ3B" = "CD318-APT", "XQ3B_APT" = "CD318-APT",
    "EPCAM" = "EPCAM-APT", "EPCAM_APT" = "EPCAM-APT", "SYL3C" = "EPCAM-APT", "SYL3C_APT" = "EPCAM-APT",
    "ZNRF3" = "ZNRF3-APT", "ZNRF3_APT" = "ZNRF3-APT",
    "SPLIT_APT" = "Scramble-APT", "SCRAMBLE" = "Scramble-APT", "SCRAMBLE_APT" = "Scramble-APT", "CTRL" = "Scramble-APT", "CTRL_APT" = "Scramble-APT",
    .default = as.character(x)
  )
}

get_assay_matrix <- function(object, assay, slot = "data") {
  tryCatch(
    GetAssayData(object, assay = assay, slot = slot),
    error = function(e) GetAssayData(object, assay = assay, layer = slot)
  )
}

save_plot <- function(plot, filename, width, height, dpi = 300) {
  pdf_file <- file.path(dirname(filename), paste0(tools::file_path_sans_ext(basename(filename)), ".pdf"))
  ggsave(pdf_file, plot = plot, width = width, height = height, useDingbats = FALSE)
  if (SAVE_PNG) {
    png_file <- file.path(dirname(filename), paste0(tools::file_path_sans_ext(basename(filename)), ".png"))
    ggsave(png_file, plot = plot, width = width, height = height, dpi = dpi)
  }
  invisible(pdf_file)
}

jacs_theme <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )
}

patient_cols <- c("P1" = "#2C6B73", "P2" = "#8BBE9F")
celltype_cols <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(celltype_levels)), celltype_levels)
apt_cols <- setNames(colorRampPalette(RColorBrewer::brewer.pal(9, "Set3"))(length(aptamer_order)), aptamer_order)

# ----------------------------- load RNA matrices -----------------------------
message2("RNA input directory: ", normalizePath(RNA_DIR, mustWork = FALSE))
message2("APT input directory: ", normalizePath(APT_DIR, mustWork = FALSE))
message2("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE))

rna_sample_dirs <- list.dirs(RNA_DIR, recursive = FALSE, full.names = TRUE)
rna_sample_dirs <- rna_sample_dirs[grepl("^AS", basename(rna_sample_dirs))]
if (length(rna_sample_dirs) == 0) stop("No AS* sample directories were found in AS_RNA_DIR: ", RNA_DIR)

# Prefer manuscript samples if present.
rna_sample_ids <- standardize_sample_id(basename(rna_sample_dirs))
keep <- rna_sample_ids %in% names(sample_to_patient)
rna_sample_dirs <- rna_sample_dirs[keep]
rna_sample_ids <- rna_sample_ids[keep]
if (length(rna_sample_dirs) == 0) stop("RNA sample folders do not match expected sample IDs: ", paste(names(sample_to_patient), collapse = ", "))

# Locate count matrix folders.
locate_10x_dir <- function(sample_dir) {
  candidates <- c(
    file.path(sample_dir, "filter_matrix"),
    file.path(sample_dir, "filtered_feature_bc_matrix"),
    file.path(sample_dir, "filtered_matrix"),
    sample_dir
  )
  candidates[dir.exists(candidates)][1]
}
rna_count_dirs <- vapply(rna_sample_dirs, locate_10x_dir, character(1))
names(rna_count_dirs) <- rna_sample_ids

message2("Loading RNA matrices for samples: ", paste(names(rna_count_dirs), collapse = ", "))
rna_list <- list()
for (sample_id in names(rna_count_dirs)) {
  counts <- Read10X(data.dir = rna_count_dirs[[sample_id]], gene.column = 1)
  obj <- CreateSeuratObject(counts = counts, min.cells = 1, project = sample_id)
  obj <- RenameCells(obj, add.cell.id = sample_id)
  obj$orig.ident <- sample_id
  obj$patient <- sample_to_patient[sample_id]
  obj$patient_label <- patient_labels[obj$patient]
  rna_list[[sample_id]] <- obj
}

combined <- Reduce(function(x, y) merge(x, y), rna_list)
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
combined$cell_key <- vapply(colnames(combined), rna_cell_key_from_colname, character(1))

# RNA QC filtering.
combined <- subset(
  combined,
  subset = nFeature_RNA > MIN_FEATURE_RNA &
    nFeature_RNA < MAX_FEATURE_RNA &
    percent.mt < MAX_PERCENT_MT
)

# Optional doublet removal. Disabled by default for GitHub reproducibility because
# DoubletFinder function names differ across Seurat versions.
if (RUN_DOUBLET_FINDER) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    warning("AS_RUN_DOUBLET_FINDER=true, but DoubletFinder is not installed. Skipping doublet detection.")
  } else {
    message2("DoubletFinder is installed, but automatic doublet removal is not run in this GitHub script because DoubletFinder APIs vary by version. Please run it separately if needed.")
  }
}

# RNA normalization, integration, clustering, and UMAP.
message2("Running RNA normalization, PCA, Harmony integration, clustering, and UMAP.")
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
combined <- RunHarmony(combined, group.by.vars = "orig.ident", reduction = "pca", assay.use = "RNA", dims.use = 1:30, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = RNA_CLUSTER_RES, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, reduction.name = "umap", verbose = FALSE)
combined$rna_seurat_clusters <- as.character(combined$seurat_clusters)

# Apply  cell-type annotation.

pbmc.markers <- FindAllMarkers(combined, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
new.cluster.ids <- c("Tumor cell-1","Naive CD4T","CD14+ Monocyte","CD16+ Monocyte","Effect CD8T","Tumor cell-2","Naive CD8T","Fibrolast-1","Fibroblast-2","DC-1","DC-2","Fibroblast-3","B cell","Fibroblast-4",
                     "Cycling T","NK","Plasma","pDC","Fubroblast-5","Cycling myeloid cell")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
DimPlot(combined, reduction = "harmony", label = TRUE, pt.size = 0.5) + NoLegend()
combined$Celltypes_new <- combined1@active.ident
head(combined@meta.data)
table(combined$Celltypes_new)
my_levels <- c("Naive CD4T","Naive CD8T","Effect CD8T","Cycling T","NK","CD14+ Monocyte","CD16+ Monocyte","Cycling myeloid cell","DC-1","DC-2","pDC","B cell","Plasma","Fibrolast-1","Fibroblast-2","Fibroblast-3","Fibroblast-4","Fubroblast-5","Tumor cell-1","Tumor cell-2")
combined$Celltypes_new <- factor(combined$Celltypes_new,levels = my_levels)


# ------------------------------ load APT counts ------------------------------
message2("Loading aptamer count tables.")
apt_sample_dirs <- list.dirs(APT_DIR, recursive = FALSE, full.names = TRUE)
apt_sample_dirs <- apt_sample_dirs[grepl("^AS", basename(apt_sample_dirs))]
if (length(apt_sample_dirs) == 0) stop("No AS* sample directories were found in AS_APT_DIR: ", APT_DIR)

read_apt_wide <- function(sample_dir) {
  raw_sample <- basename(sample_dir)
  sample_id <- standardize_sample_id(raw_sample)
  cell_file <- file.path(sample_dir, "ALL_APT_cell.csv")
  if (!file.exists(cell_file)) {
    stop("Missing ALL_APT_cell.csv in ", sample_dir)
  }
  df <- as.data.frame(data.table::fread(cell_file))
  cell_col <- intersect(c("cell", "Cell", "barcode", "Barcode", "cell_barcode"), names(df))[1]
  if (is.na(cell_col)) cell_col <- names(df)[1]
  df$sample_id <- sample_id
  df$cell_raw <- df[[cell_col]]
  feature_cols <- setdiff(names(df), c(cell_col, "sample_id", "cell_raw", "sample", "cellnumber", "apt", "sum"))
  # Keep numeric feature columns only.
  numeric_cols <- feature_cols[vapply(df[feature_cols], function(z) suppressWarnings(all(!is.na(as.numeric(z)))), logical(1))]
  df[numeric_cols] <- lapply(df[numeric_cols], function(z) as.numeric(z))
  feature_std <- standardize_aptamer_feature(numeric_cols)
  names(df)[match(numeric_cols, names(df))] <- feature_std
  # If different columns map to the same standardized feature, sum them.
  df$cell_key <- make_cell_key(sample_id, df$cell_raw)
  feature_std_unique <- unique(feature_std)
  out <- df[, c("cell_key", "sample_id", feature_std_unique), drop = FALSE]
  out <- out %>%
    group_by(cell_key, sample_id) %>%
    summarise(across(all_of(feature_std_unique), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  out
}

apt_wide <- lapply(apt_sample_dirs, read_apt_wide) %>% bind_rows()
apt_wide <- apt_wide[apt_wide$sample_id %in% names(sample_to_patient), , drop = FALSE]
if (nrow(apt_wide) == 0) stop("No aptamer count rows matched expected sample IDs.")

apt_features_found <- intersect(c(aptamer_order, "Scramble-APT"), names(apt_wide))
if (length(apt_features_found) < 2) stop("Too few aptamer features were detected. Check ALL_APT_cell.csv column names.")

rna_key_to_cell <- data.frame(cell_name = colnames(combined), cell_key = combined$cell_key, stringsAsFactors = FALSE)
apt_mat_df <- apt_wide[, c("cell_key", apt_features_found), drop = FALSE]
apt_mat_df <- apt_mat_df %>% group_by(cell_key) %>% summarise(across(all_of(apt_features_found), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
merged_keys <- inner_join(rna_key_to_cell, apt_mat_df, by = "cell_key")
message2("Matched ", nrow(merged_keys), " RNA cells with aptamer-count cells out of ", ncol(combined), " RNA cells after QC.")
if (nrow(merged_keys) < 100) warning("Fewer than 100 matched cells were found. Please check RNA/APT barcode conventions.")

combined <- combined[, merged_keys$cell_name]
apt_counts <- as.matrix(t(merged_keys[, apt_features_found, drop = FALSE]))
colnames(apt_counts) <- merged_keys$cell_name
rownames(apt_counts) <- apt_features_found
apt_counts <- Matrix::Matrix(apt_counts, sparse = TRUE)
combined[["APT"]] <- CreateAssayObject(counts = apt_counts)
combined$nCount_APT <- Matrix::colSums(get_assay_matrix(combined, "APT", "counts"))
combined$nFeature_APT <- Matrix::colSums(get_assay_matrix(combined, "APT", "counts") > 0)

DefaultAssay(combined) <- "APT"
combined <- NormalizeData(combined, normalization.method = "CLR", margin = 2, verbose = FALSE)
DefaultAssay(combined) <- "RNA"

# Save QC metrics for Table S5 consistency checks.
qc_summary <- data.frame(
  dataset = "Ascites Fluid Samples",
  cell_number = ncol(combined),
  median_gene_number = median(combined$nFeature_RNA),
  median_umi = median(combined$nCount_RNA),
  median_mitochondrial_percentage = median(combined$percent.mt),
  median_aptamer_counts = median(combined$nCount_APT),
  stringsAsFactors = FALSE
)
write.csv(qc_summary, file.path(tab_dir, "ascites_QC_summary_for_TableS5.csv"), row.names = FALSE)

# ------------------------- aptamer-based clustering --------------------------
message2("Running aptamer-based clustering for Figure S7c-d.")
apt_features_for_clustering <- intersect(aptamer_order, rownames(combined[["APT"]]))
if (length(apt_features_for_clustering) >= 3) {
  DefaultAssay(combined) <- "APT"
  combined <- ScaleData(combined, features = apt_features_for_clustering, verbose = FALSE)
  npcs_apt <- max(2, min(8, length(apt_features_for_clustering) - 1))
  combined <- RunPCA(
    combined,
    features = apt_features_for_clustering,
    npcs = npcs_apt,
    reduction.name = "apt_pca",
    reduction.key = "APTPC_",
    verbose = FALSE
  )
  apt_dims <- 1:npcs_apt
  combined <- FindNeighbors(combined, reduction = "apt_pca", dims = apt_dims, graph.name = "APT_snn", verbose = FALSE)
  combined <- FindClusters(combined, graph.name = "APT_snn", resolution = APT_CLUSTER_RES, verbose = FALSE)
  combined$apt_clusters <- as.character(Idents(combined))
  combined <- RunUMAP(
    combined,
    reduction = "apt_pca",
    dims = apt_dims,
    reduction.name = "apt_umap",
    reduction.key = "APTUMAP_",
    verbose = FALSE
  )
  Idents(combined) <- combined$rna_seurat_clusters
  DefaultAssay(combined) <- "RNA"
} else {
  warning("Not enough aptamer features were found for aptamer-based clustering.")
  combined$apt_clusters <- NA_character_
}

# --------------------------- Figure S7a: QC plots ----------------------------
qc_df <- combined@meta.data
qc_df$orig.ident <- factor(qc_df$orig.ident, levels = names(sample_to_patient))
qc_long <- bind_rows(
  data.frame(sample = qc_df$orig.ident, metric = "APT counts", value = qc_df$nCount_APT),
  data.frame(sample = qc_df$orig.ident, metric = "Gene counts", value = qc_df$nFeature_RNA),
  data.frame(sample = qc_df$orig.ident, metric = "UMI counts", value = qc_df$nCount_RNA),
  data.frame(sample = qc_df$orig.ident, metric = "Percent MT gene", value = qc_df$percent.mt)
)
qc_plot <- ggplot(qc_long, aes(x = sample, y = value, fill = sample)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.2) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(levels(qc_df$orig.ident)))) +
  labs(x = "Sample", y = "Counts / percentage") +
  jacs_theme(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_plot(qc_plot, file.path(supp_dir, "FigureS7a_ascites_single_cell_QC_violin.pdf"), width = 10, height = 3.2)

# ------------------------ Figure 5a and Figure S7b ---------------------------
fig5a_celltype <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "Celltypes_new",
  cols = celltype_cols,
  raster = FALSE,
  label = FALSE
) +
  labs(title = "Transcriptome-based clusters") +
  jacs_theme(base_size = 9) +
  theme(legend.position = "right")
fig5a_patient <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "patient",
  cols = patient_cols,
  raster = FALSE,
  shuffle = TRUE
) +
  labs(title = "Sample identity") +
  jacs_theme(base_size = 9) +
  theme(legend.position = "right")
fig5a <- plot_grid(fig5a_celltype, fig5a_patient, nrow = 1, rel_widths = c(1.25, 1))
save_plot(fig5a, file.path(fig5_dir, "Figure5a_transcriptome_UMAP_celltype_and_patient.pdf"), width = 10, height = 5)

marker_genes_present <- intersect(marker_genes, rownames(combined[["RNA"]]))
if (length(marker_genes_present) > 0) {
  marker_dot <- DotPlot(
    combined,
    features = marker_genes_present,
    group.by = "Celltypes_new",
    assay = "RNA",
    cols = c("lightgrey", "red")
  ) +
    coord_flip() +
    labs(title = "Canonical marker genes", x = NULL, y = NULL) +
    jacs_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(marker_dot, file.path(supp_dir, "FigureS7b_canonical_marker_gene_dotplot.pdf"), width = 9, height = 10)
}

# ------------------------ Figure S7c-d: APT UMAP -----------------------------
if ("apt_umap" %in% names(combined@reductions)) {
  figS7c <- DimPlot(
    combined,
    reduction = "apt_umap",
    group.by = "apt_clusters",
    label = TRUE,
    raster = FALSE
  ) + labs(title = "Aptamer-based clustering") + jacs_theme(base_size = 9)
  figS7d_patient <- DimPlot(
    combined,
    reduction = "apt_umap",
    group.by = "patient",
    cols = patient_cols,
    raster = FALSE,
    shuffle = TRUE
  ) + labs(title = "Patient identity") + jacs_theme(base_size = 9)
  figS7d_sample <- DimPlot(
    combined,
    reduction = "apt_umap",
    group.by = "orig.ident",
    raster = FALSE,
    shuffle = TRUE
  ) + labs(title = "Library identity") + jacs_theme(base_size = 9)
  save_plot(figS7c, file.path(supp_dir, "FigureS7c_aptamer_based_UMAP_clusters.pdf"), width = 5.5, height = 4.8)
  save_plot(plot_grid(figS7d_patient, figS7d_sample, nrow = 1), file.path(supp_dir, "FigureS7d_aptamer_based_UMAP_patient_and_sample.pdf"), width = 10, height = 4.8)
}

# -------------------------- Figure 5b: heatmap -------------------------------
DefaultAssay(combined) <- "APT"
apt_data <- get_assay_matrix(combined, "APT", "data")
apt_features_present <- intersect(aptamer_order, rownames(apt_data))
patient_means <- sapply(c("P1", "P2"), function(pat) {
  cells <- colnames(combined)[combined$patient == pat]
  Matrix::rowMeans(apt_data[apt_features_present, cells, drop = FALSE])
})
rownames(patient_means) <- apt_features_present
patient_means_scaled <- t(scale(t(as.matrix(patient_means))))
patient_means_scaled[is.na(patient_means_scaled)] <- 0

pdf(file.path(fig5_dir, "Figure5b_patient_mean_aptamer_signal_heatmap.pdf"), width = 4.5, height = 5.5, useDingbats = FALSE)
pheatmap(
  patient_means_scaled,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey70",
  main = "Mean aptamer signal",
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 0
)
dev.off()
if (SAVE_PNG) {
  png(file.path(fig5_dir, "Figure5b_patient_mean_aptamer_signal_heatmap.png"), width = 1500, height = 1800, res = 300)
  pheatmap(
    patient_means_scaled,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "grey70",
    main = "Mean aptamer signal",
    fontsize_row = 9,
    fontsize_col = 10,
    angle_col = 0
  )
  dev.off()
}
write.csv(patient_means, file.path(tab_dir, "Figure5b_patient_mean_CLR_aptamer_signal.csv"))

# Optional sample-level heatmap to document the libraries underlying P1/P2.
sample_means <- sapply(names(sample_to_patient), function(samp) {
  cells <- colnames(combined)[combined$orig.ident == samp]
  if (length(cells) == 0) return(rep(NA_real_, length(apt_features_present)))
  Matrix::rowMeans(apt_data[apt_features_present, cells, drop = FALSE])
})
rownames(sample_means) <- apt_features_present
sample_means_scaled <- t(scale(t(as.matrix(sample_means))))
sample_means_scaled[is.na(sample_means_scaled)] <- 0
pdf(file.path(supp_dir, "FigureS7e_sample_level_mean_aptamer_signal_heatmap.pdf"), width = 6, height = 5.5, useDingbats = FALSE)
pheatmap(
  sample_means_scaled,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey70",
  main = "Mean aptamer signal by library",
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45
)
dev.off()
write.csv(sample_means, file.path(tab_dir, "FigureS7e_sample_mean_CLR_aptamer_signal.csv"))

# ------------------------- correlation utilities -----------------------------
cluster_average <- function(object, assay, feature, group_by, slot = "data") {
  mat <- get_assay_matrix(object, assay = assay, slot = slot)
  if (!feature %in% rownames(mat)) stop("Feature not found in ", assay, ": ", feature)
  groups <- object@meta.data[[group_by]]
  if (is.null(groups)) stop("group_by column not found in metadata: ", group_by)
  tapply(as.numeric(mat[feature, colnames(object)]), groups, mean, na.rm = TRUE)
}

correlation_df_for_pair <- function(object, gene, aptamer_feature, display_name, aptamer_name, group_by = CORRELATION_GROUP_BY) {
  rna_avg <- cluster_average(object, "RNA", gene, group_by = group_by, slot = "data")
  apt_avg <- cluster_average(object, "APT", aptamer_feature, group_by = group_by, slot = "data")
  common <- intersect(names(rna_avg), names(apt_avg))
  df <- data.frame(
    group = common,
    RNA = as.numeric(rna_avg[common]),
    APT = as.numeric(apt_avg[common]),
    display_name = display_name,
    gene = gene,
    aptamer_feature = aptamer_feature,
    aptamer_name = aptamer_name,
    stringsAsFactors = FALSE
  )
  if (nrow(df) >= 3 && stats::sd(df$RNA) > 0 && stats::sd(df$APT) > 0) {
    ct <- suppressWarnings(cor.test(df$RNA, df$APT, method = "pearson"))
    df$r <- as.numeric(ct$estimate)
    df$p_value <- ct$p.value
  } else {
    df$r <- NA_real_
    df$p_value <- NA_real_
  }
  df
}

plot_correlation <- function(df) {
  r_lab <- ifelse(is.na(df$r[1]), "r = NA", paste0("r = ", sprintf("%.3f", df$r[1])))
  ggplot(df, aes(x = RNA, y = APT)) +
    geom_point(size = 2.2, color = "black", alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70", linewidth = 0.4) +
    geom_text(aes(label = group), vjust = -0.7, size = 2.5) +
    labs(
      title = paste0(unique(df$display_name), ": mRNA vs aptamer"),
      subtitle = r_lab,
      x = paste0(unique(df$gene), " RNA (cluster mean)"),
      y = paste0(unique(df$aptamer_name), " / ", unique(df$aptamer_feature), " (cluster mean)")
    ) +
    jacs_theme(base_size = 9)
}

all_correlation_dfs <- lapply(seq_len(nrow(target_pairs)), function(i) {
  pair <- target_pairs[i, ]
  if (!pair$gene %in% rownames(combined[["RNA"]]) || !pair$aptamer_feature %in% rownames(combined[["APT"]])) {
    warning("Skipping missing target pair: ", pair$gene, " / ", pair$aptamer_feature)
    return(NULL)
  }
  correlation_df_for_pair(
    combined,
    gene = pair$gene,
    aptamer_feature = pair$aptamer_feature,
    display_name = pair$display_name,
    aptamer_name = pair$aptamer_name,
    group_by = CORRELATION_GROUP_BY
  )
})
all_correlation_dfs <- Filter(Negate(is.null), all_correlation_dfs)
correlation_long <- bind_rows(all_correlation_dfs)
correlation_summary <- correlation_long %>%
  distinct(display_name, gene, aptamer_name, aptamer_feature, r, p_value) %>%
  arrange(match(aptamer_feature, aptamer_order))
write.csv(correlation_summary, file.path(tab_dir, "FigureS8_cluster_level_mRNA_aptamer_correlation_summary.csv"), row.names = FALSE)
write.csv(correlation_long, file.path(tab_dir, "FigureS8_cluster_level_mRNA_aptamer_correlation_long.csv"), row.names = FALSE)

cor_plots <- lapply(split(correlation_long, correlation_long$aptamer_feature), plot_correlation)
# Order correlation plots by manuscript target order.
cor_plots <- cor_plots[intersect(aptamer_order, names(cor_plots))]
figS8 <- cowplot::plot_grid(plotlist = cor_plots, ncol = 3)
save_plot(figS8, file.path(supp_dir, "FigureS8_all_cluster_level_mRNA_aptamer_correlations.pdf"), width = 12, height = 10)
for (nm in names(cor_plots)) {
  safe_nm <- gsub("[^A-Za-z0-9]+", "_", nm)
  save_plot(cor_plots[[nm]], file.path(supp_dir, paste0("FigureS8_", safe_nm, "_correlation.pdf")), width = 4, height = 3.5)
}

# ------------------------- Figure 5c-d: key pairs ----------------------------
feature_plot_pair <- function(object, gene, aptamer_feature, display_name, aptamer_name, out_prefix) {
  DefaultAssay(object) <- "RNA"
  p_rna <- FeaturePlot(
    object,
    features = gene,
    reduction = "umap",
    cols = c("lightgrey", "purple"),
    raster = FALSE
  ) +
    labs(title = paste0(display_name, " mRNA (", gene, ")")) +
    jacs_theme(base_size = 9)
  DefaultAssay(object) <- "APT"
  p_apt <- FeaturePlot(
    object,
    features = aptamer_feature,
    reduction = "umap",
    cols = c("lightgrey", "darkgreen"),
    raster = FALSE
  ) +
    labs(title = paste0(aptamer_name, " aptamer (", display_name, ")")) +
    jacs_theme(base_size = 9)
  cdf <- correlation_df_for_pair(object, gene, aptamer_feature, display_name, aptamer_name, group_by = CORRELATION_GROUP_BY)
  p_cor <- plot_correlation(cdf)
  p_combined <- plot_grid(p_rna, p_apt, p_cor, nrow = 1, rel_widths = c(1, 1, 1.05))
  save_plot(p_combined, file.path(fig5_dir, paste0(out_prefix, "_UMAP_and_correlation.pdf")), width = 12, height = 4)
  save_plot(p_rna, file.path(fig5_dir, paste0(out_prefix, "_mRNA_UMAP.pdf")), width = 4, height = 4)
  save_plot(p_apt, file.path(fig5_dir, paste0(out_prefix, "_aptamer_UMAP.pdf")), width = 4, height = 4)
  save_plot(p_cor, file.path(fig5_dir, paste0(out_prefix, "_cluster_level_correlation.pdf")), width = 4.2, height = 3.6)
  invisible(p_combined)
}

# Figure 5c: CD49C mRNA/ZAJ2c aptamer.
feature_plot_pair(
  combined,
  gene = "ITGA3",
  aptamer_feature = "CD49C-APT",
  display_name = "CD49C",
  aptamer_name = "ZAJ2c",
  out_prefix = "Figure5c_CD49C_ZAJ2c"
)

# Figure 5d: CD318 mRNA/XQ3b aptamer.
feature_plot_pair(
  combined,
  gene = "CDCP1",
  aptamer_feature = "CD318-APT",
  display_name = "CD318",
  aptamer_name = "XQ3b",
  out_prefix = "Figure5d_CD318_XQ3b"
)
DefaultAssay(combined) <- "RNA"

# ----------------------- Figure 5e-f: violin plots ---------------------------
fetch_aptamer_values <- function(object, aptamer_feature) {
  mat <- get_assay_matrix(object, "APT", "data")
  if (!aptamer_feature %in% rownames(mat)) stop("APT feature not found: ", aptamer_feature)
  as.numeric(mat[aptamer_feature, colnames(object)])
}

violin_patient_and_celltype <- function(object, aptamer_feature, aptamer_name, display_name, out_prefix) {
  df <- object@meta.data
  df$value <- fetch_aptamer_values(object, aptamer_feature)
  df$patient <- factor(df$patient, levels = c("P1", "P2"))
  df$Celltypes_new <- factor(df$Celltypes_new, levels = levels(object$Celltypes_new))
  p_patient <- ggplot(df, aes(x = patient, y = value, fill = patient)) +
    geom_violin(scale = "width", linewidth = 0.25, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    scale_fill_manual(values = patient_cols) +
    labs(title = paste0(aptamer_name, " by patient"), x = NULL, y = "CLR-transformed aptamer signal") +
    jacs_theme(base_size = 9) +
    theme(legend.position = "none")
  p_cluster <- ggplot(df, aes(x = Celltypes_new, y = value, fill = Celltypes_new)) +
    geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
    scale_fill_manual(values = celltype_cols, na.value = "grey80") +
    labs(title = paste0(aptamer_name, " across cell clusters"), x = NULL, y = "CLR-transformed aptamer signal") +
    jacs_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  p_combined <- plot_grid(p_patient, p_cluster, nrow = 1, rel_widths = c(0.8, 2.2))
  save_plot(p_combined, file.path(fig5_dir, paste0(out_prefix, "_patient_and_celltype_violin.pdf")), width = 11, height = 4)
  save_plot(p_patient, file.path(fig5_dir, paste0(out_prefix, "_patient_violin.pdf")), width = 3.6, height = 3.6)
  save_plot(p_cluster, file.path(fig5_dir, paste0(out_prefix, "_celltype_violin.pdf")), width = 8, height = 3.6)
  invisible(p_combined)
}

violin_patient_and_celltype(
  combined,
  aptamer_feature = "CD49C-APT",
  aptamer_name = "ZAJ2c",
  display_name = "CD49C",
  out_prefix = "Figure5e_ZAJ2c_CD49C_APT"
)
violin_patient_and_celltype(
  combined,
  aptamer_feature = "CD318-APT",
  aptamer_name = "XQ3b",
  display_name = "CD318",
  out_prefix = "Figure5f_XQ3b_CD318_APT"
)

# -------------------------- save final objects -------------------------------
saveRDS(combined, file.path(rds_dir, "ascites_scAptseq_Figure5_final_object.rds"))
write.csv(combined@meta.data, file.path(tab_dir, "ascites_scAptseq_Figure5_metadata.csv"), row.names = TRUE)

capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))

message2("Analysis complete.")
message2("Main Figure 5 outputs: ", normalizePath(fig5_dir, mustWork = FALSE))
message2("Supplementary outputs: ", normalizePath(supp_dir, mustWork = FALSE))
message2("Tables and metadata: ", normalizePath(tab_dir, mustWork = FALSE))
