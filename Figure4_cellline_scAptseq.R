#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Figure 4 and related supplementary figures: cell-line scApt-seq analysis
# Manuscript: Molecular Pattern Recognition by Sequencing Multiplex Aptamers
#             for Cell Phenotyping at Single-cell Resolution
#
# This script reproduces the analyses corresponding to revised Figure 4 and
# the related supplementary analyses for the HCT-8 concentration-test cell-line
# dataset and the PTK7-knockout specificity-control dataset.
#
# Main outputs
#   Figure 4a: QC violin plots for HCT-8 Test1/Test2 groups
#   Figure 4b: HCT-8 transcriptome UMAP, group identity, PTK7 RNA, sgc8c, scramble
#   Figure 4c: Relative composition of the nine aptamers plus scrambled control
#   Figure 4d: PTK7-knockout/negative-control UMAP and feature plots
#   Figure 4e: sgc8c and scrambled aptamer box plots in NC versus PTK7-KO cells
#
# Supplementary outputs
#   Figure S5-like: UMAP feature plots for aptamer-target RNA and aptamer signals
#   Figure S6-like: cluster-level RNA-aptamer correlation plots
#   CSV summaries: QC metrics, aptamer proportions, and correlation coefficients
#
# Usage example
#   HCT8_RNA_DIR=/path/to/10Aptamer_cell_line_test \
#   HCT8_APT_DIR=/path/to/10APT_cellline_result/autoanalysis \
#   KO_RNA_DIR=/path/to/PTK7_KO_RNA \
#   KO_APT_DIR=/path/to/PTK7_KO_APT \
#   FIG4_OUT_DIR=Figure4_cellline_outputs \
#   SAVE_PNG=true \
#   Rscript Figure4_cellline_scAptseq_GitHub.R
#
# Notes
#   1. RNA input folders should contain one sample directory per library, with a
#      filtered matrix folder named either filter_matrix, filtered_feature_bc_matrix,
#      or filtered_gene_bc_matrices.
#   2. Aptamer input folders should contain one sample directory per library and
#      either count/ALL_APT_cell.csv for cell-level aptamer counts or
#      count/Count_mean.txt for sample-level summaries.
#   3. The script is deliberately written without hard-coded absolute paths so it
#      can be uploaded to GitHub and reused on another system.
# -----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# ------------------------------ configuration -------------------------------
get_env <- function(name, default = NULL) {
  value <- Sys.getenv(name)
  if (identical(value, "")) default else value
}

HCT8_RNA_DIR <- get_env("HCT8_RNA_DIR", get_env("CELL_RNA_DIR", "./10Aptamer_cell_line_test"))
HCT8_APT_DIR <- get_env("HCT8_APT_DIR", get_env("CELL_APT_DIR", "./10APT_cellline_result/autoanalysis"))
KO_RNA_DIR   <- get_env("KO_RNA_DIR",   "")
KO_APT_DIR   <- get_env("KO_APT_DIR",   "")
OUT_DIR      <- get_env("FIG4_OUT_DIR", "Figure4_cellline_outputs")
SAVE_PNG     <- tolower(get_env("SAVE_PNG", "false")) %in% c("true", "1", "yes", "y")
RUN_DOUBLET  <- tolower(get_env("RUN_DOUBLET", "false")) %in% c("true", "1", "yes", "y")

MIN_GENES <- as.integer(get_env("MIN_GENES", "400"))
MAX_GENES <- as.integer(get_env("MAX_GENES", "4500"))
MAX_MT    <- as.numeric(get_env("MAX_MT", "10"))
N_PCS     <- as.integer(get_env("N_PCS", "30"))
RESOLUTION <- as.numeric(get_env("CLUSTER_RESOLUTION", "0.5"))

# ------------------------------ package setup --------------------------------
required_packages <- c(
  "Seurat", "dplyr", "tidyr", "ggplot2", "patchwork", "data.table",
  "Matrix", "pheatmap"
)
optional_packages <- c("harmony", "DoubletFinder", "ggsignif")

missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_required, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

HAS_HARMONY <- requireNamespace("harmony", quietly = TRUE)
HAS_DF <- requireNamespace("DoubletFinder", quietly = TRUE)
HAS_GGSIGNIF <- requireNamespace("ggsignif", quietly = TRUE)

if (RUN_DOUBLET && !HAS_DF) {
  warning("RUN_DOUBLET=true, but DoubletFinder is not installed. Doublet filtering will be skipped.")
  RUN_DOUBLET <- FALSE
}

# ------------------------------ output folders -------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "objects"), recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, width, height) {
  pdf_path <- file.path(OUT_DIR, "figures", paste0(name, ".pdf"))
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  if (SAVE_PNG) {
    png_path <- file.path(OUT_DIR, "figures", paste0(name, ".png"))
    ggplot2::ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = 300)
  }
}

# ------------------------------ common metadata ------------------------------
aptamer_pairs <- data.frame(
  aptamer_id = c("Sgc8c", "CLN0003", "HG19", "PD-L1", "ZAJ2c", "XQ20a", "XQ3b", "SYL3C", "ZNRF3"),
  aptamer_feature = c("Sgc8c_Aptamer", "CLN0003_Aptamer", "HG19_Aptamer", "PD-L1_Aptamer", "ZAJ2c_Aptamer", "XQ20a_Aptamer", "XQ3b_Aptamer", "SYL3C_Aptamer", "ZNRF3_Aptamer"),
  target_gene = c("PTK7", "MET", "TFRC", "CD274", "ITGA3", "PTPRF", "CDCP1", "EPCAM", "ZNRF3"),
  target_name = c("PTK7", "MET", "CD71", "PD-L1", "CD49C", "PTPRF", "CD318", "EpCAM", "ZNRF3"),
  stringsAsFactors = FALSE
)

aptamer_order <- c(aptamer_pairs$aptamer_feature, "Scramble_Aptamer")

aptamer_palette <- c(
  Sgc8c_Aptamer = "#8dd3c7",
  CLN0003_Aptamer = "#ffffb3",
  HG19_Aptamer = "#bebada",
  `PD-L1_Aptamer` = "#fccde5",
  ZAJ2c_Aptamer = "#80b1d3",
  XQ20a_Aptamer = "#fdb462",
  XQ3b_Aptamer = "#b3de69",
  SYL3C_Aptamer = "#f4cae4",
  ZNRF3_Aptamer = "#d9d9d9",
  Scramble_Aptamer = "#bc80bd"
)

group_palette_hct8 <- c(
  `Test1-1` = "#1E4F5F",
  `Test1-2` = "#4C9D80",
  `Test2-1` = "#C0E0BF",
  `Test2-2` = "#106232"
)

group_palette_ko <- c(
  WT = "#A6CEE3",
  NC = "#1E4F5F",
  KO = "#4C9D80",
  `PTK7-KO` = "#4C9D80"
)

cluster_palette <- c(
  "#33a02c", "#ff69b4", "#1f78b4", "#ff4500", "#b2df8a", "#a6cee3",
  "#a020f0", "#ff7f00", "#66cdaa", "#db7093", "#1e90ff", "#fdbf6f",
  "#6a3d9a", "#ffd700", "#b15928", "#8dd3c7", "#ffffb3", "#bebada"
)

# ------------------------------ helper functions -----------------------------
standardize_sample_label_hct8 <- function(x) {
  y <- basename(x)
  y <- gsub("^Ten1.*|^Test1.*1$|20.*1", "Test1-1", y, ignore.case = TRUE)
  y <- gsub("^Ten2.*|^Test1.*2$|20.*2", "Test1-2", y, ignore.case = TRUE)
  y <- gsub("^Ten3.*|^Test2.*1$|200.*1", "Test2-1", y, ignore.case = TRUE)
  y <- gsub("^Ten4.*|^Test2.*2$|200.*2", "Test2-2", y, ignore.case = TRUE)
  y
}

standardize_sample_label_ko <- function(x) {
  y <- basename(x)
  if (grepl("(^|[_-])WT($|[_-])|wild", y, ignore.case = TRUE)) return("WT")
  if (grepl("(^|[_-])NC($|[_-])|neg|negative|control", y, ignore.case = TRUE)) return("NC")
  if (grepl("KO|knock", y, ignore.case = TRUE)) return("KO")
  y
}

standardize_aptamer_name <- function(x) {
  y <- trimws(as.character(x))
  y <- gsub("\\s+", "", y)
  y <- gsub("-", "_", y)
  y <- gsub("\\.", "_", y)
  upper <- toupper(y)

  out <- y
  out[grepl("SGC8|PTK7", upper)] <- "Sgc8c_Aptamer"
  out[grepl("CLN0003|^MET(_|$)|METAPT", upper)] <- "CLN0003_Aptamer"
  out[grepl("HG19|CD71|TFRC", upper)] <- "HG19_Aptamer"
  out[grepl("PDL1|PD_L1|CD274", upper)] <- "PD-L1_Aptamer"
  out[grepl("ZAJ2C|CD49C|ITGA3", upper)] <- "ZAJ2c_Aptamer"
  out[grepl("XQ20A|PTPRF", upper)] <- "XQ20a_Aptamer"
  out[grepl("XQ3B|CD318|CDCP1", upper)] <- "XQ3b_Aptamer"
  out[grepl("SYL3C|EPCAM|EPC", upper)] <- "SYL3C_Aptamer"
  out[grepl("ZNRF3", upper)] <- "ZNRF3_Aptamer"
  out[grepl("CONTR|CTRL|CONTROL|SCRAM", upper)] <- "Scramble_Aptamer"
  out
}

find_matrix_dir <- function(sample_dir) {
  candidates <- c(
    file.path(sample_dir, "filter_matrix"),
    file.path(sample_dir, "filtered_feature_bc_matrix"),
    file.path(sample_dir, "filtered_gene_bc_matrices"),
    sample_dir
  )
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

list_sample_dirs <- function(root_dir, label_fun) {
  if (!dir.exists(root_dir)) stop("Input RNA directory does not exist: ", root_dir)
  dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  if (length(dirs) == 0) dirs <- root_dir
  mat_dirs <- vapply(dirs, find_matrix_dir, character(1))
  keep <- !is.na(mat_dirs)
  mat_dirs <- mat_dirs[keep]
  dirs <- dirs[keep]
  labels <- vapply(dirs, label_fun, character(1))
  names(mat_dirs) <- labels
  mat_dirs
}

read_rna_samples <- function(root_dir, label_fun) {
  sample_dirs <- list_sample_dirs(root_dir, label_fun)
  message("RNA samples detected: ", paste(names(sample_dirs), collapse = ", "))

  objs <- lapply(seq_along(sample_dirs), function(i) {
    label <- names(sample_dirs)[i]
    counts <- Seurat::Read10X(data.dir = sample_dirs[i], gene.column = 1)
    obj <- Seurat::CreateSeuratObject(counts = counts, project = label, min.cells = 1)
    obj$sample_group <- label
    obj$sample_id <- label
    colnames(obj) <- paste(label, colnames(obj), sep = "_")
    obj
  })

  if (length(objs) == 1) {
    combined <- objs[[1]]
  } else {
    combined <- merge(objs[[1]], y = objs[-1], merge.data = TRUE)
  }
  combined[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined, pattern = "^MT-")
  combined
}

read_aptamer_cell_counts <- function(apt_root, label_fun) {
  if (!dir.exists(apt_root)) stop("Input aptamer directory does not exist: ", apt_root)
  sample_dirs <- list.dirs(apt_root, recursive = FALSE, full.names = TRUE)
  if (length(sample_dirs) == 0) sample_dirs <- apt_root

  all_counts <- list()
  for (sample_dir in sample_dirs) {
    label <- label_fun(sample_dir)
    file <- file.path(sample_dir, "count", "ALL_APT_cell.csv")
    if (!file.exists(file)) {
      # Some workflows place ALL_APT_cell.csv directly under the sample folder.
      alt <- file.path(sample_dir, "ALL_APT_cell.csv")
      if (file.exists(alt)) file <- alt
    }
    if (!file.exists(file)) next

    df <- data.table::fread(file, data.table = FALSE)
    cell_col <- intersect(c("cell", "barcode", "Cell", "Barcode"), colnames(df))[1]
    if (is.na(cell_col)) {
      warning("Skipping aptamer file without a cell/barcode column: ", file)
      next
    }

    cell <- df[[cell_col]]
    apt_df <- df[, setdiff(colnames(df), c(cell_col, "sample")), drop = FALSE]
    # Retain numeric columns only; this removes auxiliary columns if present.
    numeric_cols <- vapply(apt_df, is.numeric, logical(1))
    apt_df <- apt_df[, numeric_cols, drop = FALSE]
    colnames(apt_df) <- standardize_aptamer_name(colnames(apt_df))
    apt_df <- apt_df[, !duplicated(colnames(apt_df)), drop = FALSE]

    # Prefix sample label to match Seurat cell names. If cell names already contain
    # the label, avoid double-prefixing.
    prefixed_cell <- ifelse(grepl(paste0("^", label, "_"), cell), cell, paste(label, cell, sep = "_"))
    rownames(apt_df) <- prefixed_cell
    all_counts[[label]] <- apt_df
  }

  if (length(all_counts) == 0) stop("No ALL_APT_cell.csv files were found under: ", apt_root)

  all_features <- unique(unlist(lapply(all_counts, colnames)))
  all_features <- c(intersect(aptamer_order, all_features), setdiff(all_features, aptamer_order))

  padded <- lapply(all_counts, function(df) {
    missing <- setdiff(all_features, colnames(df))
    for (m in missing) df[[m]] <- 0
    df[, all_features, drop = FALSE]
  })
  mat <- do.call(rbind, padded)
  mat <- as.matrix(mat)
  mat <- t(mat)
  Matrix::Matrix(mat, sparse = TRUE)
}

add_apt_assay <- function(obj, apt_root, label_fun) {
  apt_mat <- read_aptamer_cell_counts(apt_root, label_fun)
  common_cells <- intersect(colnames(obj), colnames(apt_mat))

  if (length(common_cells) == 0) {
    # Fallback: try matching by unprefixed barcode suffix.
    obj_suffix <- sub("^[^_]+_", "", colnames(obj))
    apt_suffix <- sub("^[^_]+_", "", colnames(apt_mat))
    idx <- match(obj_suffix, apt_suffix)
    keep <- !is.na(idx)
    if (sum(keep) == 0) {
      stop("No overlapping cell barcodes between RNA and aptamer matrices.")
    }
    apt_mat <- apt_mat[, idx[keep], drop = FALSE]
    colnames(apt_mat) <- colnames(obj)[keep]
    obj <- obj[, keep]
    common_cells <- colnames(obj)
  }

  obj <- obj[, common_cells]
  apt_mat <- apt_mat[, common_cells]
  obj[["APT"]] <- Seurat::CreateAssayObject(counts = apt_mat)
  Seurat::DefaultAssay(obj) <- "APT"
  obj <- Seurat::NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)
  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}

filter_and_cluster_rna <- function(obj, prefix) {
  obj <- subset(obj, subset = nFeature_RNA >= MIN_GENES & nFeature_RNA <= MAX_GENES & percent.mt <= MAX_MT)

  if (RUN_DOUBLET) {
    # DoubletFinder requires a preliminary clustering. This block is optional
    # because package versions vary across systems.
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = max(30, N_PCS), verbose = FALSE)
    obj <- Seurat::FindNeighbors(obj, dims = 1:N_PCS, verbose = FALSE)
    obj <- Seurat::FindClusters(obj, resolution = RESOLUTION, verbose = FALSE)
    # A conservative placeholder: users can add their DoubletFinder version-specific
    # implementation here. Keeping this script robust is preferred for GitHub use.
    warning("RUN_DOUBLET requested; implement the DoubletFinder version-specific call if needed.")
  }

  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = max(30, N_PCS), verbose = FALSE)

  if (HAS_HARMONY && length(unique(obj$sample_group)) > 1) {
    obj <- harmony::RunHarmony(obj, group.by.vars = "sample_group", reduction = "pca", verbose = FALSE)
    reduction_use <- "harmony"
  } else {
    reduction_use <- "pca"
  }

  obj <- Seurat::RunUMAP(obj, reduction = reduction_use, dims = 1:N_PCS, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, reduction = reduction_use, dims = 1:N_PCS, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = RESOLUTION, algorithm = 2, verbose = FALSE)
  obj$cluster_id <- as.character(Seurat::Idents(obj))
  saveRDS(obj, file.path(OUT_DIR, "objects", paste0(prefix, "_seurat_with_APT.rds")))
  obj
}

feature_present <- function(obj, assay, feature) {
  feature %in% rownames(obj[[assay]])
}

safe_feature_plot <- function(obj, feature, assay, title, high_color, low_color = "#D9D9D9") {
  if (!feature_present(obj, assay, feature)) {
    return(ggplot() + theme_void() + ggtitle(paste0(title, "\n(not detected)")))
  }
  Seurat::DefaultAssay(obj) <- assay
  Seurat::FeaturePlot(obj, features = feature, reduction = "umap", cols = c(low_color, high_color), pt.size = 0.35) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right")
}

make_qc_summary <- function(obj, prefix) {
  summary_df <- obj@meta.data %>%
    dplyr::group_by(sample_group) %>%
    dplyr::summarise(
      cell_number = dplyr::n(),
      median_gene_number = median(nFeature_RNA, na.rm = TRUE),
      median_umi = median(nCount_RNA, na.rm = TRUE),
      median_mitochondrial_percentage = median(percent.mt, na.rm = TRUE),
      median_aptamer_counts = if ("nCount_APT" %in% colnames(obj@meta.data)) median(nCount_APT, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  data.table::fwrite(summary_df, file.path(OUT_DIR, "tables", paste0(prefix, "_QC_summary.csv")))
  summary_df
}

# ------------------------------- Figure 4a ----------------------------------
plot_qc_violin <- function(obj) {
  meta <- obj@meta.data
  meta$sample_group <- factor(meta$sample_group, levels = intersect(names(group_palette_hct8), unique(meta$sample_group)))
  qc_features <- data.frame(
    variable = c("nCount_APT", "nFeature_RNA", "nCount_RNA", "percent.mt"),
    label = c("APT Counts", "Gene Counts", "UMI Counts", "Percent MT gene")
  )

  plots <- lapply(seq_len(nrow(qc_features)), function(i) {
    var <- qc_features$variable[i]
    label <- qc_features$label[i]
    ggplot(meta, aes(x = sample_group, y = .data[[var]], fill = sample_group)) +
      geom_violin(trim = TRUE, linewidth = 0.25) +
      scale_fill_manual(values = group_palette_hct8, drop = FALSE) +
      labs(title = label, x = NULL, y = ifelse(i == 1, "Counts", NULL)) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  })
  patchwork::wrap_plots(plots, nrow = 1) + patchwork::plot_annotation(tag_levels = NULL)
}

# ------------------------------- Figure 4b ----------------------------------
plot_hct8_umap_panel <- function(obj) {
  p_cluster <- Seurat::DimPlot(obj, reduction = "umap", group.by = "cluster_id", label = TRUE,
                               repel = TRUE, cols = cluster_palette, pt.size = 0.35) +
    ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_group <- Seurat::DimPlot(obj, reduction = "umap", group.by = "sample_group", cols = group_palette_hct8,
                             pt.size = 0.35) +
    ggtitle("Groups") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_ptk7 <- safe_feature_plot(obj, "PTK7", "RNA", "PTK7 RNA", high_color = "#24469C")
  p_sgc8c <- safe_feature_plot(obj, "Sgc8c_Aptamer", "APT", "Sgc8c Aptamer", high_color = "#1B6735")
  p_scramble <- safe_feature_plot(obj, "Scramble_Aptamer", "APT", "Scramble Aptamer", high_color = "#1B6735")

  (p_cluster / p_group) | (p_ptk7 / p_sgc8c / p_scramble)
}

# ------------------------------- Figure 4c ----------------------------------
plot_aptamer_proportions <- function(obj, prefix = "HCT8") {
  Seurat::DefaultAssay(obj) <- "APT"
  apt_counts <- as.matrix(Seurat::GetAssayData(obj, assay = "APT", slot = "counts"))
  features <- intersect(aptamer_order, rownames(apt_counts))
  if (length(features) == 0) stop("No expected aptamer features were found in the APT assay.")

  sample_group <- obj$sample_group
  df <- lapply(split(seq_along(sample_group), sample_group), function(idx) {
    sums <- rowSums(apt_counts[features, idx, drop = FALSE])
    data.frame(sample_group = sample_group[idx][1], aptamer = names(sums), counts = as.numeric(sums))
  }) %>% dplyr::bind_rows()

  df <- df %>%
    dplyr::group_by(sample_group) %>%
    dplyr::mutate(proportion = counts / sum(counts)) %>%
    dplyr::ungroup()

  df$sample_group <- factor(df$sample_group, levels = names(group_palette_hct8))
  df$aptamer <- factor(df$aptamer, levels = rev(aptamer_order))
  data.table::fwrite(df, file.path(OUT_DIR, "tables", paste0(prefix, "_aptamer_proportions.csv")))

  ggplot(df, aes(x = sample_group, y = proportion, fill = aptamer)) +
    geom_col(color = "black", linewidth = 0.25, width = 0.85) +
    scale_fill_manual(values = aptamer_palette, drop = FALSE, name = "Aptamer") +
    labs(title = "Aptamer Proportions in Each Group", x = "Groups", y = "Proportion") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
}

# ------------------------------ Supplement S5 --------------------------------
plot_all_feature_pairs <- function(obj, prefix = "HCT8") {
  pair_plots <- list()
  for (i in seq_len(nrow(aptamer_pairs))) {
    gene <- aptamer_pairs$target_gene[i]
    apt <- aptamer_pairs$aptamer_feature[i]
    target <- aptamer_pairs$target_name[i]
    pair_plots[[paste0(target, "_RNA")]] <- safe_feature_plot(obj, gene, "RNA", paste0(target, " RNA"), high_color = "#24469C")
    pair_plots[[paste0(target, "_APT")]] <- safe_feature_plot(obj, apt, "APT", paste0(aptamer_pairs$aptamer_id[i], " Aptamer"), high_color = "#1B6735")
  }
  p <- patchwork::wrap_plots(pair_plots, ncol = 6)
  save_plot(p, paste0(prefix, "_Supplement_feature_plots_all_targets"), width = 18, height = 10)
  p
}

# ------------------------------ Supplement S6 --------------------------------
plot_cluster_correlations <- function(obj, prefix = "HCT8") {
  Seurat::DefaultAssay(obj) <- "RNA"
  cluster_var <- "cluster_id"
  rna_means <- Seurat::AverageExpression(obj, assays = "RNA", group.by = cluster_var, slot = "data", verbose = FALSE)$RNA
  apt_means <- Seurat::AverageExpression(obj, assays = "APT", group.by = cluster_var, slot = "data", verbose = FALSE)$APT

  cor_summary <- list()
  plots <- list()
  for (i in seq_len(nrow(aptamer_pairs))) {
    gene <- aptamer_pairs$target_gene[i]
    apt <- aptamer_pairs$aptamer_feature[i]
    if (!(gene %in% rownames(rna_means)) || !(apt %in% rownames(apt_means))) next
    df <- data.frame(
      RNA = as.numeric(rna_means[gene, ]),
      Aptamer = as.numeric(apt_means[apt, ]),
      cluster = colnames(rna_means)
    )
    cor_value <- if (nrow(df) >= 3 && stats::sd(df$RNA) > 0 && stats::sd(df$Aptamer) > 0) {
      stats::cor(df$RNA, df$Aptamer, method = "pearson")
    } else {
      NA_real_
    }
    cor_summary[[i]] <- data.frame(
      target_gene = gene,
      aptamer = apt,
      pearson_r = cor_value,
      n_clusters = nrow(df)
    )
    plots[[apt]] <- ggplot(df, aes(x = RNA, y = Aptamer)) +
      geom_point(size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "black", fill = "#3F693B", alpha = 0.2, linewidth = 0.6) +
      labs(
        title = paste0(gene, " / ", aptamer_pairs$aptamer_id[i]),
        subtitle = paste0("r = ", ifelse(is.na(cor_value), "NA", sprintf("%.3f", cor_value))),
        x = "Mean RNA expression",
        y = "Mean CLR aptamer signal"
      ) +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  cor_df <- dplyr::bind_rows(cor_summary)
  data.table::fwrite(cor_df, file.path(OUT_DIR, "tables", paste0(prefix, "_cluster_level_RNA_APT_correlations.csv")))
  p <- patchwork::wrap_plots(plots, ncol = 3)
  save_plot(p, paste0(prefix, "_Supplement_cluster_level_correlations"), width = 12, height = 11)
  p
}

# ------------------------------- Figure 4d/e --------------------------------
plot_ko_umap_panel <- function(obj) {
  ko_palette_use <- group_palette_ko[names(group_palette_ko) %in% unique(obj$sample_group)]

  p_cluster <- Seurat::DimPlot(obj, reduction = "umap", group.by = "cluster_id", label = TRUE,
                               repel = TRUE, cols = cluster_palette, pt.size = 0.35) +
    ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_group <- Seurat::DimPlot(obj, reduction = "umap", group.by = "sample_group", cols = ko_palette_use,
                             pt.size = 0.35) +
    ggtitle("Groups") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_ptk7 <- safe_feature_plot(obj, "PTK7", "RNA", "PTK7 RNA", high_color = "#24469C")
  p_sgc8c <- safe_feature_plot(obj, "Sgc8c_Aptamer", "APT", "Sgc8c Aptamer", high_color = "#1B6735")
  p_scramble <- safe_feature_plot(obj, "Scramble_Aptamer", "APT", "Scramble Aptamer", high_color = "#1B6735")

  (p_cluster / p_group) | (p_ptk7 / p_sgc8c / p_scramble)
}

plot_ko_boxplots <- function(obj) {
  groups_to_keep <- intersect(c("NC", "KO", "PTK7-KO"), unique(obj$sample_group))
  if (length(groups_to_keep) < 2) {
    warning("NC and KO groups were not both detected; Figure 4e box plot will be skipped.")
    return(NULL)
  }
  obj2 <- obj[, obj$sample_group %in% groups_to_keep]
  obj2$sample_group <- ifelse(obj2$sample_group == "PTK7-KO", "KO", obj2$sample_group)
  obj2$sample_group <- factor(obj2$sample_group, levels = c("KO", "NC"))

  Seurat::DefaultAssay(obj2) <- "APT"
  df <- Seurat::FetchData(obj2, vars = c("Sgc8c_Aptamer", "Scramble_Aptamer", "sample_group"))
  data.table::fwrite(df, file.path(OUT_DIR, "tables", "KO_NC_sgc8c_scramble_counts.csv"))

  make_box <- function(feature, color, title, ylab) {
    p <- ggplot(df, aes(x = sample_group, y = .data[[feature]])) +
      geom_boxplot(width = 0.55, outlier.shape = NA, fill = NA, color = color, linewidth = 0.8) +
      labs(title = title, x = NULL, y = ylab) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    if (HAS_GGSIGNIF && all(c("KO", "NC") %in% df$sample_group)) {
      ymax <- max(df[[feature]], na.rm = TRUE)
      if (is.finite(ymax) && ymax > 0) {
        p <- p + ggsignif::geom_signif(
          comparisons = list(c("KO", "NC")),
          test = "wilcox.test",
          map_signif_level = TRUE,
          y_position = ymax * 1.08,
          tip_length = 0.01,
          color = color
        )
      }
    }
    p
  }

  p1 <- make_box("Sgc8c_Aptamer", "#7BC061", "Sgc8c", "Counts")
  p2 <- make_box("Scramble_Aptamer", "#80B1D3", "Scramble", "Counts")
  p1 / p2
}

# ------------------------------ main workflow --------------------------------
message("\n=== Processing HCT-8 concentration-test scApt-seq dataset ===")
hct8 <- read_rna_samples(HCT8_RNA_DIR, standardize_sample_label_hct8)
hct8 <- add_apt_assay(hct8, HCT8_APT_DIR, standardize_sample_label_hct8)
hct8 <- filter_and_cluster_rna(hct8, "HCT8_Figure4")
make_qc_summary(hct8, "HCT8_Figure4")

fig4a <- plot_qc_violin(hct8)
save_plot(fig4a, "Figure4a_HCT8_single_cell_QC", width = 13, height = 3.5)

fig4b <- plot_hct8_umap_panel(hct8)
save_plot(fig4b, "Figure4b_HCT8_UMAP_RNA_APT_scramble", width = 10, height = 8)

fig4c <- plot_aptamer_proportions(hct8, prefix = "HCT8_Figure4")
save_plot(fig4c, "Figure4c_HCT8_aptamer_proportions", width = 7, height = 5.5)

plot_all_feature_pairs(hct8, prefix = "HCT8_Figure4")
plot_cluster_correlations(hct8, prefix = "HCT8_Figure4")

# Optional PTK7-knockout dataset. This reproduces Figure 4d/e when KO_RNA_DIR and
# KO_APT_DIR are provided. The manuscript uses this dataset to support the
# target-dependent recognition of sgc8c after PTK7 deletion.
if (!identical(KO_RNA_DIR, "") && !identical(KO_APT_DIR, "") && dir.exists(KO_RNA_DIR) && dir.exists(KO_APT_DIR)) {
  message("\n=== Processing PTK7-knockout specificity-control scApt-seq dataset ===")
  ko_obj <- read_rna_samples(KO_RNA_DIR, standardize_sample_label_ko)
  ko_obj <- add_apt_assay(ko_obj, KO_APT_DIR, standardize_sample_label_ko)
  ko_obj <- filter_and_cluster_rna(ko_obj, "PTK7_KO_Figure4")
  make_qc_summary(ko_obj, "PTK7_KO_Figure4")

  fig4d <- plot_ko_umap_panel(ko_obj)
  save_plot(fig4d, "Figure4d_PTK7_KO_UMAP_RNA_APT_scramble", width = 10, height = 8)

  fig4e <- plot_ko_boxplots(ko_obj)
  if (!is.null(fig4e)) {
    save_plot(fig4e, "Figure4e_PTK7_KO_sgc8c_scramble_boxplots", width = 7, height = 7.5)
  }
} else {
  message("\nKO_RNA_DIR and/or KO_APT_DIR were not provided. Skipping Figure 4d/e outputs.")
}

# Save a compact manifest for GitHub reproducibility.
manifest <- data.frame(
  item = c("HCT8_RNA_DIR", "HCT8_APT_DIR", "KO_RNA_DIR", "KO_APT_DIR", "OUT_DIR", "MIN_GENES", "MAX_GENES", "MAX_MT", "N_PCS", "CLUSTER_RESOLUTION"),
  value = c(HCT8_RNA_DIR, HCT8_APT_DIR, KO_RNA_DIR, KO_APT_DIR, OUT_DIR, MIN_GENES, MAX_GENES, MAX_MT, N_PCS, RESOLUTION)
)
data.table::fwrite(manifest, file.path(OUT_DIR, "tables", "run_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

message("\nFinished. Outputs were written to: ", normalizePath(OUT_DIR))
