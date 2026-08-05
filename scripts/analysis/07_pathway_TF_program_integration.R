## =========================================================================
## Step 07: Pathway Enrichment, PROGENy, CollecTRI & Program Integration
## =========================================================================

.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))
options(stringsAsFactors = FALSE)

cat("========================================================================\n")
cat("Step 07: Pathway, TF & Program Integration\n")
cat("========================================================================\n\n")
flush.console()

step_start <- Sys.time()

## =========================================================================
## IV. Load Packages
## =========================================================================
cat("[IV] Loading packages...\n")
flush.console()

required_pkgs <- c("fgsea", "BiocParallel", "decoupleR", "data.table",
                    "dplyr", "tidyr", "tibble", "ggplot2", "pheatmap",
                    "msigdbr", "RColorBrewer")
lib <- path.expand("~/Library/R/arm64/4.6/library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, function(p) {
  requireNamespace(p, lib.loc = lib, quietly = TRUE)
})]

if (length(missing_pkgs) > 0) {
  cat(sprintf("[INFO] Missing packages: %s\n", paste(missing_pkgs, collapse = ", ")))
  bioc_miss <- intersect(missing_pkgs, c("fgsea", "BiocParallel", "decoupleR"))
  if (length(bioc_miss) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org",
                       type = "binary", lib = lib)
    BiocManager::install(bioc_miss, ask = FALSE, update = FALSE, lib = lib)
  }
  cran_miss <- setdiff(missing_pkgs, c("fgsea", "BiocParallel", "decoupleR"))
  if (length(cran_miss) > 0)
    install.packages(cran_miss, repos = "https://cloud.r-project.org",
                     type = "binary", lib = lib)
}

suppressPackageStartupMessages({
  library(fgsea)
  library(BiocParallel)
  library(decoupleR)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(msigdbr)
  library(RColorBrewer)
})

env_info <- c(
  sprintf("R version: %s", R.version.string),
  sprintf("R platform: %s", R.version$platform),
  sprintf("fgsea version: %s", packageVersion("fgsea")),
  sprintf("decoupleR version: %s", packageVersion("decoupleR")),
  sprintf("msigdbr version: %s", packageVersion("msigdbr")),
  sprintf("BiocParallel version: %s", packageVersion("BiocParallel")),
  sprintf("libPaths: %s", paste(.libPaths(), collapse = ", ")),
  sprintf("Date: %s", Sys.time())
)

dir.create("03_results/step07_programs/qc", recursive = TRUE, showWarnings = FALSE)
writeLines(env_info, "03_results/step07_programs/qc/GSE243013_step07_environment.txt")
cat("[INFO] Saved environment info\n")
cat(sprintf("[INFO] fgsea v%s, decoupleR v%s, msigdbr v%s\n",
            packageVersion("fgsea"), packageVersion("decoupleR"), packageVersion("msigdbr")))

## =========================================================================
## V. Freeze Step 06 Results
## =========================================================================
cat("\n[V] Freezing Step 06 results...\n")

if (!file.exists("03_results/GSE243013_step06_COMPLETE.txt")) {
  stop("[FATAL] Step 06 not complete. Cannot proceed.")
}
cat("[INFO] Step 06 COMPLETE confirmed.\n")

step06_files <- c(
  "03_results/step06_edgeR/combined/GSE243013_edgeR_model_status.csv",
  "03_results/step06_edgeR/combined/GSE243013_primary_edgeR_summary.csv",
  "03_results/step06_edgeR/combined/GSE243013_primary_edgeR_all_celltypes_all_genes.csv.gz",
  "03_results/step06_edgeR/combined/GSE243013_primary_TREAT_all_celltypes_all_genes.csv.gz"
)

missing_s6 <- step06_files[!file.exists(step06_files)]
if (length(missing_s6) > 0) {
  stop(sprintf("[FATAL] Missing Step 06 files: %s", paste(missing_s6, collapse = ", ")))
}

freeze_rows <- list()
for (f in step06_files) {
  fi <- file.info(f)
  md5 <- tryCatch(tools::md5sum(f), error = function(e) NA_character_)
  freeze_rows[[length(freeze_rows) + 1]] <- data.frame(
    file_path = f, file_size = fi$size, modification_time = as.character(fi$mtime),
    md5 = md5, exists = TRUE, validation_status = "OK", stringsAsFactors = FALSE
  )
}
freeze_df <- do.call(rbind, freeze_rows)
write.csv(freeze_df, "03_results/step07_programs/qc/GSE243013_step06_frozen_file_manifest.csv",
          row.names = FALSE)
cat("[INFO] Step 06 frozen.\n")

## =========================================================================
## VI. Select Analyzable Models
## =========================================================================
cat("\n[VI] Selecting analyzable models...\n")

model_status <- read.csv("03_results/step06_edgeR/combined/GSE243013_edgeR_model_status.csv",
                         stringsAsFactors = FALSE)
cat(sprintf("[INFO] Model status columns: %s\n", paste(colnames(model_status), collapse = ", ")))

pri_complete <- model_status[model_status$analysis == "primary_anti_PD1" &
                              model_status$status == "COMPLETE", ]
pri_skipped <- model_status[model_status$analysis == "primary_anti_PD1" &
                             model_status$status != "COMPLETE", ]

cat(sprintf("[INFO] Primary COMPLETE: %d\n", nrow(pri_complete)))
cat(sprintf("[INFO] Primary SKIPPED/FAILED: %d\n", nrow(pri_skipped)))

## Build input status
input_rows <- list()
for (i in seq_len(nrow(model_status))) {
  row <- model_status[i, ]
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", row$cell_type)
  qlf_path <- sprintf("03_results/step06_edgeR/primary_anti_PD1/%s__edgeR_QLF_all_genes.csv.gz", ct_safe)
  input_rows[[length(input_rows) + 1]] <- data.frame(
    analysis_name = row$analysis, cell_type = row$cell_type,
    cell_type_safe = ct_safe, model_status = row$status,
    n_samples = row$n_samples, n_resp = row$n_resp, n_nonresp = row$n_nonresp,
    n_genes_after = row$n_genes_after,
    qlf_file_exists = file.exists(qlf_path),
    stringsAsFactors = FALSE
  )
}
input_status_df <- do.call(rbind, input_rows)
write.csv(input_status_df, "03_results/step07_programs/qc/GSE243013_step07_model_input_status.csv",
          row.names = FALSE)

## =========================================================================
## VII. Build Gene Rankings
## =========================================================================
cat("\n[VII] Building gene rankings...\n")
flush.console()

rank_qc_rows <- list()
primary_ranks <- list()
primary_analyzed <- character(0)

for (i in seq_len(nrow(pri_complete))) {
  row <- pri_complete[i, ]
  ct <- row$cell_type
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  analysis_name <- "primary_anti_PD1"

  cat(sprintf("\n--- Ranking: %s ---\n", ct))
  flush.console()

  qlf_path <- sprintf("03_results/step06_edgeR/primary_anti_PD1/%s__edgeR_QLF_all_genes.csv.gz", ct_safe)
  if (!file.exists(qlf_path)) {
    cat(sprintf("[SKIP] QLF file not found: %s\n", qlf_path))
    rank_qc_rows[[length(rank_qc_rows) + 1]] <- data.frame(
      analysis_name = analysis_name, cell_type = ct, n_input_genes = NA_integer_,
      n_valid_genes = NA_integer_, n_removed = NA_integer_, n_duplicate_genes = NA_integer_,
      n_tied_primary_stats = NA_integer_, min_rank = NA_real_, max_rank = NA_real_,
      median_rank = NA_real_, positive_gene_count = NA_integer_, negative_gene_count = NA_integer_,
      rank_validation_status = "SKIPPED_NO_QLF_FILE", stringsAsFactors = FALSE
    )
    next
  }

  qlf <- read.csv(qlf_path, stringsAsFactors = FALSE)

  ## Validate required columns
  required_cols <- c("gene", "logFC", "F", "PValue")
  missing_cols <- setdiff(required_cols, colnames(qlf))
  if (length(missing_cols) > 0) {
    cat(sprintf("[SKIP] Missing columns: %s\n", paste(missing_cols, collapse = ", ")))
    rank_qc_rows[[length(rank_qc_rows) + 1]] <- data.frame(
      analysis_name = analysis_name, cell_type = ct, n_input_genes = nrow(qlf),
      n_valid_genes = NA_integer_, n_removed = NA_integer_, n_duplicate_genes = NA_integer_,
      n_tied_primary_stats = NA_integer_, min_rank = NA_real_, max_rank = NA_real_,
      median_rank = NA_real_, positive_gene_count = NA_integer_, negative_gene_count = NA_integer_,
      rank_validation_status = "SKIPPED_MISSING_COLUMNS", stringsAsFactors = FALSE
    )
    next
  }

  n_input <- nrow(qlf)

  ## Remove rows with missing/empty gene or non-finite stats
  valid <- !is.na(qlf$gene) & qlf$gene != "" &
           is.finite(qlf$logFC) & is.finite(qlf$F) & is.finite(qlf$PValue)
  qlf <- qlf[valid, ]

  ## Remove duplicate genes (keep first)
  dup_genes <- qlf$gene[duplicated(qlf$gene)]
  if (length(dup_genes) > 0) {
    cat(sprintf("[WARNING] %d duplicate genes found\n", length(dup_genes)))
    qlf <- qlf[!duplicated(qlf$gene), ]
  }

  ## Create rank statistics
  qlf$rank_signed_sqrtF <- sign(qlf$logFC) * sqrt(pmax(qlf$F, 0))
  qlf$rank_logFC_logP <- qlf$logFC * -log10(pmax(qlf$PValue, 1e-300))

  ## Sort by primary rank (descending), then logFC, then gene alphabetically
  qlf <- qlf[order(-qlf$rank_signed_sqrtF, -qlf$logFC, qlf$gene), ]

  ## Check for ties
  n_tied <- sum(duplicated(qlf$rank_signed_sqrtF))

  ## Create rank vectors
  rank_primary <- setNames(qlf$rank_signed_sqrtF, qlf$gene)
  rank_sensitivity <- setNames(qlf$rank_logFC_logP, qlf$gene)

  ## Save ranks
  rank_df <- data.frame(
    gene = qlf$gene, rank_signed_sqrtF = qlf$rank_signed_sqrtF,
    rank_logFC_logP = qlf$rank_logFC_logP, logFC = qlf$logFC,
    F = qlf$F, PValue = qlf$PValue,
    stringsAsFactors = FALSE
  )
  rank_path <- sprintf("03_results/step07_programs/ranks/%s__gene_ranks.csv.gz", ct_safe)
  con <- gzfile(rank_path, "w")
  write.csv(rank_df, con, row.names = FALSE)
  close(con)

  primary_ranks[[ct]] <- list(
    primary = rank_primary, sensitivity = rank_sensitivity,
    analysis_name = analysis_name, cell_type = ct, cell_type_safe = ct_safe
  )
  primary_analyzed <- c(primary_analyzed, ct)

  rank_qc_rows[[length(rank_qc_rows) + 1]] <- data.frame(
    analysis_name = analysis_name, cell_type = ct,
    n_input_genes = n_input, n_valid_genes = nrow(qlf),
    n_removed = n_input - nrow(qlf), n_duplicate_genes = length(dup_genes),
    n_tied_primary_stats = n_tied,
    min_rank = min(qlf$rank_signed_sqrtF), max_rank = max(qlf$rank_signed_sqrtF),
    median_rank = median(qlf$rank_signed_sqrtF),
    positive_gene_count = sum(qlf$rank_signed_sqrtF > 0),
    negative_gene_count = sum(qlf$rank_signed_sqrtF < 0),
    rank_validation_status = "OK", stringsAsFactors = FALSE
  )
  cat(sprintf("[OK] %s: %d genes ranked (%d positive, %d negative)\n",
              ct, nrow(qlf), sum(qlf$rank_signed_sqrtF > 0), sum(qlf$rank_signed_sqrtF < 0)))
}

rank_qc_df <- do.call(rbind, rank_qc_rows)
write.csv(rank_qc_df, "03_results/step07_programs/qc/GSE243013_gene_rank_qc.csv", row.names = FALSE)
cat(sprintf("\n[INFO] Gene rankings built for %d cell types\n", length(primary_analyzed)))

## =========================================================================
## VIII. Prepare MSigDB Gene Sets
## =========================================================================
cat("\n[VIII] Preparing MSigDB gene sets...\n")
flush.console()

hs <- "Homo sapiens"

## Hallmark
hallmark_path <- "00_config/step07_resources/MSigDB_Hallmark.rds"
if (file.exists(hallmark_path)) {
  hallmark_gsets <- readRDS(hallmark_path)
  cat("[INFO] Loaded cached Hallmark gene sets\n")
} else {
  cat("[INFO] Downloading Hallmark gene sets from MSigDB via msigdbr...\n")
  hallmark_raw <- msigdbr(species = hs, collection = "H")
  hallmark_df <- hallmark_raw[, c("gs_id", "gs_name", "gene_symbol", "gs_description")]
  hallmark_df$gs_collection <- "H"
  hallmark_df$gs_subcollection <- "HALLMARK"
  hallmark_df <- hallmark_df[!is.na(hallmark_df$gene_symbol) & !is.na(hallmark_df$gs_name), ]
  hallmark_df <- hallmark_df[!duplicated(paste(hallmark_df$gs_name, hallmark_df$gene_symbol)), ]
  hallmark_gsets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
  saveRDS(hallmark_gsets, hallmark_path)
  cat(sprintf("[INFO] Hallmark: %d gene sets, %d memberships\n",
              length(hallmark_gsets), nrow(hallmark_df)))
}

## Reactome
reactome_path <- "00_config/step07_resources/MSigDB_Reactome.rds"
if (file.exists(reactome_path)) {
  reactome_gsets <- readRDS(reactome_path)
  cat("[INFO] Loaded cached Reactome gene sets\n")
} else {
  cat("[INFO] Downloading Reactome gene sets from MSigDB via msigdbr...\n")
  reactome_raw <- msigdbr(species = hs, collection = "C2", subcollection = "CP:REACTOME")
  reactome_df <- reactome_raw[, c("gs_id", "gs_name", "gene_symbol", "gs_description")]
  reactome_df$gs_collection <- "C2"
  reactome_df$gs_subcollection <- "CP:REACTOME"
  reactome_df <- reactome_df[!is.na(reactome_df$gene_symbol) & !is.na(reactome_df$gs_name), ]
  reactome_df <- reactome_df[!duplicated(paste(reactome_df$gs_name, reactome_df$gene_symbol)), ]
  reactome_gsets <- split(reactome_df$gene_symbol, reactome_df$gs_name)
  saveRDS(reactome_gsets, reactome_path)
  cat(sprintf("[INFO] Reactome: %d gene sets, %d memberships\n",
              length(reactome_gsets), nrow(reactome_df)))
}

## Save metadata
msigdb_meta <- data.frame(
  db_version = "2024.1.Hs",
  collection = c("H", "C2"),
  subcollection = c("HALLMARK", "CP:REACTOME"),
  gene_set_count = c(length(hallmark_gsets), length(reactome_gsets)),
  gene_membership_count = c(
    sum(sapply(hallmark_gsets, length)),
    sum(sapply(reactome_gsets, length))
  ),
  retrieval_time = as.character(Sys.time()),
  msigdbr_version = as.character(packageVersion("msigdbr")),
  stringsAsFactors = FALSE
)
write.csv(msigdb_meta, "00_config/step07_resources/MSigDB_resource_metadata.csv", row.names = FALSE)

## =========================================================================
## IX. Run fgsea (Primary)
## =========================================================================
cat("\n[IX] Running fgsea on primary analysis...\n")
flush.console()

BPPARAM <- tryCatch({
  BiocParallel::MulticoreParam(workers = 2, RNGseed = 20260804, progressbar = FALSE)
}, error = function(e) {
  cat(sprintf("[WARNING] MulticoreParam failed, using SerialParam: %s\n", e$message))
  BiocParallel::SerialParam()
})

run_one_fgsea <- function(rank_vector, pathways, analysis_name, cell_type,
                          collection_name, output_path, rank_type) {
  if (length(rank_vector) < 50) {
    cat(sprintf("[SKIP] Too few ranked genes (%d)\n", length(rank_vector)))
    return(data.frame())
  }

  ## Filter pathways to those with overlap
  pathway_sizes <- sapply(pathways, function(ps) sum(ps %in% names(rank_vector)))
  valid_pathways <- pathways[pathway_sizes >= 15 & pathway_sizes <= 500]

  if (length(valid_pathways) == 0) {
    cat(sprintf("[SKIP] No valid pathways for %s/%s\n", cell_type, collection_name))
    return(data.frame())
  }

  cat(sprintf("[INFO] Running fgsea: %d pathways, %d genes\n", length(valid_pathways), length(rank_vector)))

  fgsea_res <- fgsea::fgseaMultilevel(
    pathways = valid_pathways,
    stats = rank_vector,
    minSize = 15,
    maxSize = 500,
    eps = 1e-50,
    scoreType = "std",
    gseaParam = 1,
    BPPARAM = BPPARAM
  )

  if (nrow(fgsea_res) == 0) return(data.frame())

  fgsea_res <- as.data.frame(fgsea_res)
  fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, function(x) paste(x, collapse = "; "))
  fgsea_res$analysis_name <- analysis_name
  fgsea_res$cell_type <- cell_type
  fgsea_res$collection <- collection_name
  fgsea_res$rank_type <- rank_type
  fgsea_res$positive_NES_direction <- "Higher_in_Responder"
  fgsea_res$negative_NES_direction <- "Higher_in_Non_responder"
  fgsea_res$n_ranked_genes <- length(rank_vector)
  fgsea_res$fgsea_version <- as.character(packageVersion("fgsea"))
  fgsea_res$MSigDB_version <- "2024.1.Hs"

  ## Save
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(output_path, "w")
  write.csv(fgsea_res, con, row.names = FALSE)
  close(con)

  n_sig <- sum(fgsea_res$padj < 0.05, na.rm = TRUE)
  cat(sprintf("[OK] %s/%s: %d pathways tested, %d sig (padj<0.05)\n",
              cell_type, collection_name, nrow(fgsea_res), n_sig))
  fgsea_res
}

## Run Hallmark + Reactome for each primary cell type
all_hallmark_primary <- list()
all_reactome_primary <- list()
all_hallmark_sensitivity <- list()

for (ct in primary_analyzed) {
  rk <- primary_ranks[[ct]]
  ct_safe <- rk$cell_type_safe

  cat(sprintf("\n--- fgsea: %s ---\n", ct))

  ## Hallmark primary
  hm_path <- sprintf("03_results/step07_programs/fgsea_primary/%s__Hallmark_fgsea.csv.gz", ct_safe)
  hm_res <- run_one_fgsea(
    rank_vector = rk$primary, pathways = hallmark_gsets,
    analysis_name = "primary_anti_PD1", cell_type = ct,
    collection_name = "Hallmark", output_path = hm_path,
    rank_type = "rank_signed_sqrtF"
  )
  if (nrow(hm_res) > 0) all_hallmark_primary[[ct]] <- hm_res

  ## Reactome primary
  rx_path <- sprintf("03_results/step07_programs/fgsea_primary/%s__Reactome_fgsea.csv.gz", ct_safe)
  rx_res <- run_one_fgsea(
    rank_vector = rk$primary, pathways = reactome_gsets,
    analysis_name = "primary_anti_PD1", cell_type = ct,
    collection_name = "Reactome", output_path = rx_path,
    rank_type = "rank_signed_sqrtF"
  )
  if (nrow(rx_res) > 0) all_reactome_primary[[ct]] <- rx_res

  ## Hallmark sensitivity ranking
  hs_path <- sprintf("03_results/step07_programs/fgsea_primary/%s__Hallmark_rank_sensitivity.csv.gz", ct_safe)
  hs_res <- run_one_fgsea(
    rank_vector = rk$sensitivity, pathways = hallmark_gsets,
    analysis_name = "primary_anti_PD1", cell_type = ct,
    collection_name = "Hallmark", output_path = hs_path,
    rank_type = "rank_logFC_logP"
  )
  if (nrow(hs_res) > 0) all_hallmark_sensitivity[[ct]] <- hs_res

  ## Reactome sensitivity (only for significant pathways)
  if (nrow(rx_res) > 0) {
    sig_rx <- rx_res[rx_res$padj < 0.10 | abs(rx_res$NES) >= quantile(abs(rx_res$NES), 0.95, na.rm = TRUE), ]
    if (nrow(sig_rx) > 0) {
      rx_sens_path <- sprintf("03_results/step07_programs/fgsea_primary/%s__Reactome_rank_sensitivity.csv.gz", ct_safe)
      rx_sens_res <- run_one_fgsea(
        rank_vector = rk$sensitivity, pathways = reactome_gsets[sig_rx$pathway],
        analysis_name = "primary_anti_PD1", cell_type = ct,
        collection_name = "Reactome", output_path = rx_sens_path,
        rank_type = "rank_logFC_logP"
      )
    }
  }
}

## =========================================================================
## X. Strict Cohort GSEA
## =========================================================================
cat("\n[X] Running strict cohort GSEA...\n")
flush.console()

strict_dir <- "03_results/step06_edgeR/strict_chemoimmunotherapy"
strict_complete <- model_status[model_status$analysis == "strict_chemoimmunotherapy" &
                                 model_status$status == "COMPLETE", ]

all_hallmark_strict <- list()
all_reactome_strict <- list()

if (nrow(strict_complete) > 0) {
  for (i in seq_len(nrow(strict_complete))) {
    row <- strict_complete[i, ]
    ct <- row$cell_type
    ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)

    qlf_path <- sprintf("%s/%s__edgeR_QLF_all_genes.csv.gz", strict_dir, ct_safe)
    if (!file.exists(qlf_path)) {
      cat(sprintf("[SKIP] Strict QLF not found: %s\n", ct))
      next
    }

    cat(sprintf("\n--- Strict fgsea: %s ---\n", ct))
    qlf <- read.csv(qlf_path, stringsAsFactors = FALSE)

    if (!all(c("gene", "logFC", "F", "PValue") %in% colnames(qlf))) {
      cat("[SKIP] Missing required columns\n")
      next
    }

    valid <- !is.na(qlf$gene) & qlf$gene != "" & is.finite(qlf$logFC) & is.finite(qlf$F) & is.finite(qlf$PValue)
    qlf <- qlf[valid, ]
    qlf <- qlf[!duplicated(qlf$gene), ]
    rank_vec <- setNames(sign(qlf$logFC) * sqrt(pmax(qlf$F, 0)), qlf$gene)
    rank_vec <- sort(rank_vec, decreasing = TRUE)

    hm_path <- sprintf("03_results/step07_programs/fgsea_strict/%s__Hallmark_fgsea.csv.gz", ct_safe)
    hm_res <- run_one_fgsea(rank_vec, hallmark_gsets, "strict_chemoimmunotherapy", ct,
                            "Hallmark", hm_path, "rank_signed_sqrtF")
    if (nrow(hm_res) > 0) all_hallmark_strict[[ct]] <- hm_res

    rx_path <- sprintf("03_results/step07_programs/fgsea_strict/%s__Reactome_fgsea.csv.gz", ct_safe)
    rx_res <- run_one_fgsea(rank_vec, reactome_gsets, "strict_chemoimmunotherapy", ct,
                            "Reactome", rx_path, "rank_signed_sqrtF")
    if (nrow(rx_res) > 0) all_reactome_strict[[ct]] <- rx_res
  }
} else {
  cat("[INFO] No strict COMPLETE models found\n")
}

## =========================================================================
## XI. PROGENy and CollecTRI Resources
## =========================================================================
cat("\n[XI] Loading PROGENy and CollecTRI resources...\n")

progeny_path <- "00_config/step07_resources/PROGENy_human_top500.rds"
collectri_path <- "00_config/step07_resources/CollecTRI_human.rds"

progeny_net <- NULL
collectri_net <- NULL
progeny_available <- FALSE
collectri_available <- FALSE

## PROGENy
if (file.exists(progeny_path)) {
  progeny_net <- readRDS(progeny_path)
  progeny_available <- TRUE
  cat("[INFO] Loaded cached PROGENy network\n")
} else {
  for (attempt in 1:3) {
    cat(sprintf("[INFO] PROGENy download attempt %d/3...\n", attempt))
    progeny_net <- tryCatch(decoupleR::get_progeny(organism = "human", top = 500),
                            error = function(e) { cat(sprintf("[WARNING] %s\n", e$message)); NULL })
    if (!is.null(progeny_net)) {
      saveRDS(progeny_net, progeny_path)
      progeny_available <- TRUE
      cat("[INFO] PROGENy downloaded and cached\n")
      break
    }
    if (attempt < 3) Sys.sleep(10)
  }
}

## CollecTRI
if (file.exists(collectri_path)) {
  collectri_net <- readRDS(collectri_path)
  collectri_available <- TRUE
  cat("[INFO] Loaded cached CollecTRI network\n")
} else {
  for (attempt in 1:3) {
    cat(sprintf("[INFO] CollecTRI download attempt %d/3...\n", attempt))
    collectri_net <- tryCatch(decoupleR::get_collectri(organism = "human", split_complexes = FALSE),
                              error = function(e) { cat(sprintf("[WARNING] %s\n", e$message)); NULL })
    if (!is.null(collectri_net)) {
      saveRDS(collectri_net, collectri_path)
      collectri_available <- TRUE
      cat("[INFO] CollecTRI downloaded and cached\n")
      break
    }
    if (attempt < 3) Sys.sleep(10)
  }
}

## Audit resources
audit_rows <- list()
if (!is.null(progeny_net)) {
  audit_rows[[length(audit_rows) + 1]] <- data.frame(
    resource = "PROGENy", columns = paste(colnames(progeny_net), collapse = ", "),
    n_rows = nrow(progeny_net), status = "OK", stringsAsFactors = FALSE
  )
} else {
  audit_rows[[length(audit_rows) + 1]] <- data.frame(
    resource = "PROGENy", columns = NA, n_rows = NA_integer_,
    status = ifelse(file.exists(progeny_path), "CACHE_EXISTS_BUT_LOAD_FAILED", "RESOURCE_UNAVAILABLE"),
    stringsAsFactors = FALSE
  )
}
if (!is.null(collectri_net)) {
  audit_rows[[length(audit_rows) + 1]] <- data.frame(
    resource = "CollecTRI", columns = paste(colnames(collectri_net), collapse = ", "),
    n_rows = nrow(collectri_net), status = "OK", stringsAsFactors = FALSE
  )
} else {
  audit_rows[[length(audit_rows) + 1]] <- data.frame(
    resource = "CollecTRI", columns = NA, n_rows = NA_integer_,
    status = ifelse(file.exists(collectri_path), "CACHE_EXISTS_BUT_LOAD_FAILED", "RESOURCE_UNAVAILABLE"),
    stringsAsFactors = FALSE
  )
}
audit_df <- do.call(rbind, audit_rows)
write.csv(audit_df, "03_results/step07_programs/qc/GSE243013_decoupleR_resource_audit.csv", row.names = FALSE)

## =========================================================================
## XII. PROGENy Pathway Activity
## =========================================================================
cat("\n[XII] Running PROGENy pathway activity...\n")
flush.console()

all_progeny <- list()

if (progeny_available) {
  ## Determine parameter names
  prog_formals <- names(formals(decoupleR::run_mlm))
  use_network <- "network" %in% prog_formals

  for (ct in primary_analyzed) {
    rk <- primary_ranks[[ct]]
    ct_safe <- rk$cell_type_safe
    cat(sprintf("\n--- PROGENy: %s ---\n", ct))

    ## Build statistic matrix (gene x 1 condition)
    stat_mat <- matrix(rk$primary, ncol = 1, dimnames = list(names(rk$primary), ct_safe))

    prog_res <- tryCatch({
      if (use_network) {
        decoupleR::run_mlm(mat = stat_mat, network = progeny_net,
                           .source = "source", .target = "target", .mor = "weight", minsize = 5)
      } else {
        decoupleR::run_mlm(mat = stat_mat, net = progeny_net,
                           .source = "source", .target = "target", .mor = "weight", minsize = 5)
      }
    }, error = function(e) {
      cat(sprintf("[WARNING] PROGENy failed: %s\n", e$message))
      data.frame()
    })

    if (nrow(prog_res) > 0) {
      prog_res <- as.data.frame(prog_res)
      prog_res$analysis_name <- "primary_anti_PD1"
      prog_res$cell_type <- ct
      prog_res$n_targets_tested <- NA_integer_
      prog_res$positive_score_direction <- "More_active_in_Responder"
      prog_res$negative_score_direction <- "More_active_in_Non_responder"
      prog_res$FDR_within_celltype <- p.adjust(prog_res$p_value, method = "BH")

      out_path <- sprintf("03_results/step07_programs/progeny/%s__PROGENy_activity.csv", ct_safe)
      write.csv(prog_res, out_path, row.names = FALSE)
      all_progeny[[ct]] <- prog_res

      n_sig <- sum(prog_res$FDR_within_celltype < 0.05, na.rm = TRUE)
      cat(sprintf("[OK] %s: %d pathways, %d sig\n", ct, nrow(prog_res), n_sig))
    }
  }
} else {
  cat("[INFO] PROGENy not available, skipping\n")
}

## =========================================================================
## XIII. CollecTRI TF Activity
## =========================================================================
cat("\n[XIII] Running CollecTRI TF activity...\n")
flush.console()

all_collectri <- list()

if (collectri_available) {
  ctr_formals <- names(formals(decoupleR::run_ulm))
  use_network_ct <- "network" %in% ctr_formals

  for (ct in primary_analyzed) {
    rk <- primary_ranks[[ct]]
    ct_safe <- rk$cell_type_safe
    cat(sprintf("\n--- CollecTRI: %s ---\n", ct))

    stat_mat <- matrix(rk$primary, ncol = 1, dimnames = list(names(rk$primary), ct_safe))

    ct_res <- tryCatch({
      if (use_network_ct) {
        decoupleR::run_ulm(mat = stat_mat, network = collectri_net,
                           .source = "source", .target = "target", .mor = "mor", minsize = 5)
      } else {
        decoupleR::run_ulm(mat = stat_mat, net = collectri_net,
                           .source = "source", .target = "target", .mor = "mor", minsize = 5)
      }
    }, error = function(e) {
      cat(sprintf("[WARNING] CollecTRI failed: %s\n", e$message))
      data.frame()
    })

    if (nrow(ct_res) > 0) {
      ct_res <- as.data.frame(ct_res)
      ct_res$analysis_name <- "primary_anti_PD1"
      ct_res$cell_type <- ct
      ct_res$n_targets_tested <- NA_integer_
      ct_res$positive_score_direction <- "More_active_in_Responder"
      ct_res$negative_score_direction <- "More_active_in_Non_responder"
      ct_res$FDR_within_celltype <- p.adjust(ct_res$p_value, method = "BH")

      out_path <- sprintf("03_results/step07_programs/collectri/%s__CollecTRI_TF_activity.csv", ct_safe)
      write.csv(ct_res, out_path, row.names = FALSE)
      all_collectri[[ct]] <- ct_res

      n_sig <- sum(ct_res$FDR_within_celltype < 0.05, na.rm = TRUE)
      cat(sprintf("[OK] %s: %d TFs, %d sig\n", ct, nrow(ct_res), n_sig))
    }
  }
} else {
  cat("[INFO] CollecTRI not available, skipping\n")
}

## =========================================================================
## XIV. Multiple Testing Correction (Global FDR)
## =========================================================================
cat("\n[XIV] Computing global FDR across cell types...\n")

## A. Hallmark global FDR
if (length(all_hallmark_primary) > 0) {
  hm_all <- do.call(rbind, all_hallmark_primary)
  hm_all$FDR_global_Hallmark <- p.adjust(hm_all$pval, method = "BH")
  con <- gzfile("03_results/step07_programs/combined/GSE243013_Hallmark_all_celltypes.csv.gz", "w")
  write.csv(hm_all, con, row.names = FALSE)
  close(con)
  cat(sprintf("[INFO] Hallmark combined: %d rows, %d global sig\n",
              nrow(hm_all), sum(hm_all$FDR_global_Hallmark < 0.05, na.rm = TRUE)))
} else {
  hm_all <- data.frame()
  cat("[WARNING] No Hallmark results to combine\n")
}

## B. Reactome global FDR
if (length(all_reactome_primary) > 0) {
  rx_all <- do.call(rbind, all_reactome_primary)
  rx_all$FDR_global_Reactome <- p.adjust(rx_all$pval, method = "BH")
  con <- gzfile("03_results/step07_programs/combined/GSE243013_Reactome_all_celltypes.csv.gz", "w")
  write.csv(rx_all, con, row.names = FALSE)
  close(con)
  cat(sprintf("[INFO] Reactome combined: %d rows, %d global sig\n",
              nrow(rx_all), sum(rx_all$FDR_global_Reactome < 0.05, na.rm = TRUE)))
} else {
  rx_all <- data.frame()
  cat("[WARNING] No Reactome results to combine\n")
}

## C. PROGENy global FDR
if (length(all_progeny) > 0) {
  prog_all <- do.call(rbind, all_progeny)
  prog_all$FDR_global_PROGENy <- p.adjust(prog_all$p_value, method = "BH")
  write.csv(prog_all, "03_results/step07_programs/combined/GSE243013_PROGENy_all_celltypes.csv",
            row.names = FALSE)
  cat(sprintf("[INFO] PROGENy combined: %d rows, %d global sig\n",
              nrow(prog_all), sum(prog_all$FDR_global_PROGENy < 0.05, na.rm = TRUE)))
} else {
  prog_all <- data.frame()
}

## D. CollecTRI global FDR
if (length(all_collectri) > 0) {
  ct_all <- do.call(rbind, all_collectri)
  ct_all$FDR_global_CollecTRI <- p.adjust(ct_all$p_value, method = "BH")
  con <- gzfile("03_results/step07_programs/combined/GSE243013_CollecTRI_all_celltypes.csv.gz", "w")
  write.csv(ct_all, con, row.names = FALSE)
  close(con)
  cat(sprintf("[INFO] CollecTRI combined: %d rows, %d global sig\n",
              nrow(ct_all), sum(ct_all$FDR_global_CollecTRI < 0.05, na.rm = TRUE)))
} else {
  ct_all <- data.frame()
}

## =========================================================================
## XVI. Primary vs Sensitivity Rank Concordance
## =========================================================================
cat("\n[XVI] Primary vs sensitivity rank concordance...\n")

if (length(all_hallmark_primary) > 0 && length(all_hallmark_sensitivity) > 0) {
  concordance_rows <- list()
  for (ct in names(all_hallmark_primary)) {
    if (!(ct %in% names(all_hallmark_sensitivity))) next
    prim <- all_hallmark_primary[[ct]]
    sens <- all_hallmark_sensitivity[[ct]]
    merged <- merge(
      data.frame(pathway = prim$pathway, NES_primary = prim$NES, padj_primary = prim$padj,
                 stringsAsFactors = FALSE),
      data.frame(pathway = sens$pathway, NES_sensitivity = sens$NES, padj_sensitivity = sens$padj,
                 stringsAsFactors = FALSE),
      by = "pathway", all = TRUE
    )
    merged$cell_type <- ct
    merged$direction_consistent <- sign(merged$NES_primary) == sign(merged$NES_sensitivity)
    merged$both_significant <- merged$padj_primary < 0.05 & merged$padj_sensitivity < 0.05
    merged$primary_sig_same_dir <- merged$padj_primary < 0.05 & merged$direction_consistent
    merged$rank_method_status <- ifelse(merged$direction_consistent, "CONCORDANT",
                                        ifelse(is.na(merged$direction_consistent), "NA", "RANK_METHOD_DIRECTION_CONFLICT"))
    concordance_rows[[length(concordance_rows) + 1]] <- merged
  }
  if (length(concordance_rows) > 0) {
    rank_conc <- do.call(rbind, concordance_rows)
    write.csv(rank_conc, "03_results/step07_programs/combined/GSE243013_Hallmark_rank_method_concordance.csv",
              row.names = FALSE)
    n_conflict <- sum(rank_conc$rank_method_status == "RANK_METHOD_DIRECTION_CONFLICT", na.rm = TRUE)
    cat(sprintf("[INFO] Rank concordance: %d conflicts out of %d\n", n_conflict, nrow(rank_conc)))
  }
}

## =========================================================================
## XVII. Primary vs Strict Concordance
## =========================================================================
cat("\n[XVII] Primary vs strict concordance...\n")

for (coll_name in c("Hallmark", "Reactome")) {
  prim_list <- if (coll_name == "Hallmark") all_hallmark_primary else all_reactome_primary
  strict_list <- if (coll_name == "Hallmark") all_hallmark_strict else all_reactome_strict

  if (length(prim_list) > 0 && length(strict_list) > 0) {
    conc_rows <- list()
    for (ct in names(prim_list)) {
      if (!(ct %in% names(strict_list))) next
      prim <- prim_list[[ct]]
      strt <- strict_list[[ct]]
      merged <- merge(
        data.frame(pathway = prim$pathway, primary_NES = prim$NES, primary_padj = prim$padj,
                   stringsAsFactors = FALSE),
        data.frame(pathway = strt$pathway, strict_NES = strt$NES, strict_padj = strt$padj,
                   stringsAsFactors = FALSE),
        by = "pathway", all = TRUE
      )
      merged$cell_type <- ct
      merged$direction_concordant <- sign(merged$primary_NES) == sign(merged$strict_NES)
      merged$primary_significant <- merged$primary_padj < 0.05
      merged$strict_significant <- merged$strict_padj < 0.05
      merged$both_significant <- merged$primary_significant & merged$strict_significant
      merged$primary_significant_strict_same_dir <- merged$primary_significant & merged$direction_concordant
      conc_rows[[length(conc_rows) + 1]] <- merged
    }
    if (length(conc_rows) > 0) {
      conc_df <- do.call(rbind, conc_rows)
      out_path <- sprintf("03_results/step07_programs/combined/GSE243013_primary_vs_strict_%s_concordance.csv", coll_name)
      write.csv(conc_df, out_path, row.names = FALSE)
      n_concordant <- sum(conc_df$direction_concordant, na.rm = TRUE)
      cat(sprintf("[INFO] %s primary vs strict: %d concordant out of %d\n",
                  coll_name, n_concordant, nrow(conc_df)))
    }
  }
}

## =========================================================================
## XVIII. Cross-Cell-Type Recurrence
## =========================================================================
cat("\n[XVIII] Cross-cell-type recurrence...\n")

compute_pathway_recurrence <- function(all_results, fdr_col, name_label) {
  if (length(all_results) == 0) return(data.frame())
  combined <- do.call(rbind, all_results)
  pathways <- unique(combined$pathway)
  rec_rows <- list()
  for (pw in pathways) {
    sub <- combined[combined$pathway == pw, ]
    n_tested <- nrow(sub)
    n_sig <- sum(sub[[fdr_col]] < 0.05, na.rm = TRUE)
    n_resp <- sum(sub[[fdr_col]] < 0.05 & sub$NES > 0, na.rm = TRUE)
    n_nonresp <- sum(sub[[fdr_col]] < 0.05 & sub$NES < 0, na.rm = TRUE)
    direction_consistent <- if (n_sig >= 2) {
      n_resp == 0 || n_nonresp == 0
    } else NA

    ## Leading edge genes
    sig_sub <- sub[sub[[fdr_col]] < 0.05 & !is.na(sub[[fdr_col]]), ]
    all_le <- unlist(strsplit(paste(sig_sub$leadingEdge, collapse = ";"), ";"))
    all_le <- all_le[all_le != ""]
    le_union <- unique(all_le)
    le_intersect <- if (nrow(sig_sub) >= 2) {
      le_lists <- strsplit(sig_sub$leadingEdge, ";")
      Reduce(intersect, le_lists)
    } else if (nrow(sig_sub) == 1) {
      unlist(strsplit(sig_sub$leadingEdge, ";"))
    } else character(0)

    rec_rows[[length(rec_rows) + 1]] <- data.frame(
      pathway = pw, n_celltypes_tested = n_tested, n_significant = n_sig,
      n_Responder_direction = n_resp, n_Non_responder_direction = n_nonresp,
      direction_consistent_pct = ifelse(n_tested > 0, n_sig / n_tested, 0),
      direction_consistent = direction_consistent,
      min_padj = min(sub[[fdr_col]], na.rm = TRUE),
      median_NES = median(sub$NES, na.rm = TRUE),
      max_abs_NES = max(abs(sub$NES), na.rm = TRUE),
      significant_celltypes = paste(sub$cell_type[sub[[fdr_col]] < 0.05], collapse = "; "),
      leading_edge_union = paste(unique(le_union), collapse = "; "),
      leading_edge_intersect = paste(unique(le_intersect), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rec_rows)
}

## Hallmark recurrence
hm_rec <- compute_pathway_recurrence(all_hallmark_primary, "FDR_global_Hallmark", "Hallmark")

## Reactome recurrence
rx_rec <- compute_pathway_recurrence(all_reactome_primary, "FDR_global_Reactome", "Reactome")

## TF recurrence
compute_tf_recurrence <- function(all_ct_results, fdr_col) {
  if (length(all_ct_results) == 0) return(data.frame())
  combined <- do.call(rbind, all_ct_results)
  tfs <- unique(combined$source)
  rec_rows <- list()
  for (tf in tfs) {
    sub <- combined[combined$source == tf, ]
    n_tested <- nrow(sub)
    n_sig <- sum(sub[[fdr_col]] < 0.05, na.rm = TRUE)
    n_resp <- sum(sub[[fdr_col]] < 0.05 & sub$score > 0, na.rm = TRUE)
    n_nonresp <- sum(sub[[fdr_col]] < 0.05 & sub$score < 0, na.rm = TRUE)
    direction_consistent <- if (n_sig >= 2) n_resp == 0 || n_nonresp == 0 else NA

    rec_rows[[length(rec_rows) + 1]] <- data.frame(
      TF = tf, n_celltypes_tested = n_tested, n_significant = n_sig,
      n_Responder_active = n_resp, n_Non_responder_active = n_nonresp,
      direction_consistent_pct = ifelse(n_tested > 0, n_sig / n_tested, 0),
      direction_consistent = direction_consistent,
      min_FDR = min(sub[[fdr_col]], na.rm = TRUE),
      median_score = median(sub$score, na.rm = TRUE),
      max_abs_score = max(abs(sub$score), na.rm = TRUE),
      significant_celltypes = paste(sub$cell_type[sub[[fdr_col]] < 0.05], collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rec_rows)
}

tf_rec <- compute_tf_recurrence(all_collectri, "FDR_global_CollecTRI")

## Save recurrence tables
if (nrow(hm_rec) > 0) {
  write.csv(hm_rec, "03_results/step07_programs/combined/GSE243013_Hallmark_recurrence.csv", row.names = FALSE)
}
if (nrow(rx_rec) > 0) {
  write.csv(rx_rec, "03_results/step07_programs/combined/GSE243013_Reactome_recurrence.csv", row.names = FALSE)
}

## Combined recurrence
all_rec <- rbind(
  if (nrow(hm_rec) > 0) data.frame(collection = "Hallmark", hm_rec, stringsAsFactors = FALSE),
  if (nrow(rx_rec) > 0) data.frame(collection = "Reactome", rx_rec, stringsAsFactors = FALSE)
)
if (nrow(all_rec) > 0) {
  write.csv(all_rec, "03_results/step07_programs/combined/GSE243013_pathway_recurrence_across_celltypes.csv",
            row.names = FALSE)
}

if (nrow(tf_rec) > 0) {
  write.csv(tf_rec, "03_results/step07_programs/combined/GSE243013_TF_recurrence_across_celltypes.csv",
            row.names = FALSE)
}

## =========================================================================
## XIX. Cell-Type Specific Programs
## =========================================================================
cat("\n[XIX] Cell-type specific programs...\n")

## Pathway-specific
ct_specific_pathways <- list()
if (nrow(all_rec) > 0) {
  for (i in seq_len(nrow(all_rec))) {
    row <- all_rec[i, ]
    if (row$n_significant == 1 && !is.na(row$n_significant)) {
      ct_name <- unlist(strsplit(row$significant_celltypes, "; "))[1]
      global_fdr <- if (row$collection == "Hallmark") {
        sub <- hm_all[hm_all$pathway == row$pathway & hm_all$cell_type == ct_name, ]
        if (nrow(sub) > 0) sub$FDR_global_Hallmark[1] else NA
      } else {
        sub <- rx_all[rx_all$pathway == row$pathway & rx_all$cell_type == ct_name, ]
        if (nrow(sub) > 0) sub$FDR_global_Reactome[1] else NA
      }
      if (!is.na(global_fdr) && global_fdr < 0.10) {
        ct_specific_pathways[[length(ct_specific_pathways) + 1]] <- data.frame(
          cell_type = ct_name, collection = row$collection, pathway = row$pathway,
          NES = row$median_NES, global_FDR = global_fdr, stringsAsFactors = FALSE
        )
      }
    }
  }
}
if (length(ct_specific_pathways) > 0) {
  ct_spec_df <- do.call(rbind, ct_specific_pathways)
  write.csv(ct_spec_df, "03_results/step07_programs/combined/GSE243013_celltype_specific_pathways.csv",
            row.names = FALSE)
  cat(sprintf("[INFO] Cell-type specific pathways: %d\n", nrow(ct_spec_df)))
}

## TF-specific
ct_specific_tfs <- list()
if (nrow(tf_rec) > 0) {
  for (i in seq_len(nrow(tf_rec))) {
    row <- tf_rec[i, ]
    if (row$n_significant == 1 && !is.na(row$n_significant)) {
      ct_name <- unlist(strsplit(row$significant_celltypes, "; "))[1]
      sub <- ct_all[ct_all$source == row$TF & ct_all$cell_type == ct_name, ]
      global_fdr <- if (nrow(sub) > 0) sub$FDR_global_CollecTRI[1] else NA
      if (!is.na(global_fdr) && global_fdr < 0.10) {
        ct_specific_tfs[[length(ct_specific_tfs) + 1]] <- data.frame(
          cell_type = ct_name, TF = row$TF, score = row$median_score,
          global_FDR = global_fdr, stringsAsFactors = FALSE
        )
      }
    }
  }
}
if (length(ct_specific_tfs) > 0) {
  ct_spec_tf_df <- do.call(rbind, ct_specific_tfs)
  write.csv(ct_spec_tf_df, "03_results/step07_programs/combined/GSE243013_celltype_specific_TFs.csv",
            row.names = FALSE)
  cat(sprintf("[INFO] Cell-type specific TFs: %d\n", nrow(ct_spec_tf_df)))
}

## =========================================================================
## XX. TF-Pathway-Leading-Edge Network
## =========================================================================
cat("\n[XX] Building TF-pathway-leading-edge network...\n")

tf_pw_links <- list()

if (nrow(hm_all) > 0 && nrow(ct_all) > 0 && collectri_available) {
  ## Get significant pathways (padj < 0.05)
  sig_hm <- hm_all[hm_all$padj < 0.05 & !is.na(hm_all$padj), ]
  ## Get significant TFs
  sig_tfs <- ct_all[ct_all$FDR_within_celltype < 0.05 & !is.na(ct_all$FDR_within_celltype), ]

  for (i in seq_len(nrow(sig_hm))) {
    hm_row <- sig_hm[i, ]
    ct <- hm_row$cell_type
    le_genes <- unlist(strsplit(hm_row$leadingEdge, ";"))
    le_genes <- le_genes[le_genes != ""]

    ct_tfs <- sig_tfs[sig_tfs$cell_type == ct, ]
    for (j in seq_len(nrow(ct_tfs))) {
      tf_row <- ct_tfs[j, ]
      ## Get TF targets from collectri
      tf_targets <- collectri_net$target[collectri_net$source == tf_row$source]
      overlap <- intersect(le_genes, tf_targets)
      if (length(overlap) >= 3) {
        direction_concordant <- sign(tf_row$score) == sign(hm_row$NES)
        tf_pw_links[[length(tf_pw_links) + 1]] <- data.frame(
          cell_type = ct, TF = tf_row$source, TF_score = tf_row$score,
          TF_FDR = tf_row$FDR_within_celltype,
          pathway = hm_row$pathway, collection = "Hallmark",
          NES = hm_row$NES, pathway_padj = hm_row$padj,
          pathway_global_FDR = hm_row$FDR_global_Hallmark,
          overlap_gene_count = length(overlap),
          overlap_genes = paste(overlap, collapse = "; "),
          direction_concordant = direction_concordant,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(tf_pw_links) > 0) {
  tf_pw_df <- do.call(rbind, tf_pw_links)
  tf_pw_df <- tf_pw_df[order(-tf_pw_df$overlap_gene_count), ]
  con <- gzfile("03_results/step07_programs/combined/GSE243013_TF_pathway_leading_edge_links.csv.gz", "w")
  write.csv(tf_pw_df, con, row.names = FALSE)
  close(con)
  cat(sprintf("[INFO] TF-pathway links: %d associations\n", nrow(tf_pw_df)))
} else {
  tf_pw_df <- data.frame()
  cat("[INFO] No TF-pathway links found\n")
}

## =========================================================================
## XXI. Prioritized Immune Programs
## =========================================================================
cat("\n[XXI] Building prioritized immune programs...\n")

priority_rows <- list()

if (nrow(hm_all) > 0) {
  for (i in seq_len(nrow(hm_all))) {
    row <- hm_all[i, ]
    ct <- row$cell_type
    pw <- row$pathway

    ## Check strict concordance
    strict_dir <- TRUE
    if (length(all_hallmark_strict) > 0 && ct %in% names(all_hallmark_strict)) {
      strt <- all_hallmark_strict[[ct]]
      strt_pw <- strt[strt$pathway == pw, ]
      if (nrow(strt_pw) > 0) {
        strict_dir <- sign(row$NES) == sign(strt_pw$NES)
      }
    }

    ## Check rank method concordance
    rank_dir <- TRUE
    if (length(all_hallmark_sensitivity) > 0 && ct %in% names(all_hallmark_sensitivity)) {
      sens <- all_hallmark_sensitivity[[ct]]
      sens_pw <- sens[sens$pathway == pw, ]
      if (nrow(sens_pw) > 0) {
        rank_dir <- sign(row$NES) == sign(sens_pw$NES)
      }
    }

    ## Check TF support
    sig_tf_count <- 0
    tf_overlap_count <- 0
    top_tfs <- ""
    if (nrow(tf_pw_df) > 0) {
      tf_links <- tf_pw_df[tf_pw_df$cell_type == ct & tf_pw_df$pathway == pw &
                             tf_pw_df$direction_concordant, ]
      sig_tf_count <- nrow(tf_links)
      if (nrow(tf_links) > 0) {
        tf_overlap_count <- max(tf_links$overlap_gene_count)
        top_tfs <- paste(tf_links$TF[1:min(5, nrow(tf_links))], collapse = "; ")
      }
    }

    ## Tier assignment
    is_sig <- row$padj < 0.05 & !is.na(row$padj)
    is_global_sig <- !is.na(row$FDR_global_Hallmark) && row$FDR_global_Hallmark < 0.05
    direction <- if (row$NES > 0) "Responder" else "Non_responder"

    tier <- "Tier 3"
    evidence <- "exploratory"

    if (is_sig && strict_dir && rank_dir && sig_tf_count > 0 && tf_overlap_count >= 3) {
      if (is_global_sig) {
        tier <- "Tier 1"
        evidence <- "strong"
      } else {
        tier <- "Tier 2"
        evidence <- "moderate"
      }
    } else if (is_sig && strict_dir && rank_dir) {
      tier <- "Tier 2"
      evidence <- "moderate"
    } else if (is_sig) {
      tier <- "Tier 3"
      evidence <- "exploratory"
    }

    priority_rows[[length(priority_rows) + 1]] <- data.frame(
      priority_tier = tier, cell_type = ct, collection = "Hallmark",
      pathway = pw, NES = row$NES, padj = row$padj,
      global_FDR = row$FDR_global_Hallmark, direction = direction,
      primary_strict_direction_concordant = strict_dir,
      rank_method_direction_concordant = rank_dir,
      top_supporting_TFs = top_tfs,
      leading_edge_genes = row$leadingEdge,
      TF_leading_edge_overlap = tf_overlap_count,
      evidence_summary = evidence,
      stringsAsFactors = FALSE
    )
  }
}

if (length(priority_rows) > 0) {
  priority_df <- do.call(rbind, priority_rows)
  priority_df <- priority_df[order(match(priority_df$priority_tier, c("Tier 1", "Tier 2", "Tier 3")),
                                   -abs(priority_df$NES)), ]
  write.csv(priority_df, "03_results/step07_programs/combined/GSE243013_prioritized_immune_programs.csv",
            row.names = FALSE)
  cat(sprintf("[INFO] Priority programs: Tier 1=%d, Tier 2=%d, Tier 3=%d\n",
              sum(priority_df$priority_tier == "Tier 1"),
              sum(priority_df$priority_tier == "Tier 2"),
              sum(priority_df$priority_tier == "Tier 3")))
} else {
  priority_df <- data.frame()
}

## =========================================================================
## XXII-XXIII. Figures
## =========================================================================
cat("\n[XXII-XXIII] Generating figures...\n")
flush.console()

with_pdf <- function(filename, width = 7, height = 6, plot_fun) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file = filename, width = width, height = height, onefile = TRUE)
  opened_device <- grDevices::dev.cur()
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && opened_device %in% open_devices)
      try(grDevices::dev.off(which = opened_device), silent = TRUE)
  }, add = TRUE)
  force(plot_fun())
  invisible(filename)
}

plot_recovery <- list()
log_plot <- function(plot_type, path, status, error_msg = "") {
  dev_before <- grDevices::dev.cur()
  plot_recovery[[length(plot_recovery) + 1]] <<- data.frame(
    plot_type = plot_type, path = path, status = status,
    error_message = error_msg, device_before = dev_before,
    device_after = grDevices::dev.cur(), stringsAsFactors = FALSE
  )
  if (grDevices::dev.cur() != 1L) {
    cat(sprintf("[WARNING] Device leak after %s\n", plot_type))
    graphics.off()
  }
}

fig_dir <- "04_figures/step07_programs/combined"

## 1. Hallmark NES heatmap
tryCatch({
  if (nrow(hm_all) > 0) {
    sig_pathways <- unique(hm_all$pathway[hm_all$padj < 0.05 & !is.na(hm_all$padj)])
    if (length(sig_pathways) > 0) {
      sig_hm <- hm_all[hm_all$pathway %in% sig_pathways, ]
      mat <- tapply(sig_hm$NES, list(sig_hm$pathway, sig_hm$cell_type), identity)
      if (is.list(mat)) mat <- do.call(rbind, mat)
      if (!is.matrix(mat) || nrow(mat) == 0) stop("Empty matrix")
      mat[is.na(mat)] <- 0

      with_pdf(file.path(fig_dir, "Hallmark_NES_across_celltypes.pdf"),
               width = max(10, ncol(mat) * 0.8), height = max(8, nrow(mat) * 0.3),
               plot_fun = function() {
                 pheatmap::pheatmap(mat, main = "Hallmark NES across cell types",
                                    cluster_rows = TRUE, cluster_cols = TRUE,
                                    fontsize_row = 5, fontsize_col = 6,
                                    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
               })
      log_plot("Hallmark_NES_heatmap", file.path(fig_dir, "Hallmark_NES_across_celltypes.pdf"), "OK")
    }
  }
}, error = function(e) {
  log_plot("Hallmark_NES_heatmap", "", "FAILED_PLOT", e$message)
  cat(sprintf("[WARNING] Hallmark NES heatmap failed: %s\n", e$message))
})

## 2. Reactome top NES heatmap
tryCatch({
  if (nrow(rx_all) > 0) {
    rx_all$abs_NES <- abs(rx_all$NES)
    top_rx <- rx_all[order(rx_all$FDR_global_Reactome), ]
    top_pathways <- unique(top_rx$pathway[1:min(40, nrow(top_rx))])
    sig_rx <- rx_all[rx_all$pathway %in% top_pathways, ]
    mat <- tapply(sig_rx$NES, list(sig_rx$pathway, sig_rx$cell_type), identity)
    if (is.list(mat)) mat <- do.call(rbind, mat)
    if (!is.matrix(mat) || nrow(mat) == 0) stop("Empty matrix")
    mat[is.na(mat)] <- 0

    if (nrow(mat) > 0) {
      with_pdf(file.path(fig_dir, "Reactome_top_NES_across_celltypes.pdf"),
               width = max(10, ncol(mat) * 0.8), height = max(8, nrow(mat) * 0.25),
               plot_fun = function() {
                 pheatmap::pheatmap(mat, main = "Reactome Top NES across cell types",
                                    cluster_rows = TRUE, cluster_cols = TRUE,
                                    fontsize_row = 4, fontsize_col = 6,
                                    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
               })
      log_plot("Reactome_NES_heatmap", file.path(fig_dir, "Reactome_top_NES_across_celltypes.pdf"), "OK")
    }
  }
}, error = function(e) {
  log_plot("Reactome_NES_heatmap", "", "FAILED_PLOT", e$message)
  cat(sprintf("[WARNING] Reactome heatmap failed: %s\n", e$message))
})

## 3. PROGENy activity heatmap
tryCatch({
  if (nrow(prog_all) > 0) {
    mat <- tapply(prog_all$score, list(prog_all$source, prog_all$cell_type), identity)
    if (is.list(mat)) mat <- do.call(rbind, mat)
    if (!is.matrix(mat) || nrow(mat) == 0) stop("Empty matrix")
    mat[is.na(mat)] <- 0

    with_pdf(file.path(fig_dir, "PROGENy_activity_across_celltypes.pdf"),
             width = max(10, ncol(mat) * 0.8), height = max(6, nrow(mat) * 0.4),
             plot_fun = function() {
               pheatmap::pheatmap(mat, main = "PROGENy Pathway Activity",
                                  cluster_rows = TRUE, cluster_cols = TRUE,
                                  fontsize_row = 6, fontsize_col = 6,
                                  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
             })
    log_plot("PROGENy_heatmap", file.path(fig_dir, "PROGENy_activity_across_celltypes.pdf"), "OK")
  }
}, error = function(e) {
  log_plot("PROGENy_heatmap", "", "FAILED_PLOT", e$message)
  cat(sprintf("[WARNING] PROGENy heatmap failed: %s\n", e$message))
})

## 4. CollecTRI TF heatmap
tryCatch({
  if (nrow(ct_all) > 0) {
    ct_all$abs_score <- abs(ct_all$score)
    top_tfs <- ct_all[order(ct_all$FDR_global_CollecTRI), ]
    tf_names <- unique(top_tfs$source[1:min(40, nrow(top_tfs))])
    sig_ct <- ct_all[ct_all$source %in% tf_names, ]
    mat <- tapply(sig_ct$score, list(sig_ct$source, sig_ct$cell_type), identity)
    if (is.list(mat)) mat <- do.call(rbind, mat)
    if (!is.matrix(mat) || nrow(mat) == 0) stop("Empty matrix")
    mat[is.na(mat)] <- 0

    if (nrow(mat) > 0) {
      with_pdf(file.path(fig_dir, "CollecTRI_top_TF_activity_across_celltypes.pdf"),
               width = max(10, ncol(mat) * 0.8), height = max(8, nrow(mat) * 0.25),
               plot_fun = function() {
                 pheatmap::pheatmap(mat, main = "CollecTRI TF Activity",
                                    cluster_rows = TRUE, cluster_cols = TRUE,
                                    fontsize_row = 5, fontsize_col = 6,
                                    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
               })
      log_plot("CollecTRI_heatmap", file.path(fig_dir, "CollecTRI_top_TF_activity_across_celltypes.pdf"), "OK")
    }
  }
}, error = function(e) {
  log_plot("CollecTRI_heatmap", "", "FAILED_PLOT", e$message)
  cat(sprintf("[WARNING] CollecTRI heatmap failed: %s\n", e$message))
})

## 5. Pathway recurrence top 20
tryCatch({
  if (nrow(all_rec) > 0) {
    top20 <- all_rec[order(-all_rec$n_significant), ][1:min(20, nrow(all_rec)), ]
    top20$label <- sprintf("%s (%s)", top20$pathway, top20$collection)

    with_pdf(file.path(fig_dir, "Pathway_recurrence_top20.pdf"),
             width = 10, height = 7,
             plot_fun = function() {
               print(ggplot(top20, aes(x = reorder(label, n_significant), y = n_significant, fill = collection)) +
                 geom_bar(stat = "identity") +
                 coord_flip() +
                 theme_minimal() +
                 labs(title = "Top 20 Pathways by Cross-Cell-Type Recurrence",
                      x = "Pathway", y = "Number of Significant Cell Types") +
                 scale_fill_manual(values = c("Hallmark" = "#E41A1C", "Reactome" = "#377EB8")))
             })
    log_plot("Pathway_recurrence", file.path(fig_dir, "Pathway_recurrence_top20.pdf"), "OK")
  }
}, error = function(e) {
  log_plot("Pathway_recurrence", "", "FAILED_PLOT", e$message)
  cat(sprintf("[WARNING] Pathway recurrence plot failed: %s\n", e$message))
})

## 6. TF recurrence top 20
tryCatch({
  if (nrow(tf_rec) > 0) {
    top20tf <- tf_rec[order(-tf_rec$n_significant), ][1:min(20, nrow(tf_rec)), ]

    with_pdf(file.path(fig_dir, "TF_recurrence_top20.pdf"),
             width = 10, height = 7,
             plot_fun = function() {
               print(ggplot(top20tf, aes(x = reorder(TF, n_significant), y = n_significant)) +
                 geom_bar(stat = "identity", fill = "#4DAF4A") +
                 coord_flip() +
                 theme_minimal() +
                 labs(title = "Top 20 TFs by Cross-Cell-Type Recurrence",
                      x = "Transcription Factor", y = "Number of Significant Cell Types"))
             })
    log_plot("TF_recurrence", file.path(fig_dir, "TF_recurrence_top20.pdf"), "OK")
  }
}, error = function(e) {
  log_plot("TF_recurrence", "", "FAILED_PLOT", e$message)
  cat(sprintf("[WARNING] TF recurrence plot failed: %s\n", e$message))
})

## 7. Primary vs Strict NES scatter
for (coll_name in c("Hallmark", "Reactome")) {
  tryCatch({
    conc_path <- sprintf("03_results/step07_programs/combined/GSE243013_primary_vs_strict_%s_concordance.csv", coll_name)
    if (file.exists(conc_path)) {
      conc <- read.csv(conc_path, stringsAsFactors = FALSE)
      conc <- conc[!is.na(conc$primary_NES) & !is.na(conc$strict_NES), ]

      with_pdf(file.path(fig_dir, sprintf("Primary_vs_strict_%s_NES.pdf", coll_name)),
               width = 8, height = 7,
               plot_fun = function() {
                 print(ggplot(conc, aes(x = primary_NES, y = strict_NES)) +
                   geom_point(aes(color = direction_concordant), alpha = 0.5, size = 1.5) +
                   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
                   theme_minimal() +
                   labs(title = sprintf("Primary vs Strict Cohort NES (%s)", coll_name),
                        x = "Primary NES", y = "Strict NES",
                        color = "Direction Concordant") +
                   scale_color_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")))
               })
      log_plot(sprintf("Primary_vs_strict_%s", coll_name),
               file.path(fig_dir, sprintf("Primary_vs_strict_%s_NES.pdf", coll_name)), "OK")
    }
  }, error = function(e) {
    log_plot(sprintf("Primary_vs_strict_%s", coll_name), "", "FAILED_PLOT", e$message)
    cat(sprintf("[WARNING] %s scatter failed: %s\n", coll_name, e$message))
  })
}

## 8. Tier 1 enrichment plots (max 20)
tryCatch({
  if (nrow(priority_df) > 0) {
    tier1 <- priority_df[priority_df$priority_tier == "Tier 1", ]
    n_plots <- min(20, nrow(tier1))
    if (n_plots > 0) {
      for (p_idx in seq_len(n_plots)) {
        t1 <- tier1[p_idx, ]
        ct <- t1$cell_type
        pw <- t1$pathway
        rk <- primary_ranks[[ct]]
        if (is.null(rk)) next

        pathway_genes <- hallmark_gsets[[pw]]
        if (is.null(pathway_genes)) next

        rank_vec <- rk$primary
        genes_in_rank <- pathway_genes[pathway_genes %in% names(rank_vec)]

        if (length(genes_in_rank) < 5) next

        plot_df <- data.frame(
          gene = names(rank_vec), rank = rank_vec,
          in_pathway = names(rank_vec) %in% genes_in_rank,
          stringsAsFactors = FALSE
        )
        plot_df$order <- seq_len(nrow(plot_df))

        ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
        pw_safe <- gsub("[^A-Za-z0-9_]", "_", pw)
        enrichment_path <- sprintf("04_figures/step07_programs/fgsea/%s__%s__enrichment.pdf",
                                   ct_safe, pw_safe)

        with_pdf(enrichment_path, width = 10, height = 5,
                 plot_fun = function() {
                   par(mfrow = c(1, 2))
                   ## Rank plot
                   colors <- ifelse(plot_df$in_pathway, "#E41A1C", "grey70")
                   plot(plot_df$order, plot_df$rank, col = colors, pch = 16, cex = 0.3,
                        main = sprintf("GSEA: %s", pw),
                        xlab = "Gene rank", ylab = "Rank statistic")
                   ## Running enrichment
                   in_idx <- which(plot_df$in_pathway)
                   running_es <- numeric(nrow(plot_df))
                   running_es[in_idx] <- 1 / length(in_idx)
                   running_es[!plot_df$in_pathway] <- -1 / (nrow(plot_df) - length(in_idx))
                   plot_df$cum_es <- cumsum(running_es)
                   plot(plot_df$order, plot_df$cum_es, type = "l", lwd = 2,
                        col = ifelse(t1$NES > 0, "#377EB8", "#E41A1C"),
                        main = sprintf("Running ES (NES=%.2f)", t1$NES),
                        xlab = "Gene rank", ylab = "Enrichment score")
                   abline(h = 0, lty = 2, col = "grey50")
                 })
      }
      log_plot("Tier1_enrichment", "multiple", "OK")
    }
  }
}, error = function(e) {
  log_plot("Tier1_enrichment", "", "FAILED_PLOT", e$message)
  cat(sprintf("[WARNING] Tier 1 enrichment plots failed: %s\n", e$message))
})

## 9. TF-pathway bubble plot
tryCatch({
  if (nrow(tf_pw_df) > 0) {
    top_links <- tf_pw_df[1:min(30, nrow(tf_pw_df)), ]

    with_pdf(file.path(fig_dir, "Top_TF_pathway_links.pdf"),
             width = 12, height = 8,
             plot_fun = function() {
               top_links$label <- sprintf("%s -> %s", top_links$TF, top_links$pathway)
               print(ggplot(top_links, aes(x = TF, y = pathway, size = overlap_gene_count,
                                           color = direction_concordant)) +
                 geom_point(alpha = 0.7) +
                 theme_minimal() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                 labs(title = "Top TF-Pathway Leading-Edge Links",
                      x = "Transcription Factor", y = "Hallmark Pathway",
                      size = "Overlap genes", color = "Direction concordant") +
                 scale_color_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")))
             })
    log_plot("TF_pathway_bubble", file.path(fig_dir, "Top_TF_pathway_links.pdf"), "OK")
  }
}, error = function(e) {
  log_plot("TF_pathway_bubble", "", "FAILED_PLOT", e$message)
  cat(sprintf("[WARNING] TF-pathway bubble plot failed: %s\n", e$message))
})

## Save plot recovery
if (length(plot_recovery) > 0) {
  pr_df <- do.call(rbind, plot_recovery)
  write.csv(pr_df, "03_results/step07_programs/qc/GSE243013_step07_plot_recovery.csv", row.names = FALSE)
  n_ok <- sum(pr_df$status == "OK", na.rm = TRUE)
  n_fail <- sum(pr_df$status == "FAILED_PLOT", na.rm = TRUE)
  cat(sprintf("[INFO] Plots: %d OK, %d FAILED\n", n_ok, n_fail))
}

## =========================================================================
## XXIV. Result Index and Status
## =========================================================================
cat("\n[XXIV] Creating result index...\n")

result_files <- list(
  list(type = "ranks", path = "03_results/step07_programs/ranks", pattern = "*.csv.gz"),
  list(type = "fgsea_primary", path = "03_results/step07_programs/fgsea_primary", pattern = "*.csv.gz"),
  list(type = "fgsea_strict", path = "03_results/step07_programs/fgsea_strict", pattern = "*.csv.gz"),
  list(type = "progeny", path = "03_results/step07_programs/progeny", pattern = "*.csv"),
  list(type = "collectri", path = "03_results/step07_programs/collectri", pattern = "*.csv"),
  list(type = "combined", path = "03_results/step07_programs/combined", pattern = "*.csv*")
)

idx_rows <- list()
for (rf in result_files) {
  files <- list.files(rf$path, pattern = glob2rx(rf$pattern), full.names = TRUE)
  for (f in files) {
    fi <- file.info(f)
    idx_rows[[length(idx_rows) + 1]] <- data.frame(
      result_type = rf$type, file_path = f, file_size = fi$size,
      status = "OK", error_message = "", stringsAsFactors = FALSE
    )
  }
}
if (length(idx_rows) > 0) {
  idx_df <- do.call(rbind, idx_rows)
  write.csv(idx_df, "03_results/step07_programs/combined/GSE243013_step07_result_index.csv", row.names = FALSE)
}

## Model status
model_status_rows <- list()
for (ct in primary_analyzed) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  has_hm <- any(sapply(all_hallmark_primary, function(x) x$cell_type[1] == ct))
  has_rx <- any(sapply(all_reactome_primary, function(x) x$cell_type[1] == ct))
  has_prog <- ct %in% names(all_progeny)
  has_ctri <- ct %in% names(all_collectri)

  model_status_rows[[length(model_status_rows) + 1]] <- data.frame(
    cell_type = ct, cell_type_safe = ct_safe,
    hallmarks_complete = has_hm, reactome_complete = has_rx,
    progeny_complete = has_prog, collectri_complete = has_ctri,
    status = "COMPLETE", stringsAsFactors = FALSE
  )
}
model_status_step07 <- do.call(rbind, model_status_rows)
write.csv(model_status_step07, "03_results/step07_programs/combined/GSE243013_step07_model_status.csv",
          row.names = FALSE)

## =========================================================================
## XXV. Analysis Definition
## =========================================================================
cat("\n[XXV] Writing analysis definition...\n")

def_text <- c(
  "GSE243013 Step 07 Analysis Definition",
  "======================================",
  "",
  "Input: patient-level edgeR differential expression complete gene rankings from Step 06",
  "Not limited to significant DE genes",
  "Primary ranking: sign(logFC) * sqrt(F)",
  "Positive NES or activity score: Responder relatively higher/more active",
  "Negative NES or activity score: Non_responder relatively higher/more active",
  "Hallmark and Reactome used for primary pathway analysis",
  "PROGENy used for signal pathway activity inference",
  "CollecTRI used for transcription factor activity inference",
  "Primary and strict chemoimmunotherapy cohorts used for robustness check",
  "Global FDR used for cross-cell-type strict screening",
  "Leading-edge genes used to identify core enrichment drivers",
  "TF-pathway overlap is regulatory support evidence only",
  "TF-pathway associations not stated as causal",
  "Current results are post-treatment pathological response-associated immune programs",
  "Not referred to as pre-treatment predictive biomarkers",
  "TCGA, CPTAC, and DepMap external validation not yet completed"
)
writeLines(def_text, "00_config/GSE243013_step07_analysis_definition.txt")

## =========================================================================
## XXVI. Completion Marker
## =========================================================================
cat("\n[XXVI] Creating completion marker...\n")

total_runtime <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))

n_hallmark_sig <- if (nrow(hm_all) > 0) sum(hm_all$padj < 0.05 & !is.na(hm_all$padj)) else 0
n_reactome_sig <- if (nrow(rx_all) > 0) sum(rx_all$padj < 0.05 & !is.na(rx_all$padj)) else 0
n_hallmark_global_sig <- if (nrow(hm_all) > 0) sum(hm_all$FDR_global_Hallmark < 0.05 & !is.na(hm_all$FDR_global_Hallmark)) else 0
n_progeny_sig <- if (nrow(prog_all) > 0) sum(prog_all$FDR_within_celltype < 0.05 & !is.na(prog_all$FDR_within_celltype)) else 0
n_collectri_sig <- if (nrow(ct_all) > 0) sum(ct_all$FDR_within_celltype < 0.05 & !is.na(ct_all$FDR_within_celltype)) else 0

n_tier1 <- sum(priority_df$priority_tier == "Tier 1", na.rm = TRUE)
n_tier2 <- sum(priority_df$priority_tier == "Tier 2", na.rm = TRUE)
n_tier3 <- sum(priority_df$priority_tier == "Tier 3", na.rm = TRUE)

n_plot_ok <- sum(sapply(plot_recovery, function(x) x$status == "OK"), na.rm = TRUE)
n_plot_fail <- sum(sapply(plot_recovery, function(x) x$status == "FAILED_PLOT"), na.rm = TRUE)

resource_warnings <- character(0)
if (!progeny_available) resource_warnings <- c(resource_warnings, "PROGENy unavailable")
if (!collectri_available) resource_warnings <- c(resource_warnings, "CollecTRI unavailable")

disk_after <- tryCatch({
  as.numeric(system("df -Pk . | tail -1 | awk '{print $4}'", intern = TRUE)) / (1024 * 1024)
}, error = function(e) NA_real_)

all_checks <- (
  length(primary_analyzed) >= 1 &&
  nrow(hm_all) > 0 &&
  nrow(rx_all) > 0 &&
  file.exists("03_results/step07_programs/combined/GSE243013_Hallmark_all_celltypes.csv.gz") &&
  file.exists("03_results/step07_programs/combined/GSE243013_Reactome_all_celltypes.csv.gz") &&
  file.exists("03_results/step07_programs/combined/GSE243013_pathway_recurrence_across_celltypes.csv") &&
  file.exists("03_results/step07_programs/combined/GSE243013_prioritized_immune_programs.csv")
)

if (all_checks) {
  complete_text <- c(
    "GSE243013 Step 07 COMPLETE",
    "===========================",
    "",
    sprintf("Completion time: %s", Sys.time()),
    sprintf("R version: %s", R.version.string),
    sprintf("fgsea version: %s", packageVersion("fgsea")),
    sprintf("decoupleR version: %s", packageVersion("decoupleR")),
    sprintf("msigdbr version: %s", packageVersion("msigdbr")),
    sprintf("MSigDB version: 2024.1.Hs"),
    sprintf("Primary cell types analyzed: %d", length(primary_analyzed)),
    sprintf("Hallmark gene sets: %d", length(hallmark_gsets)),
    sprintf("Reactome gene sets: %d", length(reactome_gsets)),
    sprintf("Significant Hallmark pathway-celltype: %d", n_hallmark_sig),
    sprintf("Significant Reactome pathway-celltype: %d", n_reactome_sig),
    sprintf("Global FDR significant pathways: %d", n_hallmark_global_sig),
    sprintf("Significant PROGENy pathway-celltype: %d", n_progeny_sig),
    sprintf("Significant CollecTRI TF-celltype: %d", n_collectri_sig),
    sprintf("Tier 1 programs: %d", n_tier1),
    sprintf("Tier 2 programs: %d", n_tier2),
    sprintf("Tier 3 programs: %d", n_tier3),
    sprintf("Resource warnings: %s", ifelse(length(resource_warnings) > 0,
                                            paste(resource_warnings, collapse = "; "), "none")),
    sprintf("Plots OK: %d, Failed: %d", n_plot_ok, n_plot_fail),
    sprintf("Disk available: %.1f GB", disk_after),
    sprintf("Total runtime: %.1f seconds", total_runtime),
    "",
    "No edgeR models re-run",
    "No Step 06 results modified"
  )
  writeLines(complete_text, "03_results/GSE243013_step07_COMPLETE.txt")
  cat("[INFO] Step 07 COMPLETE marker written\n")
} else {
  failed_text <- c(
    "GSE243013 Step 07 FAILED",
    sprintf("Time: %s", Sys.time()),
    sprintf("Primary cell types: %d", length(primary_analyzed)),
    sprintf("Hallmark tests: %d", nrow(hm_all)),
    sprintf("Reactome tests: %d", nrow(rx_all))
  )
  writeLines(failed_text, "03_results/GSE243013_step07_FAILED.txt")
  cat("[WARNING] Step 07 FAILED marker written\n")
}

## Final device check
if (grDevices::dev.cur() != 1L) {
  cat("[WARNING] Closing remaining graphics devices\n")
  graphics.off()
}

cat(sprintf("\nStep 07 status: %s\n", ifelse(all_checks, "COMPLETE", "FAILED")))
cat(sprintf("Runtime: %.1f seconds\n", total_runtime))
cat("========================================================================\n")
cat("Step 07 completed.\n")
cat("========================================================================\n")
