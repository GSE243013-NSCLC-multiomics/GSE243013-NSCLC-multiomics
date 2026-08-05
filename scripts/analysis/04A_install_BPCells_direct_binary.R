## =========================================================================
## Step 04A: Direct Binary Install of BPCells (bypass R-universe index)
## =========================================================================

cat("========================================================================\n")
cat("Step 04A: Direct Binary Install of BPCells\n")
cat("========================================================================\n\n")

RESULTS <- "03_results"
CACHE   <- "00_config/package_cache"
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

step_start <- Sys.time()

## =========================================================================
## III. System Check
## =========================================================================
cat("[III] System check...\n")

preflight <- c(
  sprintf("Date: %s", Sys.time()),
  sprintf("R version: %s", R.version.string),
  sprintf("R major.minor: %s.%s", R.version$major, R.version$minor),
  sprintf("R platform: %s", R.version$platform),
  sprintf("lib Paths: %s", paste(.libPaths(), collapse = ", "))
)

tryCatch({
  os <- Sys.info()
  preflight <- c(preflight, sprintf("OS: %s %s", os["sysname"], os["release"]))
  preflight <- c(preflight, sprintf("Machine: %s", os["machine"]))
}, error = function(e) {
  preflight <- c(preflight, sprintf("OS error: %s", conditionMessage(e)))
})

r_major <- as.integer(R.version$major)
r_minor <- as.integer(strsplit(R.version$minor, "\\.")[[1]][1])
is_r46 <- r_major == 4 && r_minor >= 6
is_arm64 <- grepl("aarch64|arm64", R.version$platform, ignore.case = TRUE)
is_mac <- grepl("darwin", R.version$platform, ignore.case = TRUE) ||
          Sys.info()["sysname"] == "Darwin"

preflight <- c(preflight, sprintf("R >= 4.6: %s", is_r46))
preflight <- c(preflight, sprintf("ARM64: %s", is_arm64))
preflight <- c(preflight, sprintf("macOS: %s", is_mac))

cat(paste(preflight, collapse = "\n"), "\n\n")

if (!is_r46 || !is_arm64 || !is_mac) {
  cat("[FATAL] Requirements not met. Stopping.\n")
  writeLines(preflight, file.path(RESULTS, "GSE243013_BPCells_binary_preflight.txt"))
  stop("System requirements not met")
}

writeLines(preflight, file.path(RESULTS, "GSE243013_BPCells_binary_preflight.txt"))
cat("[INFO] Saved: GSE243013_BPCells_binary_preflight.txt\n")

## =========================================================================
## IV. Confirm curl
## =========================================================================
cat("\n[IV] Checking curl...\n")

if (!requireNamespace("curl", quietly = TRUE)) {
  cat("[INFO] Installing curl from CRAN...\n")
  install.packages("curl", repos = "https://cloud.r-project.org", type = "binary")
}

cat(sprintf("[INFO] curl version: %s\n", packageVersion("curl")))

## =========================================================================
## V. Install BPCells Dependencies
## =========================================================================
cat("\n[V] Installing BPCells dependencies...\n")

required_deps <- c(
  "Rcpp", "RcppEigen", "magrittr", "Matrix", "rlang", "vctrs",
  "lifecycle", "stringr", "tibble", "dplyr", "tidyr", "readr",
  "ggplot2", "scales", "patchwork", "scattermore", "ggrepel",
  "RColorBrewer", "hexbin"
)

dep_status <- data.frame(
  package = required_deps,
  installed_before = logical(length(required_deps)),
  installation_attempted = logical(length(required_deps)),
  installed_after = logical(length(required_deps)),
  version = character(length(required_deps)),
  status = character(length(required_deps)),
  error_message = character(length(required_deps)),
  stringsAsFactors = FALSE
)

missing_pkgs <- character()
for (i in seq_along(required_deps)) {
  pkg <- required_deps[i]
  has_it <- requireNamespace(pkg, quietly = TRUE)
  dep_status$installed_before[i] <- has_it
  dep_status$version[i] <- if (has_it) as.character(packageVersion(pkg)) else ""
  if (!has_it) missing_pkgs <- c(missing_pkgs, pkg)
}

cat(sprintf("[INFO] Already installed: %d/%d\n", sum(dep_status$installed_before), length(required_deps)))

if (length(missing_pkgs) > 0) {
  cat(sprintf("[INFO] Installing %d missing packages...\n", length(missing_pkgs)))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", type = "binary")
}

for (i in seq_along(required_deps)) {
  pkg <- required_deps[i]
  if (!dep_status$installed_before[i]) {
    dep_status$installation_attempted[i] <- TRUE
    dep_status$installed_after[i] <- requireNamespace(pkg, quietly = TRUE)
    dep_status$version[i] <- if (dep_status$installed_after[i]) as.character(packageVersion(pkg)) else ""
    dep_status$status[i] <- if (dep_status$installed_after[i]) "installed" else "FAILED"
  } else {
    dep_status$installed_after[i] <- TRUE
    dep_status$status[i] <- "already_installed"
  }
}

all_deps_ok <- all(dep_status$installed_after)
cat(sprintf("[INFO] All dependencies installed: %s\n", all_deps_ok))
write.csv(dep_status, file.path(RESULTS, "GSE243013_BPCells_dependency_status.csv"), row.names = FALSE)

## =========================================================================
## VI. Check if BPCells Already Installed
## =========================================================================
cat("\n[VI] Checking if BPCells is already installed...\n")

if (requireNamespace("BPCells", quietly = TRUE)) {
  cat(sprintf("[INFO] BPCells already installed: %s\n", packageVersion("BPCells")))
  cat(sprintf("[INFO] BPCells path: %s\n", find.package("BPCells")))

  ## Test functions
  cat("[INFO] Testing key functions...\n")
  all_funcs <- all(sapply(c("import_matrix_market", "open_matrix_dir", "storage_order"),
                          function(fn) exists(fn, where = asNamespace("BPCells"), mode = "function")))
  cat(sprintf("[INFO] All key functions exist: %s\n", all_funcs))

  if (all(all_deps_ok, all_funcs)) {
    cat("\n[INFO] BPCells is ready. Proceeding to Step 04 rerun.\n")
    install_ok <- TRUE
  } else {
    install_ok <- FALSE
  }
} else {
  cat("[INFO] BPCells not installed. Attempting direct binary download...\n")
  install_ok <- FALSE
}

## =========================================================================
## VII. Direct Binary Download (if needed)
## =========================================================================
if (!install_ok) {
  cat("\n[VII] Attempting direct binary download...\n")

  tgz_path <- file.path(CACHE, "BPCells_0.3.1_R4.6_sonoma_arm64.tgz")

  urls <- c(
    "https://bnprks.r-universe.dev/bin/macosx/sonoma-arm64/contrib/4.6/BPCells_0.3.1.tgz",
    "https://r2.ropensci.org/2998333f5f14dbbd0fa6b73cd8d1252a05c243362f24349636b392411fe189ca"
  )

  download_ok <- FALSE
  used_url <- NA_character_

  for (url in urls) {
    cat(sprintf("[INFO] Trying: %s\n", url))
    tryCatch({
      h <- curl::new_handle()
      curl::handle_setopt(h,
        followlocation = TRUE,
        failonerror = TRUE,
        connecttimeout = 60L,
        timeout = 600L
      )
      curl::curl_download(url = url, destfile = tgz_path, handle = h)

      fsize <- file.info(tgz_path)$size
      if (!is.na(fsize) && fsize > 1024 * 1024) {
        con <- file(tgz_path, "rb")
        magic <- readBin(con, "raw", n = 2)
        close(con)
        if (magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)) {
          download_ok <- TRUE
          used_url <- url
          cat(sprintf("[INFO] Downloaded: %d bytes\n", fsize))
          break
        }
      }
      cat("[WARNING] Downloaded file invalid, trying next URL\n")
    }, error = function(e) {
      cat(sprintf("[ERROR] Download failed: %s\n", conditionMessage(e)))
    })
  }

  if (!download_ok) {
    cat("\n[FATAL] Cannot download BPCells binary automatically.\n")
    cat("[INFO] Cloudflare protection is blocking automated downloads from r-universe.\n\n")
    cat("========================================\n")
    cat("MANUAL INSTALLATION REQUIRED\n")
    cat("========================================\n\n")
    cat("Please run the following commands in R or RStudio:\n\n")
    cat('  install.packages("BPCells",\n')
    cat('    repos = "https://bnprks.r-universe.dev",\n')
    cat('    type = "binary"\n')
    cat('  )\n\n')
    cat("If that fails, try from GitHub source:\n\n")
    cat('  install.packages("remotes")\n')
    cat('  remotes::install_github("bnprks/BPCells")\n\n')
    cat("After installation, re-run Step 04:\n\n")
    cat('  /usr/local/bin/Rscript 01_scripts/04_download_and_import_GSE243013_counts.R\n\n')
    cat("========================================\n")

    writeLines(c(
      "BPCells binary download failed",
      sprintf("Time: %s", Sys.time()),
      "Reason: Cloudflare protection blocks automated downloads from r-universe",
      "",
      "Manual installation required:",
      '  install.packages("BPCells", repos = "https://bnprks.r-universe.dev", type = "binary")',
      "",
      "After manual installation, re-run Step 04."
    ), file.path(RESULTS, "GSE243013_BPCells_INSTALL_MANUAL_REQUIRED.txt"))
    stop("Manual installation required")
  }

  ## Verify and install
  cat("\n[INFO] Verifying downloaded binary...\n")
  tryCatch({
    install.packages(normalizePath(tgz_path), repos = NULL, type = "binary")
    if (requireNamespace("BPCells", quietly = TRUE)) {
      cat(sprintf("[INFO] BPCells installed: %s\n", packageVersion("BPCells")))
      install_ok <- TRUE
    }
  }, error = function(e) {
    cat("[ERROR] Installation failed:", conditionMessage(e), "\n")
  })
}

## =========================================================================
## VIII. Final Verification
## =========================================================================
if (install_ok) {
  cat("\n[VIII] Final verification...\n")

  cat(sprintf("[INFO] BPCells version: %s\n", packageVersion("BPCells")))
  cat(sprintf("[INFO] BPCells path: %s\n", find.package("BPCells")))

  tryCatch({
    library(BPCells)
    cat("[INFO] BPCells loaded successfully\n")
  }, error = function(e) {
    cat("[ERROR] library(BPCells) failed:", conditionMessage(e), "\n")
  })

  ## Small test
  cat("[INFO] Running small matrix test...\n")
  test_dir <- "02_data/tmp/BPCells_install_test"
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)

  test_passed <- FALSE
  tryCatch({
    small_mat <- Matrix::sparseMatrix(
      i = c(1, 1, 2, 3, 3), j = c(1, 3, 2, 1, 4),
      x = c(5, 3, 8, 2, 7), dims = c(3, 4),
      dimnames = list(c("gene_A", "gene_B", "gene_C"),
                      c("cell_1", "cell_2", "cell_3", "cell_4"))
    )

    bp_test_dir <- file.path(test_dir, "small_test")
    dir.create(bp_test_dir, showWarnings = FALSE)
    BPCells::write_matrix_dir(small_mat, dir = bp_test_dir)
    bp_reopen <- BPCells::open_matrix_dir(bp_test_dir)

    test_passed <- all(
      nrow(bp_reopen) == 3, ncol(bp_reopen) == 4,
      identical(rownames(bp_reopen), rownames(small_mat)),
      identical(colnames(bp_reopen), colnames(small_mat)),
      identical(as.matrix(bp_reopen), as.matrix(small_mat))
    )
    cat(sprintf("[INFO] Small matrix test: %s\n", ifelse(test_passed, "PASSED", "FAILED")))
  }, error = function(e) {
    cat("[ERROR] Test failed:", conditionMessage(e), "\n")
  })

  if (test_passed) {
    cat("\n[INFO] Creating installation complete marker...\n")
    writeLines(c(
      "GSE243013 BPCells INSTALL COMPLETE",
      sprintf("Time: %s", Sys.time()),
      sprintf("R: %s", R.version.string),
      sprintf("BPCells: %s", packageVersion("BPCells")),
      sprintf("Path: %s", find.package("BPCells")),
      "Small matrix test: PASSED"
    ), file.path(RESULTS, "GSE243013_BPCells_INSTALL_COMPLETE.txt"))
    cat("[INFO] Saved: GSE243013_BPCells_INSTALL_COMPLETE.txt\n")
  }
}

## =========================================================================
## IX. Re-run Step 04
## =========================================================================
cat("\n[IX] Re-running Step 04...\n")

if (install_ok) {
  step04_script <- "01_scripts/04_download_and_import_GSE243013_counts.R"
  step04_log <- "logs/04_download_and_import_GSE243013_counts_rerun.log"

  if (file.exists(step04_script)) {
    cat(sprintf("[INFO] Running: %s\n", step04_script))
    exit_code <- system2("/usr/local/bin/Rscript", step04_script,
                         stdout = step04_log, stderr = step04_log)
    cat(sprintf("[INFO] Step 04 rerun exit code: %d\n", exit_code))
  }
}

## =========================================================================
## Summary
## =========================================================================
total_time <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
cat("\n========================================================================\n")
cat("STEP 04A SUMMARY\n")
cat("========================================================================\n")
cat(sprintf("BPCells installed: %s\n", install_ok))
if (install_ok) {
  cat(sprintf("BPCells version: %s\n", packageVersion("BPCells")))
  cat(sprintf("BPCells path: %s\n", find.package("BPCells")))
}
cat(sprintf("Total runtime: %.1f seconds\n", total_time))
cat(sprintf("Status: %s\n", ifelse(install_ok, "COMPLETE", "MANUAL_INSTALL_REQUIRED")))
cat("========================================================================\n")
