## =========================================================================
## Step 00: R Environment Check & GSE243013 File Manifest
## =========================================================================

options(timeout = 1200)

cat("========================================================================\n")
cat("Step 00: R Environment Check & GSE243013 File Manifest\n")
cat("========================================================================\n\n")

## --- 01. Working Directory ------------------------------------------------
cat("[INFO] Current working directory:\n")
tryCatch({
  cat(getwd(), "\n\n")
}, error = function(e) {
  cat("[ERROR] Failed to get working directory:", conditionMessage(e), "\n")
})

## --- 02. R Version --------------------------------------------------------
cat("[INFO] R version:\n")
tryCatch({
  cat(R.version.string, "\n\n")
}, error = function(e) {
  cat("[ERROR] Failed to get R version:", conditionMessage(e), "\n")
})

## --- 03. OS Info ----------------------------------------------------------
cat("[INFO] Operating system info:\n")
tryCatch({
  os_info <- Sys.info()
  cat(paste(names(os_info), os_info, sep = ": "), sep = "\n")
  cat("\n")
}, error = function(e) {
  cat("[ERROR] Failed to get OS info:", conditionMessage(e), "\n")
})

## --- 04. R Library Paths --------------------------------------------------
cat("[INFO] R package library paths:\n")
tryCatch({
  cat(.libPaths(), sep = "\n")
  cat("\n")
}, error = function(e) {
  cat("[ERROR] Failed to get library paths:", conditionMessage(e), "\n")
})

## --- 05. System Memory (macOS) --------------------------------------------
cat("[INFO] System memory (macOS):\n")
tryCatch({
  mem_bytes <- system("sysctl -n hw.memsize", intern = TRUE)
  if (length(mem_bytes) > 0 && !is.na(mem_bytes)) {
    mem_gb <- as.numeric(mem_bytes) / (1024^3)
    cat(sprintf("Total memory: %s bytes (%.1f GB)\n\n", mem_bytes, mem_gb))
  } else {
    cat("Could not read memory info\n\n")
  }
}, error = function(e) {
  cat("[WARNING] Could not read memory info:", conditionMessage(e), "\n\n")
})

## --- 06. Disk Space -------------------------------------------------------
cat("[INFO] Current disk space:\n")
tryCatch({
  disk_info <- system("df -h .", intern = TRUE)
  cat(disk_info, sep = "\n")
  cat("\n")
}, error = function(e) {
  cat("[WARNING] Could not read disk space:", conditionMessage(e), "\n\n")
})

## --- 07. Install BiocManager if needed ------------------------------------
cat("[INFO] Checking BiocManager...\n")
tryCatch({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("[INFO] Installing BiocManager from CRAN...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
    cat("[INFO] BiocManager installed successfully.\n")
  } else {
    cat("[INFO] BiocManager already installed.\n")
  }
}, error = function(e) {
  cat("[ERROR] Failed to install BiocManager:", conditionMessage(e), "\n")
})

## --- 08. Install GEOquery if needed ---------------------------------------
cat("[INFO] Checking GEOquery...\n")
tryCatch({
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    cat("[INFO] Installing GEOquery via BiocManager...\n")
    BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
    cat("[INFO] GEOquery installed successfully.\n")
  } else {
    cat("[INFO] GEOquery already installed.\n")
  }
}, error = function(e) {
  cat("[ERROR] Failed to install GEOquery:", conditionMessage(e), "\n")
})

## --- 09. Get GSE243013 Supplementary File Manifest -------------------------
cat("[INFO] Fetching GSE243013 supplementary file manifest (NO downloads)...\n")
gse243013_manifest <- NULL
tryCatch({
  gse243013_manifest <- GEOquery::getGEOSuppFiles(
    "GSE243013",
    fetch_files = FALSE
  )
  cat("[INFO] GEOquery manifest result:\n")
  if (is.null(gse243013_manifest) || nrow(gse243013_manifest) == 0) {
    cat("[INFO] GEOquery returned empty manifest (FTP access may be restricted).\n")
    cat("[INFO] Falling back to SOFT metadata for file list.\n")
  } else {
    cat(sprintf("[INFO] GEOquery found %d files.\n", nrow(gse243013_manifest)))
  }
}, error = function(e) {
  cat("[WARNING] GEOquery getGEOSuppFiles failed:", conditionMessage(e), "\n")
  cat("[INFO] Falling back to SOFT metadata for file list.\n")
})

## --- 10. Fallback: Parse SOFT metadata for file list ----------------------
cat("[INFO] Parsing SOFT metadata for supplementary file URLs...\n")
soft_urls <- NULL
tryCatch({
  soft_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013&targ=self&form=text"
  soft_text <- readLines(soft_url, warn = FALSE)
  supp_lines <- soft_text[grepl("^!Series_supplementary_file", soft_text)]
  soft_urls <- sub("^!Series_supplementary_file\\s*=\\s*", "", supp_lines)
  cat(sprintf("[INFO] Found %d supplementary file URLs from SOFT metadata.\n", length(soft_urls)))
}, error = function(e) {
  cat("[WARNING] Failed to parse SOFT metadata:", conditionMessage(e), "\n")
})

## --- 11. Build combined manifest ------------------------------------------
cat("[INFO] Building combined manifest...\n")
manifest_rows <- list()

tryCatch({
  if (!is.null(gse243013_manifest) && nrow(gse243013_manifest) > 0) {
    for (i in seq_len(nrow(gse243013_manifest))) {
      fname <- rownames(gse243013_manifest)[i]
      manifest_rows[[fname]] <- data.frame(
        file_name = fname,
        source = "GEOquery",
        url = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
}, error = function(e) {
  cat("[WARNING] Error processing GEOquery results:", conditionMessage(e), "\n")
})

tryCatch({
  if (!is.null(soft_urls) && length(soft_urls) > 0) {
    for (url in soft_urls) {
      fname <- basename(url)
      if (!(fname %in% names(manifest_rows))) {
        manifest_rows[[fname]] <- data.frame(
          file_name = fname,
          source = "SOFT_metadata",
          url = url,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}, error = function(e) {
  cat("[WARNING] Error processing SOFT URLs:", conditionMessage(e), "\n")
})

## --- 12. Check file sizes via HEAD requests --------------------------------
cat("[INFO] Checking file sizes via HTTP HEAD requests...\n")
if (length(manifest_rows) > 0) {
  manifest_df <- do.call(rbind, manifest_rows)
  rownames(manifest_df) <- NULL
  manifest_df$size_bytes <- NA_character_
  manifest_df$status <- NA_character_

  for (i in seq_len(nrow(manifest_df))) {
    fname <- manifest_df$file_name[i]
    # Construct HTTPS download URL
    https_url <- paste0(
      "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243013&format=file&file=",
      fname
    )
    tryCatch({
      h <- curl::new_handle()
      curl::handle_setopt(h, nobody = TRUE, followlocation = TRUE)
      resp <- curl::curl_fetch_memory(https_url, handle = h)
      headers <- rawToChar(resp$headers)
      cl_match <- regmatches(headers, regexpr("(?i)content-length:\\s*\\d+", headers))
      if (length(cl_match) > 0) {
        size_str <- sub("(?i)content-length:\\s*", "", cl_match)
        manifest_df$size_bytes[i] <- size_str
        manifest_df$status[i] <- paste0("HTTP ", resp$status_code)
      } else {
        manifest_df$status[i] <- paste0("HTTP ", resp$status_code, " (no Content-Length)")
      }
      cat(sprintf("  %-55s  %s  %s\n", fname, manifest_df$status[i],
                  ifelse(!is.na(manifest_df$size_bytes[i]),
                         sprintf("(%s bytes)", manifest_df$size_bytes[i]),
                         "")))
    }, error = function(e) {
      manifest_df$status[i] <- paste0("ERROR: ", conditionMessage(e))
      cat(sprintf("  %-55s  %s\n", fname, manifest_df$status[i]))
    })
  }
} else {
  manifest_df <- data.frame(
    file_name = character(0),
    source = character(0),
    url = character(0),
    size_bytes = character(0),
    status = character(0),
    stringsAsFactors = FALSE
  )
  cat("[WARNING] No files found in either source.\n")
}

## --- 13. Save Manifest to CSV ---------------------------------------------
cat("\n[INFO] Saving manifest to CSV...\n")
tryCatch({
  if (nrow(manifest_df) > 0) {
    csv_path <- "02_data/manifest/GSE243013_supplementary_files.csv"
    write.csv(manifest_df, csv_path, row.names = FALSE)
    cat(sprintf("[INFO] Manifest saved to: %s\n", csv_path))
    cat(sprintf("[INFO] Total files in manifest: %d\n", nrow(manifest_df)))
    cat("\n--- File List ---\n")
    for (i in seq_len(nrow(manifest_df))) {
      size_str <- ifelse(!is.na(manifest_df$size_bytes[i]),
                         sprintf("%s bytes", manifest_df$size_bytes[i]),
                         "unknown")
      cat(sprintf("  %s  [%s]  %s\n", manifest_df$file_name[i],
                  manifest_df$status[i], size_str))
    }
    cat("--- End File List ---\n\n")
  } else {
    cat("[WARNING] Manifest is empty. No CSV saved.\n")
  }
}, error = function(e) {
  cat("[ERROR] Failed to save manifest CSV:", conditionMessage(e), "\n")
})

## --- 14. Save sessionInfo -------------------------------------------------
cat("[INFO] Saving sessionInfo()...\n")
tryCatch({
  si <- capture.output(sessionInfo())
  writeLines(si, "logs/sessionInfo.txt")
  cat("[INFO] sessionInfo saved to: logs/sessionInfo.txt\n\n")
}, error = function(e) {
  cat("[ERROR] Failed to save sessionInfo:", conditionMessage(e), "\n")
})

cat("========================================================================\n")
cat("Step 00 completed.\n")
cat("========================================================================\n")
