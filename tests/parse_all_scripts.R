#!/usr/bin/env Rscript
# Test: Parse all R scripts in the repository
cat("=== Parse All Scripts Test ===\n")

repo <- Sys.getenv("GSE243013_PROJECT_ROOT", unset = normalizePath("..", mustWork = FALSE))
if (!dir.exists(file.path(repo, "scripts"))) {
  repo <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."))
}

dirs <- c("scripts/analysis", "scripts/manuscript")
scripts <- c()
for (d in dirs) {
  fp <- file.path(repo, d)
  if (dir.exists(fp)) {
    scripts <- c(scripts, list.files(fp, pattern="\\.R$", full.names=TRUE))
  }
}
# Also check run_pipeline.R
rp <- file.path(repo, "scripts/run_pipeline.R")
if (file.exists(rp)) scripts <- c(scripts, rp)

cat("Scripts found:", length(scripts), "\n")
errors <- data.frame(file=character(), error=character(), stringsAsFactors=FALSE)

for (s in scripts) {
  result <- tryCatch(parse(file=s), error=function(e) e)
  if (inherits(result, "error")) {
    errors <- rbind(errors, data.frame(file=basename(s), error=result$message, stringsAsFactors=FALSE))
    cat("  FAIL:", basename(s), "-", result$message, "\n")
  }
}

if (nrow(errors) > 0) {
  cat("\nFAILED:", nrow(errors), "scripts\n")
  write.csv(errors, file.path(repo, "docs", "parse_errors.csv"), row.names=FALSE)
  quit(status=1)
} else {
  cat("\nPASS: All", length(scripts), "scripts parse successfully\n")
  quit(status=0)
}
