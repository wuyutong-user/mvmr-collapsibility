if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE)
}
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1]], mustWork = TRUE)
dir.create(file.path(root, "environment"), showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(root, "environment", "current_local_sessionInfo.txt"))
pkgs <- c("TwoSampleMR", "MVMR", "MRSL", "data.table", "ieugwasr", "igraph",
          "ggplot2", "dplyr", "tidyr", "readr", "patchwork", "scales", "svglite",
          "jsonlite", "ggrepel", "ragg")
versions <- vapply(pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p))
  else "not installed in current local default library"
}, character(1))
write.csv(data.frame(package = pkgs, current_local_version = versions),
          file.path(root, "environment", "current_local_package_versions.csv"), row.names = FALSE)
rfiles <- list.files(root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
for (f in rfiles) parse(file = f)
cat(sprintf("PARSE PASS: %d archived R files.\n", length(rfiles)))
cat("Current local environment recorded separately from historical application environment.\n")
