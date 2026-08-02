.libPaths(c(file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6"), .libPaths()))
suppressMessages(library(highfrequency))
cat("=== highfrequency version ===\n")
cat(as.character(packageVersion("highfrequency")), "\n\n")

fns <- c("rAVGCov", "rTSCov", "rRTSCov", "rMRCov", "rHYCov", "HARmodel", "HEAVYmodel")
cat("=== function availability ===\n")
for (fn in fns) {
  cat(sprintf("%-15s: %s\n", fn, ifelse(exists(fn, where=asNamespace("highfrequency")), "OK", "MISSING")))
}

cat("\n=== function signatures ===\n")
for (fn in fns) {
  if (exists(fn, where=asNamespace("highfrequency"))) {
    cat(sprintf("\n--- %s ---\n", fn))
    fn_obj <- get(fn, envir=asNamespace("highfrequency"))
    cat(deparse(args(fn_obj)), sep="\n")
  }
}
