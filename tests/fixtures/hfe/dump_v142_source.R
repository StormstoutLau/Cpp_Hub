.libPaths(c(file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6"), .libPaths()))
suppressMessages(library(highfrequency))

fns <- c("rAVGCov", "rTSCov", "rRTSCov", "rMRCov", "rHYCov", "HARmodel", "HEAVYmodel")

for (fn in fns) {
  cat(sprintf("\n%s\n", paste(rep("=", 80), collapse="")))
  cat(sprintf("=== %s ===\n", fn))
  cat(paste(rep("=", 80), collapse=""), "\n")
  obj <- get(fn, envir=asNamespace("highfrequency"))
  # 方法: getAnywhere 找所有定义
  ga <- getAnywhere(fn)
  # 打印非命名空间隐藏的对象
  for (i in seq_along(ga$objs)) {
    o <- ga$objs[[i]]
    if (is.function(o)) {
      cat(sprintf("\n--- [%s] from: %s ---\n", i, ga$where[i]))
      # 打印完整源码
      src <- deparse(o)
      cat(src, sep="\n")
      cat(sprintf("\n--- end [%s] (length: %d lines) ---\n", i, length(src)))
    }
  }
}
