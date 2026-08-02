## 下载 highfrequency 源码 tarball, 提取 kernelEstimator.cpp
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))

pkgVersion <- as.character(packageDescription("highfrequency")$Version)
cat(sprintf("highfrequency version: %s\n", pkgVersion))

tarball_url <- sprintf("https://cran.r-project.org/src/contrib/highfrequency_%s.tar.gz", pkgVersion)
cat(sprintf("Downloading: %s\n", tarball_url))

dest_dir <- file.path(getwd(), "tests", "fixtures", "hfe", "hf_source")
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
tarball_path <- file.path(dest_dir, basename(tarball_url))

download.file(tarball_url, tarball_path, mode = "wb", quiet = FALSE)
cat(sprintf("Downloaded to: %s (%.1f KB)\n", tarball_path, file.size(tarball_path)/1024))

## 解压 src/ 目录
cat("\n=== Untar src/ only ===\n")
untar(tarball_path, exdir = dest_dir, files = "highfrequency/src/")
src_dir <- file.path(dest_dir, "highfrequency", "src")
cat(sprintf("src/ files:\n"))
print(list.files(src_dir))

## 读取 kernelEstimator.cpp
cpp_path <- file.path(src_dir, "kernelEstimator.cpp")
if (file.exists(cpp_path)) {
  cat(sprintf("\n=== kernelEstimator.cpp (%.1f KB) ===\n", file.size(cpp_path)/1024))
  cpp_lines <- readLines(cpp_path)
  cat(sprintf("Total lines: %d\n\n", length(cpp_lines)))
  cat("--- First 250 lines ---\n")
  cat(paste(cpp_lines[1:min(250, length(cpp_lines))], collapse = "\n"))
  cat("\n")
} else {
  cat("kernelEstimator.cpp NOT FOUND\n")
  cat("Available files:\n")
  print(list.files(src_dir, pattern = "\\.cpp$"))
}

cat("\nDONE.\n")
