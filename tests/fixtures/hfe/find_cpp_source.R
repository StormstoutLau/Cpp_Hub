## 查找 highfrequency 包的 C++ 源码
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
pkgPath <- system.file(package = "highfrequency")
cat(sprintf("Package path: %s\n", pkgPath))
cat(sprintf("libs/ subdirs:\n"))
libs_dir <- file.path(pkgPath, "libs")
if (dir.exists(libs_dir)) print(list.dirs(libs_dir, recursive = FALSE, full.names = FALSE))

## 检查是否含 src/ 源码 (R 包安装时通常不含)
src_dir <- file.path(pkgPath, "src")
cat(sprintf("\nsrc/ exists: %s\n", dir.exists(src_dir)))
if (dir.exists(src_dir)) {
  print(list.files(src_dir, recursive = TRUE))
}

## 检查 R 源码文件
R_dir <- file.path(pkgPath, "R")
cat(sprintf("\nR/ files (first 20):\n"))
print(head(list.files(R_dir), 20))

## 尝试获取包源码 tarball URL
cat("\nPackage Description:\n")
print(packageDescription("highfrequency")[c("Package", "Version", "Repository", "URL")])

## 在 R 环境中搜索 .cpp 或 .c 文件
cat("\nSearching for source files in package dir:\n")
all_files <- list.files(pkgPath, recursive = TRUE, full.names = FALSE)
src_files <- all_files[grep("\\.(cpp|c|f|h|hpp)$", all_files, ignore.case = TRUE)]
cat(sprintf("Found %d C/C++/Fortran source files:\n", length(src_files)))
if (length(src_files) > 0) print(head(src_files, 30))

cat("\nDONE.\n")
