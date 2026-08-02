## =============================================================================
## verify_v141_functions.R
## v1.4.1 R 函数可用性验证 (排除幻觉的第一步)
##
## 目标: 验证 highfrequency 1.0.3 中以下函数是否存在, 并打印签名
##   - rKernelCov (Realized Kernel)
##   - sparseSampling (稀疏抽样)
##   - noiseBPM (噪声方差估计)
##   - noiseAC (噪声自相关)
##   - optimFrequ (最优抽样频率)
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

cat("=== highfrequency version:", as.character(packageVersion("highfrequency")), "===\n\n")

## ---- 函数存在性检查 ----
funcs <- c("rKernelCov", "sparseSampling", "noiseBPM", "noiseAC",
           "optimFrequ", "rKernelCov", "rAVGCov", "rTSCov",
           "noise vars", "QLRV")
cat("=== Function existence check ===\n")
for (f in c("rKernelCov", "sparseSampling", "noiseBPM", "noiseAC", "optimFrequ",
            "rAVGCov", "rTSCov", "rHYCov", "rRVar", "rBPCov")) {
  exists <- exists(f, where=asNamespace("highfrequency"), inherits=FALSE)
  cat(sprintf("  %-20s %s\n", f, ifelse(exists, "EXISTS", "MISSING")))
}

## ---- 函数签名打印 ----
cat("\n=== Function signatures (args) ===\n")
for (f in c("rKernelCov", "sparseSampling", "noiseBPM", "noiseAC", "optimFrequ")) {
  if (exists(f, where=asNamespace("highfrequency"), inherits=FALSE)) {
    cat(sprintf("\n--- %s ---\n", f))
    cat(deparse(args(getFromNamespace(f, "highfrequency"))), "\n")
  }
}

## ---- 实际调用测试 (用 sampleTData) ----
cat("\n=== Functional test with sampleTData ===\n")
data(sampleTData)

## sparseSampling 测试
cat("\n--- sparseSampling ---\n")
tryCatch({
  ss <- sparseSampling(sampleTData, alignBy = "minutes", alignPeriod = 5)
  cat(sprintf("  sparseSampling OK: %d obs -> %d obs\n",
              nrow(sampleTData), nrow(ss)))
  cat(sprintf("  first PRICE: %.4f  last PRICE: %.4f\n",
              as.numeric(ss$PRICE[1]), as.numeric(ss$PRICE[nrow(ss)])))
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

## rKernelCov 测试 (需要收益率输入)
cat("\n--- rKernelCov ---\n")
tryCatch({
  ret <- makeReturns(sampleTData$PRICE)
  rk <- rKernelCov(ret, kernel = "Bartlett", bandWidth = 5)
  cat(sprintf("  rKernelCov OK: value = %.10f\n", as.numeric(rk)))
  cat(sprintf("  class: %s\n", paste(class(rk), collapse=", ")))
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

## noiseBPM 测试
cat("\n--- noiseBPM ---\n")
tryCatch({
  nbpm <- noiseBPM(sampleTData)
  cat(sprintf("  noiseBPM OK:\n"))
  print(nbpm)
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

## noiseAC 测试
cat("\n--- noiseAC ---\n")
tryCatch({
  nac <- noiseAC(sampleTData, lags = 5)
  cat(sprintf("  noiseAC OK:\n"))
  print(nac)
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

## optimFrequ 测试
cat("\n--- optimFrequ ---\n")
tryCatch({
  of <- optimFrequ(sampleTData, k = 5, freq = "sec")
  cat(sprintf("  optimFrequ OK:\n"))
  print(of)
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("\nDONE.\n")
