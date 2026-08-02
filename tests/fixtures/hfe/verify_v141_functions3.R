## =============================================================================
## verify_v141_functions3.R
## 深度搜索 highfrequency 1.0.3 中所有 r* 系列估计量 + HAR/HEAVY + liquidity 签名
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

funcs <- ls("package:highfrequency")
cat(sprintf("Total exported: %d\n\n", length(funcs)))

## ---- 所有 r 开头函数 (realized 系列估计量) ----
cat("=== All r* functions ===\n")
r_funcs <- funcs[grep("^r[A-Z]", funcs)]
print(r_funcs)

## ---- 所有 R 开头函数 ----
cat("\n=== All R* functions ===\n")
R_funcs <- funcs[grep("^R[A-Z]", funcs)]
print(R_funcs)

## ---- realized / Realized 相关 ----
cat("\n=== realized-related ===\n")
print(funcs[grep("ealized|REALIZED", funcs)])

## ---- jump test 全部签名 ----
cat("\n=== Jump test signatures ===\n")
for (fn in c("BNSjumpTest", "AJjumpTest", "JOjumpTest", "rankJumpTest", "intradayJumpTest")) {
  cat(sprintf("\n--- %s ---\n", fn))
  print(args(get(fn)))
}

## ---- rKernelCov 完整签名 + listAvailableKernels ----
cat("\n--- rKernelCov ---\n")
print(args(rKernelCov))
cat("\n--- listAvailableKernels ---\n")
print(listAvailableKernels())

## ---- rAVGCov / rTSCov / rBPCov / rCov / rQuar / rSV / rSVar 签名 ----
cat("\n=== r* covariance/quarticity signatures ===\n")
for (fn in c("rCov","rBPCov","rQuar","rRTSCov","rTSCov","rAVGCov","rKernelCov",
             "rSVar","rSV","rRVar","rBPVar","rTPQuar","rOPVar")) {
  if (exists(fn, where = "package:highfrequency", inherits = FALSE)) {
    cat(sprintf("\n--- %s ---\n", fn))
    print(args(get(fn, envir = asNamespace <- getNamespace("highfrequency"))))
  } else {
    cat(sprintf("\n--- %s : NOT FOUND ---\n", fn))
  }
}

## ---- HARmodel / HEAVYmodel 签名 ----
cat("\n=== HARmodel signature ===\n")
print(args(HARmodel))
cat("\n=== HEAVYmodel signature ===\n")
print(args(HEAVYmodel))

## ---- getLiquidityMeasures 签名 ----
cat("\n=== getLiquidityMeasures signature ===\n")
print(args(getLiquidityMeasures))

## ---- 示例数据集 ----
cat("\n=== Sample datasets ===\n")
print(data(package = "highfrequency")$results[, "Item"])

cat("\nDONE.\n")
