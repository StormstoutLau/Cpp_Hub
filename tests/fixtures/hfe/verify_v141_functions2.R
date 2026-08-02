## =============================================================================
## verify_v141_functions2.R
## 搜索 highfrequency 1.0.3 中实际的噪声/抽样相关函数
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

cat("=== All highfrequency exported functions ===\n")
funcs <- ls("package:highfrequency")
cat(sprintf("Total: %d functions\n\n", length(funcs)))

## ---- 搜索噪声相关 ----
cat("=== Noise-related functions ===\n")
noise_funcs <- funcs[grep("noise|Noise|NOISE", funcs)]
print(noise_funcs)

## ---- 搜索抽样相关 ----
cat("\n=== Sampling-related functions ===\n")
samp_funcs <- funcs[grep("samp|Samp|SAMP|sparse|Sparse|SPARSE|aggreg|Aggreg", funcs)]
print(samp_funcs)

## ---- 搜索频率相关 ----
cat("\n=== Frequency-related functions ===\n")
freq_funcs <- funcs[grep("freq|Freq|FREQ|optim|Optim|OPTIM", funcs)]
print(freq_funcs)

## ---- 搜索 Kernel 相关 ----
cat("\n=== Kernel-related functions ===\n")
kern_funcs <- funcs[grep("ernel|ERNEL", funcs)]
print(kern_funcs)

## ---- 搜索 liquidity 相关 (v1.4.3 预查) ----
cat("\n=== Liquidity-related functions ===\n")
liq_funcs <- funcs[grep("liqu|Liqu|LIQU|spread|Spread|SPREAD|amihud|Amihud|AMIHUD", funcs)]
print(liq_funcs)

## ---- 搜索 HAR/HEAVY 相关 (v1.4.2 预查) ----
cat("\n=== HAR/HEAVY-related functions ===\n")
model_funcs <- funcs[grep("HAR|HEAVY|har|heavy|forecast|Forecast", funcs)]
print(model_funcs)

## ---- 搜索 jump 相关 (已有 + 扩展) ----
cat("\n=== Jump-related functions ===\n")
jump_funcs <- funcs[grep("ump|UMP", funcs)]
print(jump_funcs)

## ---- rKernelCov 正确调用 ----
cat("\n=== rKernelCov correct call ===\n")
data(sampleTData)
tryCatch({
  ret <- makeReturns(sampleTData$PRICE)
  rk <- rKernelCov(rData = ret, kernelType = "Bartlett", kernelParam = 5)
  cat(sprintf("  rKernelCov OK: value = %.10f\n", as.numeric(rk)))
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

## ---- rAVGCov / rTSCov 签名 ----
cat("\n=== rAVGCov signature ===\n")
cat(deparse(args(rAVGCov)), "\n")
cat("\n=== rTSCov signature ===\n")
cat(deparse(args(rTSCov)), "\n")

cat("\nDONE.\n")
