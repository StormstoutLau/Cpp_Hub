## =============================================================================
## verify_kernel_sources.R
## 查 highfrequency 1.0.3 核函数内部实现, 排除 spec 公式幻觉
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

cat("=== rKernelCov 源码 (前 80 行) ===\n")
print(getAnywhere("rKernelCov"))

cat("\n\n=== 搜索核函数相关内部函数 ===\n")
ns <- asNamespace("highfrequency")
all_internal <- ls(ns)
kernel_funcs <- all_internal[grep("ernel|kfunc|kVal|kernelWeights|rectangular|Bartlett|Parzen|Tukey|Epanechnikov|Cubic|Fifth|Sixth|Seventh|Eighth|Second",
                                   all_internal, ignore.case = TRUE)]
print(kernel_funcs)

cat("\n\n=== 逐个查核函数源码 ===\n")
for (fn in kernel_funcs) {
  obj <- get(fn, envir = ns, inherits = FALSE)
  cat(sprintf("\n--- %s (class: %s) ---\n", fn, class(obj)[1]))
  if (is.function(obj)) {
    cat("ARGS: "); print(args(obj))
    cat("BODY:\n")
    print(body(obj))
  } else {
    print(obj)
  }
}

## ---- 用各核函数对 x=0/0.25/0.5/0.75/1.0 直接采样数值 ----
cat("\n\n=== 各核函数数值采样 (用 rKernelCov 内部测试) ===\n")
set.seed(123)
prices <- 100 + cumsum(rnorm(50) * 0.001)
ret <- diff(log(prices))
cat(sprintf("ret[1:5]: %s\n", paste(round(ret[1:5], 6), collapse=", ")))
cat(sprintf("RV = %.10f, n = %d\n\n", sum(ret^2), length(ret)))

for (kn in listAvailableKernels()) {
  rk <- rKernelCov(rData = ret, kernelType = kn, kernelParam = 1,
                   kernelDOFadj = FALSE)
  cat(sprintf("  %-25s RK (no DOF adj) = %.10f\n", kn, as.numeric(rk)))
}

cat("\n=== WITH DOF adjustment ===\n")
for (kn in listAvailableKernels()) {
  rk <- rKernelCov(rData = ret, kernelType = kn, kernelParam = 1,
                   kernelDOFadj = TRUE)
  cat(sprintf("  %-25s RK (DOF adj)   = %.10f\n", kn, as.numeric(rk)))
}

cat("\nDONE.\n")
