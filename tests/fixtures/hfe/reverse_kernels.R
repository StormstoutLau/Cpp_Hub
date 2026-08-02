## =============================================================================
## reverse_kernels.R
## 用已知序列反推 highfrequency 1.0.3 各核函数 k(x) 数值, 排除 spec 公式幻觉
##
## 原理:
##   序列 r = [1, 1, 0, 0, ..., 0], 长度 n
##   γ_0 = 2, γ_1 = 1, γ_h = 0 (h >= 2)
##   RK(H) = γ_0 + Σ_{h=1}^{H} k(h/H) * (γ_h + γ_{-h})
##         = 2 + 2 * k(1/H) * 1   (当 H >= 2, 因 γ_h = 0 for h>=2)
##         = 2 + 2 * k(1/H)
##   => k(1/H) = (RK - 2) / 2     (no DOF adj)
##
##   特殊: H=1 时 RK = 2 + 2*k(1), 即 k(1) = (RK - 2) / 2
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

## 构造 r = [1, 1, 0, ..., 0], n = 30 (足够大, 避免边界)
n <- 30
ret <- c(1, 1, rep(0, n - 2))
cat(sprintf("序列: r = [1, 1, 0, ..., 0], n = %d\n", n))
cat(sprintf("预期: γ_0 = 2, γ_1 = 1, γ_h = 0 (h>=2)\n"))
cat(sprintf("公式: k(x) = (RK(H=1/x) - 2) / 2, no DOF adj\n\n"))

kernels <- listAvailableKernels()
x_grid <- c(1.0, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10)
H_grid <- round(1 / x_grid)
H_grid[H_grid < 1] <- 1

cat(sprintf("%-25s %s\n", "Kernel", paste(sprintf("x=%.4f", x_grid), collapse="  ")))
cat(paste(rep("-", 25 + 10*12), collapse=""), "\n")

reverse_ktable <- data.frame(x = x_grid, H = H_grid)
for (kn in kernels) {
  kvals <- numeric(length(H_grid))
  for (i in seq_along(H_grid)) {
    H <- H_grid[i]
    rk <- rKernelCov(rData = ret, kernelType = kn, kernelParam = H,
                     kernelDOFadj = FALSE)
    rk_val <- as.numeric(rk)
    kvals[i] <- (rk_val - 2) / 2
  }
  reverse_ktable[[kn]] <- kvals
  cat(sprintf("%-25s %s\n", kn,
              paste(sprintf("%10.6f", kvals), collapse="  ")))
}

cat("\n\n=== 关键检查点 ===\n")
cat("k(0) 应该 = 1 (所有核, 但 H=∞ 不可达, 用 H=10 近似 x=0.1)\n")
cat("k(1) 应该 = 0 (除 Rectangular = 1)\n")
cat("支撑 |x| > 1 应该 = 0\n\n")

## 单独验证 k(1) - 用 H=1
cat("=== k(1) 验证 (H=1, x=1/1=1.0) ===\n")
for (kn in kernels) {
  rk <- rKernelCov(rData = ret, kernelType = kn, kernelParam = 1,
                   kernelDOFadj = FALSE)
  k1 <- (as.numeric(rk) - 2) / 2
  cat(sprintf("  %-25s k(1.0) = %+.6f\n", kn, k1))
}

## 验证 k(0.5) - 用 H=2
cat("\n=== k(0.5) 验证 (H=2, x=1/2=0.5) ===\n")
for (kn in kernels) {
  rk <- rKernelCov(rData = ret, kernelType = kn, kernelParam = 2,
                   kernelDOFadj = FALSE)
  k05 <- (as.numeric(rk) - 2) / 2
  cat(sprintf("  %-25s k(0.5) = %+.6f\n", kn, k05))
}

## DOF 调整验证
cat("\n=== DOF 调整验证 (kernelParam=2, n=30) ===\n")
cat("预期: RK_adj = RK_raw * n / (n - H) = RK_raw * 30/28\n\n")
for (kn in kernels) {
  rk_raw <- as.numeric(rKernelCov(rData = ret, kernelType = kn, kernelParam = 2,
                                  kernelDOFadj = FALSE))
  rk_adj <- as.numeric(rKernelCov(rData = ret, kernelType = kn, kernelParam = 2,
                                  kernelDOFadj = TRUE))
  ratio <- rk_adj / rk_raw
  expected_ratio <- n / (n - 2)
  cat(sprintf("  %-25s raw=%.6f  adj=%.6f  ratio=%.6f  expected=%.6f\n",
              kn, rk_raw, rk_adj, ratio, expected_ratio))
}

## 第二种序列验证: r = [1, 2, 3, 0, ..., 0] 测试 γ_1 和 γ_2
cat("\n\n=== 第二序列验证: r = [1, 2, 3, 0, ..., 0] ===\n")
ret2 <- c(1, 2, 3, rep(0, n - 3))
cat(sprintf("γ_0 = 1+4+9 = 14, γ_1 = r2*r1 + r3*r2 = 2*1+3*2 = 8\n"))
cat(sprintf("γ_2 = r3*r1 = 3*1 = 3, γ_h = 0 (h>=3)\n"))
cat(sprintf("RK(H) = 14 + 2*k(1/H)*8 + 2*k(2/H)*3\n\n"))

cat(sprintf("%-25s %s\n", "Kernel",
            paste(sprintf("H=%d", H_grid), collapse="  ")))
for (kn in kernels) {
  rks <- numeric(length(H_grid))
  for (i in seq_along(H_grid)) {
    rk <- rKernelCov(rData = ret2, kernelType = kn, kernelParam = H_grid[i],
                     kernelDOFadj = FALSE)
    rks[i] <- as.numeric(rk)
  }
  cat(sprintf("%-25s %s\n", kn,
              paste(sprintf("%10.6f", rks), collapse="  ")))
}

## 用第二序列反推 k(0.5), k(1.0) 验证一致性
## H=2: RK = 14 + 2*k(0.5)*8 + 2*k(1.0)*3
## H=1: RK = 14 + 2*k(1.0)*8
## => k(1.0) = (RK(H=1) - 14) / 16
## => k(0.5) = (RK(H=2) - 14 - 2*k(1.0)*3) / 16
cat("\n=== 第二序列反推 k(0.5), k(1.0) 与第一序列对照 ===\n")
for (kn in kernels) {
  rk_H1 <- as.numeric(rKernelCov(rData = ret2, kernelType = kn, kernelParam = 1,
                                 kernelDOFadj = FALSE))
  rk_H2 <- as.numeric(rKernelCov(rData = ret2, kernelType = kn, kernelParam = 2,
                                 kernelDOFadj = FALSE))
  k1_seq2 <- (rk_H1 - 14) / 16
  k05_seq2 <- (rk_H2 - 14 - 2 * k1_seq2 * 3) / 16

  rk_H1_s1 <- as.numeric(rKernelCov(rData = ret, kernelType = kn, kernelParam = 1,
                                    kernelDOFadj = FALSE))
  rk_H2_s1 <- as.numeric(rKernelCov(rData = ret, kernelType = kn, kernelParam = 2,
                                    kernelDOFadj = FALSE))
  k1_seq1 <- (rk_H1_s1 - 2) / 2
  k05_seq1 <- (rk_H2_s1 - 2) / 2

  cat(sprintf("  %-25s k(1.0): seq1=%+.6f  seq2=%+.6f  |  k(0.5): seq1=%+.6f  seq2=%+.6f\n",
              kn, k1_seq1, k1_seq2, k05_seq1, k05_seq2))
}

cat("\nDONE.\n")
