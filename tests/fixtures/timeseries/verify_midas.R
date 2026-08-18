# verify_midas.R - MIDAS 基准生成: midasr 0.9 (Phase 7C M4)
#
# 对照点 (spec §1.3):
#   nealmon/nbeta 权重 vs midasr 逐点          1e-12
#   U-MIDAS vs midasr midas_u (纯 OLS)          1e-10
#   MIDAS NLS vs midas_r 夹具 (收紧 control)    1e-6~1e-8
#   hAh_test 三列全录 (K-Z 2012, MD6)
#
# 实证任务 (决策 22/24):
#   W-dir: 权重 w_i 与 mls 列 h 的对应方向裁决 — 构造精确线性 DGP
#     y_t = μ + Σ_i w_i(λ*)·x_{tm−(i−1)} (假设: w_1 ↔ lag0 最新)
#     midas_r 恢复 λ* 则方向定案 → probe 实测: Form A 恢复 (5.03, -0.499)
#   mls 期末对齐实测 (MD3)
#
# 用法: Rscript verify_midas.R
.libPaths(c('F:/R/win-library/4.6', file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6"), .libPaths()))
suppressMessages(library(midasr))

here <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", here[grepl("^--file=", here)])))
set.seed(42)

fmt17 <- function(x) sprintf("%.17g", x)

cat("=== W1: 权重逐点基准 (midasr 源函数直调, %.17g) ===\n")
cat("nealmon(c(1, -0.5), 4):\n")
for (v in nealmon(c(1, -0.5), 4, NULL)) cat("  ", fmt17(v), "\n")
cat("nealmon(c(0.8, 0.1, -0.02), 6):\n")
for (v in nealmon(c(0.8, 0.1, -0.02), 6, NULL)) cat("  ", fmt17(v), "\n")
cat("nbetaMT(c(1, 2, 3, 0), 5):\n")
for (v in nbetaMT(c(1, 2, 3, 0), 5, NULL)) cat("  ", fmt17(v), "\n")
cat("nbetaMT(c(0.5, 1.5, 1.5, 0.1), 7):\n")
for (v in nbetaMT(c(0.5, 1.5, 1.5, 0.1), 7, NULL)) cat("  ", fmt17(v), "\n")
cat("almonp(c(0.1, 0.05, -0.01), 5):\n")
for (v in almonp(c(0.1, 0.05, -0.01), 5, NULL)) cat("  ", fmt17(v), "\n")
cat("polystep(c(0.5, 0.2, 0.1), c(2, 5), 8):\n")
for (v in polystep(c(0.5, 0.2, 0.1), 8, NULL, c(2, 5))) cat("  ", fmt17(v), "\n")
cat("harstep(c(0.6, 0.3, 0.1), 20):\n")
hv <- harstep(c(0.6, 0.3, 0.1), 20, NULL)
cat("  ", fmt17(hv[1]), fmt17(hv[2]), fmt17(hv[6]), fmt17(hv[20]), "\n")

cat("\n=== W2: mls 期末对齐实测 (MD3) ===\n")
m12 <- mls(1:12, 0:2, 3)
cat("mls(1:12, 0:2, 3) lag0 列:", paste(m12[, 1], collapse = ","), "\n")
cat("mls(1:12, 0:2, 3) lag1 列:", paste(m12[, 2], collapse = ","), "\n")
m13 <- mls(1:12, 0:3, 3)
cat("mls(1:12, 0:3, 3) 行数 (burn-in 丢期):", nrow(m13),
    " lag0 首行:", m13[1, 1], "\n")

# ---------------------------------------------------------------------------
# DGP (固定, 写 CSV): m=3, n=100 低频期, x 高频 300 + 3 预样本
# 方向裁决: y_t = 2 + Σ_{i=1..4} w_i·x_{3t−(i−1)}, w = nealmon(c(5,−0.5),4)
#   即 w_1 ↔ h=0 (最新/期末), w_4 ↔ h=3 (最旧) — 已由 Form A 恢复 λ* 定案
# 对齐 (MD3 关键): mls(x_use, 0:3, 3) 行 t 列 h = x_use[3t−h] (期末系);
#   期 1 行 (x_use[0]) NA → na.omit 丢 → 估计样本 = 期 2..100 (99 行)。
#   y_1 用 x_ext 预样本生成真实值 (供 AR 列期 2 用), 非占位。
#   CSV 的 y/x 与 C++ MixedFreqData 完全同数组 (100 + 300), C++ 端
#   design_matrix 行 j=2..100 ↔ R 期 2..100 逐行一致。
# U-MIDAS 真值: 同结构 β_i = w_i (k_high=4), intercept=2
# NLS 夹具: 同 DGP 但 λ 待估 (应恢复 ≈ (5, −0.5))
# ---------------------------------------------------------------------------
n_lf <- 100
m <- 3
x_ext <- as.numeric(arima.sim(list(ar = 0.5), n = n_lf * m + m,
                              rand.gen = function(n) rnorm(n, 0, 0.8)))
x_use <- x_ext[(m + 1):length(x_ext)]    # 300: 期 t 期末 = x_use[3t]
w_true <- nealmon(c(5, -0.5), 4, NULL)
xwin <- function(t, h) {
  if (m * t - h >= 1) x_use[m * t - h] else x_ext[m * t - h + m]
}
e_all <- rnorm(n_lf, 0, 0.3)
y_full <- numeric(n_lf)
for (t in 1:n_lf) {
  y_full[t] <- 2 + sum(sapply(0:3, function(h) w_true[h + 1] * xwin(t, h))) + e_all[t]
}

csv <- file.path(here, "midas_smoke_data.csv")
# 低频宽表: 100 行 × 4 列 (y_t, x_{3t−2}, x_{3t−1}, x_{3t}); C++ 端逐行
# 展开 col1..3 即恢复高频序列原序 (与 MixedFreqData y/x 同数组)
xm <- matrix(x_use, ncol = m, byrow = TRUE)
lines <- character(n_lf)
for (t in 1:n_lf) {
  lines[t] <- paste(c(fmt17(y_full[t]), fmt17(xm[t, ])), collapse = ",")
}
writeLines(c(lines, "#columns=y,x1,x2,x3"), csv)
cat("\nDGP 写入:", csv, " (w_true =", paste(fmt17(w_true), collapse = ", "),
    ")\n")

cat("\n=== W3: midas_u 纯 OLS (1e-10 主锚) ===\n")
mu <- midas_u(y_full ~ mls(x_use, 0:3, m))
n_eff <- length(residuals(mu))
cat("n_eff:", n_eff, "\n")
cat("coef:", paste(fmt17(coef(mu)), collapse = ", "), "\n")
cat("sigma2:", fmt17(sum(residuals(mu)^2) / (n_eff - 5)), "\n")
cat("SSR:", fmt17(sum(residuals(mu)^2)), "\n")

cat("\n=== W4: 方向裁决 + midas_r NLS (reltol=1e-12) ===\n")
# Form A (probe_midas_form.R 实测): 权重函数作 mls 第4参, start 给纯参数
# 退化陷阱: 权重放 start → midasr 静默按 U-MIDAS (lm) 拟合
ctrl <- list(reltol = 1e-12, maxit = 10000)
mr <- midas_r(y_full ~ mls(x_use, 0:3, m, nealmon),
              start = list(x_use = c(5, -0.5)),
              control = ctrl)
cat("DL coef:", paste(fmt17(coef(mr)), collapse = ", "), "\n")
cat("DL convergence:", mr$convergence, "\n")
cat("DL SSR:", fmt17(deviance(mr)), "\n")
cat("DL midas.coef (w1..w4):",
    paste(fmt17(mr$midas_coefficients), collapse = ", "), "\n")

cat("\n=== W5: midas_r DL 随机起始 (对照 C++ 多起点逃逸) ===\n")
mr2 <- midas_r(y_full ~ mls(x_use, 0:3, m, nealmon),
               start = list(x_use = c(1, 0)),
               control = ctrl)
cat("DL2 coef:", paste(fmt17(coef(mr2)), collapse = ", "),
    " conv:", mr2$convergence, "\n")

cat("\n=== W6: MIDAS-AR (含 y 自身滞后) ===\n")
y_ar <- y_full
mar <- midas_r(y_ar ~ mls(x_use, 0:3, m, nealmon) + mls(y_ar, 1:1, 1),
               start = list(x_use = c(5, -0.5), y_ar = c(0.1)),
               control = ctrl)
cat("AR coef:", paste(fmt17(coef(mar)), collapse = ", "),
    " conv:", mar$convergence, "\n")

cat("\n=== W7: hAh_test 三列 (K-Z 2012, MD6) ===\n")
ha <- hAh_test(mr)
cat("hAh stat:", fmt17(as.numeric(ha$statistic)),
    " p:", fmt17(as.numeric(ha$p.value)),
    " df:", paste(ha$parameter, collapse = ","), "\n")
