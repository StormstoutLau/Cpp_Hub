# -*- coding: utf-8 -*-
# verify_var.R - VAR 交叉基准: R vars::VAR + VARselect (Phase 7C M2)
#
# 对照点 (spec §5.1.5):
#   IC vs R vars::VARselect 1e-8 (交叉, 非主基准)
#   系数 vs vars::VAR (逐方程 OLS 等价性交叉)
#   fevd(orth) vs statsmodels fevd (Cholesky 轨交叉)
#
# 语义差异留档 (vars vs statsmodels):
#   - vars::VARselect 同样截断到 maxlags 统一样本 (与 statsmodels
#     offset=maxlags−p 等价); 但 vars 的 FPE 定义 = det(Σ̃·(n+p·K)/n)·(n/(n−p·K−1))^(1/K)?
#     实测对照记录, 不作 1e-8 主锚 (V6 以 statsmodels 为主)
#   - vars::VAR 的 sigma_u = SSR/(T − K·p − k_trend) (df 修正, 同 statsmodels)
#   - vars fevd 用 chol(sigma_u) 上三角转置 (t(chol) = 下三角), 与
#     np.linalg.cholesky 数值一致 (V2 断言)
#
# 用法: Rscript verify_var.R  (读 var_smoke_data.csv, 输出 12 位全精度)
.libPaths(c("F:/R/win-library/4.6", .libPaths()))
suppressMessages(library(vars))
setwd("F:/Cpp_Hub/tests/fixtures/timeseries")

d <- read.csv("var_smoke_data.csv")
y <- as.matrix(d[, c("y1", "y2", "y3")])
cat(sprintf("T=%d K=%d\n", nrow(y), ncol(y)))

# --- VAR(2) trend=const ---
v2 <- VAR(y, p = 2, type = "const")
coefs <- sapply(v2$varresult, function(m) coef(m))
dimnames(coefs) <- NULL
cat("== coef matrix (K x (K*p + k_trend), cols per-eq) ==\n")
cat("coef:\n"); print(coefs, digits = 12)
sig <- crossprod(resid(v2)) / v2$obs
cat("sigma_u_mle (SSR/T):\n"); print(sig, digits = 12)
cat("loglik:", sprintf("%.12f", logLik(v2)), "\n")
# vars AIC/SC/HQ (logdet form, m = K^2*p + K*k_trend 与 statsmodels fp 同)
n <- v2$obs; K <- 3
fp <- 2 * K^2 + K * 1
ld <- determinant(sig)$modulus[1]
cat("aic:", sprintf("%.12f", ld + 2 / n * fp), "\n")
cat("bic:", sprintf("%.12f", ld + log(n) / n * fp), "\n")
cat("hqic:", sprintf("%.12f", ld + 2 * log(log(n)) / n * fp), "\n")
nstar <- n + (2 * K + 1)  # df_model 单方程含截距
cat("fpe:", sprintf("%.12f", ((n + 2 * K + 1) / (n - 2 * K - 1))^K * exp(ld)), "\n")

# --- Cholesky (V2: R t(chol) == numpy chol 下三角) ---
P <- t(chol(crossprod(resid(v2)) / (n - 2 * K - 1)))
cat("chol_sigma_u lower:\n"); print(P, digits = 12)

# --- FEVD orth (H=10) ---
fv <- fevd(v2, n.ahead = 10)
fevd10 <- t(sapply(fv, function(m) m[10, ]))
cat("fevd orth H=10:\n"); print(fevd10, digits = 12)

# --- IRF orth H=10 (Ψ_10) ---
ir <- irf(v2, n.ahead = 10, ortho = TRUE, ci = 0.95)
psi10 <- do.call(rbind, lapply(ir$irf, function(m) m[11, ]))  # H=10
cat("ortho IRF at h=10 (KxK):\n"); print(psi10, digits = 12)

# --- VARselect (maxlags=4) ---
vs <- VARselect(y, lag.max = 4, type = "const")
cat("== VARselect ==\n")
print(vs$selection)
cat("criteria rows (HQ SC FPE AIC):\n"); print(vs$criteria, digits = 12)

# --- trend=none 单点 ---
vn <- VAR(y, p = 2, type = "none")
cn <- sapply(vn$varresult, function(m) coef(m))
cat("== trend=n coef ==\n"); print(cn, digits = 12)
cat("n loglik:", sprintf("%.12f", logLik(vn)), "\n")
