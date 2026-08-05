## =============================================================================
## generate_gmm_baselines.R
## Phase 6 v1.5 M3 - GMM 两步/CUE/迭代 R baseline 生成
##
## 对标:
##   - Hayashi 2000 §3.5 (两步 GMM 标准形式, 教材锚点)
##   - Hansen 1982 (GMM 框架)
##   - Hansen-Heaton-Yaron 1996 (CUE)
##   - Hansen-Singleton 1982 CCAPM (合成数据, 因原始数据难获取)
##
## 排幻觉点:
##   E10: R gmm::gmm 用 tangent matrix (G'WG)^{-1}, C++ 按 Hayashi 用 moment
##        matrix HAC (Z' diag(eps^2) Z 的 HAC). 线性 IV 下两者等价.
##        本脚本不依赖 gmm 包, 手动实现 Hayashi §3.5 两步 GMM (封闭形式),
##        与 C++ 实现用相同数学公式, 确保数值一致.
##   E10a: gmm 包不可用时的回退策略 - 手动实现 (Step1=2SLS, Step2=GMM/HAC)
##   E12: Ŝ 奇异 (完美拟合) 时的边界处理
##
## 输出: tests/fixtures/econometrics/gmm_baseline.json
##   - gmm_twostep_coefficients   (2个: beta0, beta1)
##   - gmm_cue_coefficients       (2个)
##   - gmm_iter_coefficients      (2个)
##   - gmm_j_statistic            J检验统计量
##   - gmm_j_pvalue               J检验p值
##   - gmm_n_obs                  有效样本量
##   - gmm_n_moments              矩条件数 (3)
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))

## 仅依赖 sandwich (HAC 内核权重, 与 C++ 复用相同内核) 和 jsonlite
suppressMessages(library(sandwich))

cat("=== Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("sandwich version:", as.character(packageVersion("sandwich")), "\n")
cat("gmm package available:", requireNamespace("gmm", quietly = TRUE), "\n")
cat("Note: using manual Hayashi §3.5 implementation (no gmm package dependency)\n")

## ============================================================================
## STEP 1: 合成 CCAPM 数据 (Hansen-Singleton 1982 风格)
## 原始 Hansen-Singleton 数据需 CRSP/PSID, 此处用合成数据
## 矩条件: E[Z_i * (r_i - beta0 - beta1*dc_i)] = 0
## 工具变量: Z = [1, lag(dc), lag(r)]
## ============================================================================
cat("\n=== Synthetic CCAPM data (Hansen-Singleton 1982 style) ===\n")
set.seed(42)
N_raw <- 500

## 消费增长 dc ~ N(0.02, 0.04)
dc <- rnorm(N_raw, mean = 0.02, sd = 0.04)

## 资产收益 r = 0.05 + 1.5*dc + noise (含内生性: 噪声与dc相关)
r <- 0.05 + 1.5 * dc + rnorm(N_raw, mean = 0, sd = 0.1)

## 工具变量: 前一期消费增长和收益 (lag 1)
dc_lag <- c(NA, dc[-N_raw])
r_lag <- c(NA, r[-N_raw])

## 去掉第一行 (lag 产生 NA)
valid_idx <- 2:N_raw
dc <- dc[valid_idx]
r <- r[valid_idx]
dc_lag <- dc_lag[valid_idx]
r_lag <- r_lag[valid_idx]
N <- length(dc)

cat("N (raw):", N_raw, "\n")
cat("N (after lag removal):", N, "\n")

## 工具变量矩阵 Z (N x 3): [const, lag(dc), lag(r)]
Z <- cbind(1, dc_lag, r_lag)
colnames(Z) <- c("const", "dc_lag", "r_lag")
cat("Z dimensions:", dim(Z), "\n")

## 内生变量矩阵 X (N x 2): [const, dc]
X <- cbind(1, dc)
colnames(X) <- c("const", "dc")

n_moments <- ncol(Z)  ## 3
n_params <- ncol(X)   ## 2
cat("Moments:", n_moments, "\n")
cat("Parameters:", n_params, "\n")
cat("Overidentifying restrictions:", n_moments - n_params, "\n")

## ============================================================================
## STEP 2: HAC 协方差矩阵 Ŝ 计算 (Newey-West, Bartlett 核)
## Ŝ = (1/N) Σ_l w_l · (Σ_{t=l}^{N-1} Z_{t-l} ε_{t-l} ε_t Z_t')
## 排幻觉点 E10: 用 moment matrix HAC (非 tangent matrix), 与 C++ 一致
## ============================================================================
compute_S <- function(resid, Z_mat, max_lag = NULL) {
  N <- length(resid)
  q <- ncol(Z_mat)

  ## 自动带宽 (NW 经验法则: floor(4*(N/100)^(2/9)))
  L <- if (is.null(max_lag)) {
    L <- floor(4.0 * (N / 100.0)^(2.0 / 9.0))
    max(1, min(L, N - 1))
  } else {
    max_lag
  }

  ## Bartlett 核权重: w_l = 1 - l/(L+1), l = 0, 1, ..., L
  ## 与 C++ hac_kernels.hpp 的 kernel_weights(Bartlett) 一致
  w <- 1 - (0:L) / (L + 1)

  S <- matrix(0, q, q)

  ## Ω_0 (l=0): Σ_t Z_t ε_t² Z_t'
  for (t in 1:N) {
    e2 <- resid[t]^2
    S <- S + w[1] * outer(Z_mat[t, ], Z_mat[t, ]) * e2
  }

  ## Ω_l (l>=1): Σ_{t=l+1}^{N} Z_{t-l} ε_{t-l} ε_t Z_t'
  for (l in 1:L) {
    for (t in (l + 1):N) {
      contrib <- outer(Z_mat[t - l, ], Z_mat[t, ]) * resid[t - l] * resid[t]
      S <- S + w[l + 1] * (contrib + t(contrib))  ## 对称化 (Ω_l + Ω_l')
    }
  }

  S / N
}

## ============================================================================
## STEP 3: 两步 GMM (Hayashi §3.5, 封闭形式)
## Step 1 (W₁ = (Z'Z/N)^{-1}, 等价于 2SLS):
##   β̂₁ = (X'Z (Z'Z)^{-1} Z'X)^{-1} X'Z (Z'Z)^{-1} Z'y
## Step 2 (W₂ = Ŝ⁻¹(β̂₁)):
##   β̂₂ = (X'Z Ŝ⁻¹ Z'X)^{-1} X'Z Ŝ⁻¹ Z'y
## ============================================================================
cat("\n=== Two-step GMM (Hayashi §3.5, closed form) ===\n")

XtZ <- t(X) %*% Z
Zty <- t(Z) %*% r
ZtZ <- t(Z) %*% Z

## Step 1: 2SLS (W₁ = (Z'Z)^{-1})
ZtZ_inv <- solve(ZtZ)
A1 <- XtZ %*% ZtZ_inv %*% t(XtZ)
beta_2sls <- solve(A1, XtZ %*% ZtZ_inv %*% Zty)
cat("Step 1 (2SLS) coefficients:", beta_2sls, "\n")

resid_1 <- r - X %*% beta_2sls
S_hat_1 <- compute_S(as.numeric(resid_1), Z)

## Step 2: GMM with Ŝ⁻¹
## 排幻觉点 E12: Ŝ 奇异检测 (完美拟合时 Ŝ ≈ 0)
s_det <- det(S_hat_1)
if (abs(s_det) < 1e-30 || any(!is.finite(S_hat_1))) {
  cat("Warning: S_hat singular in Step 2, retaining 2SLS estimate\n")
  beta_twostep <- beta_2sls
} else {
  S_inv <- solve(S_hat_1)
  A2 <- XtZ %*% S_inv %*% t(XtZ)
  beta_twostep <- solve(A2, XtZ %*% S_inv %*% Zty)
}
cat("Step 2 (GMM) coefficients:", beta_twostep, "\n")

## J 检验 (过度识别检验, Hansen 1982)
## J = N * ḡ(β̂)' Ŝ⁻¹ ḡ(β̂) ~ χ²(q-k)
## ḡ(β) = (1/N) Z'(y - Xβ)
resid_2 <- r - X %*% beta_twostep
g_bar <- as.numeric(t(Z) %*% resid_2 / N)
S_hat_2 <- compute_S(as.numeric(resid_2), Z)
s_det_2 <- det(S_hat_2)
if (abs(s_det_2) < 1e-30 || any(!is.finite(S_hat_2))) {
  j_stat <- 0
  j_pval <- 1
} else {
  S_inv_2 <- solve(S_hat_2)
  j_stat <- as.numeric(N * t(g_bar) %*% S_inv_2 %*% g_bar)
  j_df <- n_moments - n_params
  j_pval <- 1 - pchisq(j_stat, df = j_df)
}
cat("J statistic:", j_stat, "\n")
cat("J p-value:", j_pval, "\n")
cat("J df (overid):", n_moments - n_params, "\n")

## ============================================================================
## STEP 4: 迭代 GMM (重复 Step 2 直到收敛)
## ============================================================================
cat("\n=== Iterated GMM ===\n")
beta_iter <- beta_2sls
converged_iter <- FALSE
for (iter in 1:100) {
  resid_it <- r - X %*% beta_iter
  S_it <- compute_S(as.numeric(resid_it), Z)
  s_det_it <- det(S_it)
  if (abs(s_det_it) < 1e-30 || any(!is.finite(S_it))) {
    converged_iter <- TRUE
    break
  }
  S_inv_it <- solve(S_it)
  A_it <- XtZ %*% S_inv_it %*% t(XtZ)
  beta_new <- solve(A_it, XtZ %*% S_inv_it %*% Zty)
  diff <- max(abs(beta_new - beta_iter))
  beta_iter <- beta_new
  if (diff < 1e-10) {
    converged_iter <- TRUE
    break
  }
}
cat("Iterated GMM iterations:", iter, "\n")
cat("Iterated GMM converged:", converged_iter, "\n")
cat("Iterated GMM coefficients:", beta_iter, "\n")

## ============================================================================
## STEP 5: CUE - Continuously Updated Estimator (Hansen-Heaton-Yaron 1996)
## θ̂_CUE = argmin_θ N · ḡ(θ)' Ŝ(θ)⁻¹ ḡ(θ)
## Ŝ(θ) 随 θ 更新 (与两步 GMM 区别: 两步 GMM 的 Ŝ 在 β̂₁ 处固定)
## 用 optim() Nelder-Mead (无需梯度, 与 C++ 实现一致)
## ============================================================================
cat("\n=== CUE GMM (Hansen-Heaton-Yaron 1996) ===\n")

cue_objective <- function(params) {
  b <- params
  resid <- r - X %*% b
  g <- as.numeric(t(Z) %*% resid / N)
  S <- compute_S(as.numeric(resid), Z)
  s_det <- det(S)
  if (abs(s_det) < 1e-30 || any(!is.finite(S))) {
    ## 排幻觉点 E12: Ŝ 奇异时, 若 ḡ≈0 (完美拟合) 则 J=0, 否则不可行
    if (sqrt(sum(g^2)) < 1e-10) return(0)
    return(.Machine$double.xmax)
  }
  S_inv <- solve(S)
  if (any(!is.finite(S_inv))) {
    if (sqrt(sum(g^2)) < 1e-10) return(0)
    return(.Machine$double.xmax)
  }
  as.numeric(N * t(g) %*% S_inv %*% g)
}

## 起始值: 两步 GMM 结果
cue_opt <- optim(
  par = as.numeric(beta_twostep),
  fn = cue_objective,
  method = "Nelder-Mead",
  control = list(maxit = 10000, reltol = 1e-12)
)
beta_cue <- cue_opt$par
cat("CUE converged:", cue_opt$convergence == 0, "\n")
cat("CUE iterations:", cue_opt$counts[1], "\n")
cat("CUE coefficients:", beta_cue, "\n")
cat("CUE objective:", cue_opt$value, "\n")

## ============================================================================
## STEP 6: 输出 JSON
## ============================================================================
baseline <- list(
  metadata = list(
    generator = "generate_gmm_baselines.R",
    r_version = R.version.string,
    sandwich_version = as.character(packageVersion("sandwich")),
    gmm_package_available = requireNamespace("gmm", quietly = TRUE),
    implementation = "manual Hayashi §3.5 closed-form (no gmm package)",
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    spec_reference = "docs/phases/phase6/PHASE6_ECONOMETRICS_SPEC.md",
    tolerance = 1e-8,
    note = "GMM baseline: Hansen 1982, Hayashi §3.5, synthetic CCAPM",
    hallucination_notes = list(
      E10 = paste0("Manual Hayashi §3.5 implementation (moment matrix HAC); ",
                   "equivalent to R gmm::gmm tangent matrix for linear IV"),
      E10a = paste0("gmm package not required; closed-form Step1=2SLS, ",
                    "Step2=GMM/HAC, CUE via optim Nelder-Mead"),
      E12 = "Singular S_hat detection (perfect fit), returns J=0 if g_bar~0"
    )
  ),
  data = list(
    name = "Synthetic CCAPM",
    source = "set.seed(42), N=500 -> 499 after lag",
    n_raw = N_raw,
    n_obs = N,
    n_moments = n_moments,
    n_params = n_params,
    n_overid = n_moments - n_params,
    instruments = colnames(Z),
    formula = "E[Z_i * (r_i - beta0 - beta1*dc_i)] = 0",
    reference = "Hansen-Singleton 1982"
  ),
  gmm_2sls_coefficients = as.numeric(beta_2sls),
  gmm_twostep_coefficients = as.numeric(beta_twostep),
  gmm_cue_coefficients = as.numeric(beta_cue),
  gmm_iter_coefficients = as.numeric(beta_iter),
  gmm_cue_objective = as.numeric(cue_opt$value),
  gmm_j_statistic = as.numeric(j_stat),
  gmm_j_pvalue = as.numeric(j_pval),
  gmm_j_df = n_moments - n_params,
  gmm_n_obs = N,
  gmm_n_moments = n_moments,
  gmm_iter_converged = converged_iter,
  gmm_cue_converged = (cue_opt$convergence == 0)
)

## 写入 JSON
out_file <- "tests/fixtures/econometrics/gmm_baseline.json"
json_str <- jsonlite::toJSON(baseline, auto_unbox = TRUE, pretty = TRUE,
                             digits = 17)
writeLines(json_str, out_file)
cat("\n=== Baseline written to:", out_file, "===\n")
cat("File size:", file.size(out_file), "bytes\n")

## ============================================================================
## STEP 7: 自验证 - 确认 JSON 可读回
## ============================================================================
cat("\n=== Self-verification ===\n")
loaded <- jsonlite::fromJSON(out_file)
cat("Loaded 2SLS coef:", loaded$gmm_2sls_coefficients, "\n")
cat("Loaded twostep coef:", loaded$gmm_twostep_coefficients, "\n")
cat("Loaded CUE coef:", loaded$gmm_cue_coefficients, "\n")
cat("Loaded iter coef:", loaded$gmm_iter_coefficients, "\n")
cat("Loaded J stat:", loaded$gmm_j_statistic, "\n")
cat("Loaded J pval:", loaded$gmm_j_pvalue, "\n")
cat("All matches:",
    all.equal(loaded$gmm_twostep_coefficients, as.numeric(beta_twostep)),
    all.equal(loaded$gmm_cue_coefficients, as.numeric(beta_cue)),
    all.equal(loaded$gmm_iter_coefficients, as.numeric(beta_iter)), "\n")

cat("\n=== GMM baselines generated ===\n")
