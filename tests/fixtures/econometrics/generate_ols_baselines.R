## =============================================================================
## generate_ols_baselines.R
## Phase 6 v1.5 M1 - OLS + HC0-HC5 R baseline 生成
##
## 对标:
##   - R sandwich 3.1-3 vcovHC (MacKinnon-White 1985, Long-Ervin 2000)
##   - R datasets::longley (16 obs x 6 regressors,Employment ~ GNP.deflator + GNP
##                          + Unemployed + Armed.Forces + Population + Year)
##   - Greene 表3.x OLS系数对照
##
## 排幻觉点:
##   E1: R lm()默认含截距,C++侧不自动添加 (C++必须显式构造 [1, x] 列)
##   E2: HC1 = N/(N-K) * HC0  (自由度调整,非简单 1/N)
##   E3: HC2/HC3 leverage h_i = x_i'(X'X)^{-1}x_i  (hat values,非残差平方)
##
## 输出: tests/fixtures/econometrics/ols_hc_baseline.json
##   - longley_coefficients (6个)
##   - longley_hc0_vcov ~ longley_hc5_vcov (6x6矩阵)
##   - longley_n_obs, longley_n_params
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))

suppressMessages(library(sandwich))

cat("=== Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("sandwich version:", as.character(packageVersion("sandwich")), "\n")

## ============================================================================
## STEP 1: 加载 Longley 数据 (R datasets 包)
## ============================================================================
data("longley", package = "datasets")
cat("\n=== Longley data ===\n")
cat("N obs:", nrow(longley), "\n")
cat("Variables:", paste(colnames(longley), collapse = ", "), "\n")
cat("Head:\n")
print(head(longley, 3))

## ============================================================================
## STEP 2: OLS 回归 Employed ~ . (含截距)
## 模型: Employed = b0 + b1*GNP.deflator + b2*GNP + b3*Unemployed
##                + b4*Armed.Forces + b5*Population + b6*Year + e
## 排幻觉点 E1: R lm()默认含截距,C++必须显式构造 [1, x] 列
## ============================================================================
cat("\n=== OLS regression: Employed ~ GNP.deflator + GNP + Unemployed +",
    "Armed.Forces + Population + Year ===\n")
fm <- lm(Employed ~ GNP.deflator + GNP + Unemployed + Armed.Forces +
           Population + Year, data = longley)
coefs <- coef(fm)
cat("Coefficients:\n")
print(coefs)
n_obs <- length(fitted(fm))
n_coef <- length(coefs)
cat("N:", n_obs, "\n")
cat("K (含截距):", n_coef, "\n")

## ============================================================================
## STEP 3: 残差 + 设计矩阵 + XtX_inv (排幻觉点 E3 leverage 计算依赖)
## ============================================================================
residuals_ols <- residuals(fm)
X <- model.matrix(fm)  ## 含截距列 (排幻觉点 E1)
XtX <- t(X) %*% X
XtX_inv <- solve(XtX)
hat_values <- diag(X %*% XtX_inv %*% t(X))  ## leverage h_i (排幻觉点 E3)
cat("\n=== Design matrix ===\n")
cat("X dimensions:", dim(X), "\n")
cat("XtX_inv dimensions:", dim(XtX_inv), "\n")
cat("Leverage h_i (前5):", head(hat_values, 5), "\n")
cat("Sum h_i (应等于 K):", sum(hat_values), "\n")

## ============================================================================
## STEP 4: HC0 - White (1980) 经典异方差稳健协方差
## V_HC0 = (X'X)^{-1} X' diag(e_i^2) X (X'X)^{-1}
## ============================================================================
cat("\n=== HC0 (White 1980) vcovHC type=HC0 ===\n")
V_hc0 <- sandwich::vcovHC(fm, type = "HC0")
cat("V_hc0:\n")
print(V_hc0)
cat("SE:", sqrt(diag(V_hc0)), "\n")

## ============================================================================
## STEP 5: HC1 - 异方差稳健 + 自由度调整
## 排幻觉点 E2: HC1 = N/(N-K) * HC0  (非简单 1/N)
## ============================================================================
cat("\n=== HC1 vcovHC type=HC1 ===\n")
V_hc1 <- sandwich::vcovHC(fm, type = "HC1")
cat("V_hc1:\n")
print(V_hc1)
cat("SE:", sqrt(diag(V_hc1)), "\n")
cat("Verify HC1 = N/(N-K) * HC0:", all.equal(V_hc1, V_hc0 * n_obs / (n_obs - n_coef)),
    "\n")

## ============================================================================
## STEP 6: HC2 - 基于 leverage h_i 调整
## 排幻觉点 E3: HC2 用 e_i^2 / (1 - h_i)
## V_HC2 = (X'X)^{-1} [sum_i (e_i^2/(1-h_i)) x_i x_i'] (X'X)^{-1}
## ============================================================================
cat("\n=== HC2 vcovHC type=HC2 ===\n")
V_hc2 <- sandwich::vcovHC(fm, type = "HC2")
cat("V_hc2:\n")
print(V_hc2)
cat("SE:", sqrt(diag(V_hc2)), "\n")

## ============================================================================
## STEP 7: HC3 - jackknife 估计 (leverage 平方放大)
## 排幻觉点 E3: HC3 用 e_i^2 / (1 - h_i)^2
## ============================================================================
cat("\n=== HC3 vcovHC type=HC3 ===\n")
V_hc3 <- sandwich::vcovHC(fm, type = "HC3")
cat("V_hc3:\n")
print(V_hc3)
cat("SE:", sqrt(diag(V_hc3)), "\n")

## ============================================================================
## STEP 8: HC4 - Cribari-Neto (2004) 基于 kurtosis 调整
## ============================================================================
cat("\n=== HC4 vcovHC type=HC4 ===\n")
V_hc4 <- sandwich::vcovHC(fm, type = "HC4")
cat("V_hc4:\n")
print(V_hc4)
cat("SE:", sqrt(diag(V_hc4)), "\n")

## ============================================================================
## STEP 9: HC5 - Cribari-Neto (2007) 改进版
## ============================================================================
cat("\n=== HC5 vcovHC type=HC5 ===\n")
V_hc5 <- sandwich::vcovHC(fm, type = "HC5")
cat("V_hc5:\n")
print(V_hc5)
cat("SE:", sqrt(diag(V_hc5)), "\n")

## ============================================================================
## STEP 10: Nerlove 数据集 (145 obs) - 可选数据集
## 如 car 包可用则加载,否则跳过 (Longley 为主数据集)
## ============================================================================
cat("\n=== Nerlove data (optional) ===\n")
nerlove_available <- requireNamespace("car", quietly = TRUE)
cat("car package available:", nerlove_available, "\n")
if (nerlove_available) {
  data("Nerlove", package = "car")
  cat("Nerlove N obs:", nrow(Nerlove), "\n")
} else {
  cat("car package not available, using Longley as primary dataset\n")
}

## ============================================================================
## STEP 11: 输出 JSON
## ============================================================================
## 将矩阵转为 list (按行展开)
mat_to_vec <- function(M) {
  as.numeric(t(M))
}

baseline <- list(
  metadata = list(
    generator = "generate_ols_baselines.R",
    r_version = R.version.string,
    sandwich_version = as.character(packageVersion("sandwich")),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    spec_reference = "docs/phases/phase6/PHASE6_ECONOMETRICS_SPEC.md",
    tolerance = 1e-8,
    note = "OLS HC0-HC5 baseline: MacKinnon-White 1985, Long-Ervin 2000",
    hallucination_notes = list(
      E1 = "R lm() defaults to intercept; C++ must explicitly add [1, x] column",
      E2 = "HC1 = N/(N-K) * HC0 (degrees-of-freedom adjustment)",
      E3 = "HC2/HC3 leverage h_i = x_i'(X'X)^{-1}x_i (hat values)"
    )
  ),
  data = list(
    name = "Longley",
    source = "datasets::longley",
    n_obs = n_obs,
    n_params = n_coef,
    formula = "Employed ~ GNP.deflator + GNP + Unemployed + Armed.Forces + Population + Year"
  ),
  longley_coefficients = as.numeric(coefs),
  longley_n_obs = n_obs,
  longley_n_params = n_coef,
  longley_coef_names = names(coefs),
  ## HC0 - White (1980)
  longley_hc0_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_hc0),
    diag = as.numeric(diag(V_hc0)),
    se = as.numeric(sqrt(diag(V_hc0)))
  ),
  ## HC1 - 自由度调整 (排幻觉点 E2)
  longley_hc1_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_hc1),
    diag = as.numeric(diag(V_hc1)),
    se = as.numeric(sqrt(diag(V_hc1)))
  ),
  ## HC2 - leverage 调整 (排幻觉点 E3)
  longley_hc2_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_hc2),
    diag = as.numeric(diag(V_hc2)),
    se = as.numeric(sqrt(diag(V_hc2)))
  ),
  ## HC3 - jackknife
  longley_hc3_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_hc3),
    diag = as.numeric(diag(V_hc3)),
    se = as.numeric(sqrt(diag(V_hc3)))
  ),
  ## HC4 - Cribari-Neto 2004
  longley_hc4_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_hc4),
    diag = as.numeric(diag(V_hc4)),
    se = as.numeric(sqrt(diag(V_hc4)))
  ),
  ## HC5 - Cribari-Neto 2007
  longley_hc5_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_hc5),
    diag = as.numeric(diag(V_hc5)),
    se = as.numeric(sqrt(diag(V_hc5)))
  )
)

## 写入 JSON
## 注: digits=17 保证 double 全精度 (R jsonlite 默认 digits=4 严重精度丢失)
out_file <- "tests/fixtures/econometrics/ols_hc_baseline.json"
json_str <- jsonlite::toJSON(baseline, auto_unbox = TRUE, pretty = TRUE,
                             digits = 17)
writeLines(json_str, out_file)
cat("\n=== Baseline written to:", out_file, "===\n")
cat("File size:", file.size(out_file), "bytes\n")

## ============================================================================
## STEP 12: 自验证 - 确认 JSON 可读回
## ============================================================================
cat("\n=== Self-verification ===\n")
loaded <- jsonlite::fromJSON(out_file)
cat("Loaded coef:", loaded$longley_coefficients, "\n")
cat("Loaded HC0 diag:", loaded$longley_hc0_vcov$diag, "\n")
cat("Loaded HC5 diag:", loaded$longley_hc5_vcov$diag, "\n")
cat("All matches:",
    all.equal(loaded$longley_coefficients, as.numeric(coefs)),
    all.equal(loaded$longley_hc0_vcov$diag, as.numeric(diag(V_hc0))),
    all.equal(loaded$longley_hc5_vcov$diag, as.numeric(diag(V_hc5))), "\n")

cat("\n=== OLS baselines generated ===\n")
