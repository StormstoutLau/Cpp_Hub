## =============================================================================
## generate_hac_baselines.R
## Phase 6 v1.5 M1 - Newey-West HAC R baseline 生成
##
## 对标:
##   - R sandwich 3.1-3 NeweyWest (Newey-West 1987, Andrews 1991)
##   - R datasets::longley (16 obs x 6 regressors)
##   - 内核: Bartlett (默认), QS (Quadratic Spectral), Parzen
##
## 排幻觉点:
##   E4: 默认 max_lag=0 触发 NW 经验法则 floor(4*(T/100)^(2/9)),
##       非 Andrews (1991) 自动带宽选择 ( sandwich::bwAndrews )
##   E5: Bartlett 权重 w[l] = 1 - l/(L+1), 非 1 - l/L
##       (l = 1..L, L = max_lag; 关键差异在分母 (L+1) vs L)
##
## 输出: tests/fixtures/econometrics/hac_baseline.json
##   - longley_nw_vcov (6x6)       Newey-West Bartlett lag=4
##   - longley_nw_lag               使用的 max_lag
##   - longley_qs_vcov (6x6)        Quadratic Spectral 内核
##   - longley_parzen_vcov (6x6)    Parzen 内核
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))

suppressMessages(library(sandwich))

cat("=== Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("sandwich version:", as.character(packageVersion("sandwich")), "\n")

## ============================================================================
## STEP 1: 加载 Longley 数据
## ============================================================================
data("longley", package = "datasets")
cat("\n=== Longley data ===\n")
cat("N obs:", nrow(longley), "\n")
cat("Head:\n")
print(head(longley, 3))

## ============================================================================
## STEP 2: OLS 回归 (排幻觉点 E1: lm()默认含截距)
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
cat("K:", n_coef, "\n")

## ============================================================================
## STEP 3: Newey-West HAC, Bartlett 内核, max_lag = 4
## 排幻觉点 E4: 默认 lag=0 触发 NW 经验法则 floor(4*(T/100)^(2/9))
##              非 Andrews 自动带宽 (bwAndrews)
## 排幻觉点 E5: Bartlett 权重 w[l] = 1 - l/(L+1), l=1..L
##              非 1 - l/L (注意分母差异)
## ============================================================================
cat("\n=== Newey-West HAC (Bartlett, max_lag=4) ===\n")
nw_lag <- 4
V_nw <- sandwich::NeweyWest(fm, lag = nw_lag, prewhite = FALSE)
cat("V_nw:\n")
print(V_nw)
cat("SE:", sqrt(diag(V_nw)), "\n")

## 验证 NW 经验法则 (排幻觉点 E4)
nw_default_lag <- floor(4 * (n_obs / 100)^(2 / 9))
cat("NW default lag (floor(4*(T/100)^(2/9))):", nw_default_lag, "\n")
cat("Verify with default lag:\n")
V_nw_default <- sandwich::NeweyWest(fm, prewhite = FALSE)
cat("Default-lag V equals user-lag V?:",
    all.equal(V_nw_default, sandwich::NeweyWest(fm, lag = nw_default_lag,
                                               prewhite = FALSE)), "\n")

## ============================================================================
## STEP 4: 验证 Bartlett 权重 (排幻觉点 E5)
## w[l] = 1 - l/(L+1), l = 1..L
## ============================================================================
cat("\n=== Bartlett kernel weights (排幻觉点 E5) ===\n")
L_max <- nw_lag
bartlett_weights <- sapply(1:L_max, function(l) 1 - l / (L_max + 1))
cat("Bartlett w[l] = 1 - l/(L+1):", bartlett_weights, "\n")
cat("Wrong formula 1 - l/L:", sapply(1:L_max, function(l) 1 - l / L_max), "\n")

## ============================================================================
## STEP 5: 在 Longley 残差上注入 AR(1) 噪声 (phi=0.7) 验证自相关场景
## y_t = X_t beta + u_t,  u_t = 0.7 * u_{t-1} + v_t
## ============================================================================
cat("\n=== AR(1) injection (phi=0.7) ===\n")
set.seed(42)
resid_ols <- residuals(fm)
T_len <- length(resid_ols)
u_ar1 <- numeric(T_len)
u_ar1[1] <- resid_ols[1]
for (t in 2:T_len) {
  u_ar1[t] <- 0.7 * u_ar1[t - 1] + rnorm(1, sd = sd(resid_ols))
}
y_ar1 <- as.numeric(fitted(fm)) + u_ar1
longley_ar1 <- longley
longley_ar1$Employed <- y_ar1
fm_ar1 <- lm(Employed ~ GNP.deflator + GNP + Unemployed + Armed.Forces +
               Population + Year, data = longley_ar1)
cat("AR(1) OLS coefficients:\n")
print(coef(fm_ar1))
V_nw_ar1 <- sandwich::NeweyWest(fm_ar1, lag = nw_lag, prewhite = FALSE)
cat("V_nw (AR(1) injected):\n")
print(V_nw_ar1)
cat("SE:", sqrt(diag(V_nw_ar1)), "\n")

## ============================================================================
## STEP 6: Quadratic Spectral (QS) 内核
## sandwich::kweights 提供 QS 内核, bwAndrews 自动带宽
## ============================================================================
cat("\n=== QS kernel (Andrews 1991) ===\n")
## vcovHC 用 kernel="QS" 通过 bwAndrews + kernHAC
V_qs <- sandwich::kernHAC(fm, kernel = "Quadratic Spectral",
                          prewhite = FALSE, adjust = FALSE)
cat("V_qs:\n")
print(V_qs)
cat("SE:", sqrt(diag(V_qs)), "\n")

## ============================================================================
## STEP 7: Parzen 内核
## ============================================================================
cat("\n=== Parzen kernel ===\n")
V_parzen <- sandwich::kernHAC(fm, kernel = "Parzen",
                              prewhite = FALSE, adjust = FALSE)
cat("V_parzen:\n")
print(V_parzen)
cat("SE:", sqrt(diag(V_parzen)), "\n")

## ============================================================================
## STEP 8: 输出 JSON
## ============================================================================
## 将矩阵转为 list (按行展开)
mat_to_vec <- function(M) {
  as.numeric(t(M))
}

baseline <- list(
  metadata = list(
    generator = "generate_hac_baselines.R",
    r_version = R.version.string,
    sandwich_version = as.character(packageVersion("sandwich")),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    spec_reference = "docs/phases/phase6/PHASE6_ECONOMETRICS_SPEC.md",
    tolerance = 1e-8,
    note = "HAC baseline: Newey-West 1987, Andrews 1991",
    hallucination_notes = list(
      E4 = paste0("Default lag=0 triggers NW rule-of-thumb ",
                  "floor(4*(T/100)^(2/9)), not Andrews (1991) auto bandwidth"),
      E5 = "Bartlett weight w[l] = 1 - l/(L+1) for l=1..L, NOT 1 - l/L"
    )
  ),
  data = list(
    name = "Longley",
    source = "datasets::longley",
    n_obs = n_obs,
    n_params = n_coef,
    formula = "Employed ~ GNP.deflator + GNP + Unemployed + Armed.Forces + Population + Year"
  ),
  ## Newey-West Bartlett lag=4
  longley_nw_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_nw),
    diag = as.numeric(diag(V_nw)),
    se = as.numeric(sqrt(diag(V_nw)))
  ),
  longley_nw_lag = nw_lag,
  longley_nw_default_lag = nw_default_lag,
  longley_nw_bartlett_weights = as.numeric(bartlett_weights),
  ## AR(1) 注入后的 Newey-West (自相关验证场景)
  longley_nw_ar1_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_nw_ar1),
    diag = as.numeric(diag(V_nw_ar1)),
    se = as.numeric(sqrt(diag(V_nw_ar1)))
  ),
  longley_ar1_phi = 0.7,
  ## Quadratic Spectral 内核
  longley_qs_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_qs),
    diag = as.numeric(diag(V_qs)),
    se = as.numeric(sqrt(diag(V_qs)))
  ),
  ## Parzen 内核
  longley_parzen_vcov = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_parzen),
    diag = as.numeric(diag(V_parzen)),
    se = as.numeric(sqrt(diag(V_parzen)))
  )
)

## 写入 JSON
## 注: digits=17 保证 double 全精度 (R jsonlite 默认 digits=4 严重精度丢失)
out_file <- "tests/fixtures/econometrics/hac_baseline.json"
json_str <- jsonlite::toJSON(baseline, auto_unbox = TRUE, pretty = TRUE,
                             digits = 17)
writeLines(json_str, out_file)
cat("\n=== Baseline written to:", out_file, "===\n")
cat("File size:", file.size(out_file), "bytes\n")

## ============================================================================
## STEP 9: 自验证 - 确认 JSON 可读回
## ============================================================================
cat("\n=== Self-verification ===\n")
loaded <- jsonlite::fromJSON(out_file)
cat("Loaded NW diag:", loaded$longley_nw_vcov$diag, "\n")
cat("Loaded QS diag:", loaded$longley_qs_vcov$diag, "\n")
cat("Loaded Parzen diag:", loaded$longley_parzen_vcov$diag, "\n")
cat("All matches:",
    all.equal(loaded$longley_nw_vcov$diag, as.numeric(diag(V_nw))),
    all.equal(loaded$longley_qs_vcov$diag, as.numeric(diag(V_qs))),
    all.equal(loaded$longley_parzen_vcov$diag, as.numeric(diag(V_parzen))), "\n")

cat("\n=== HAC baselines generated ===\n")
