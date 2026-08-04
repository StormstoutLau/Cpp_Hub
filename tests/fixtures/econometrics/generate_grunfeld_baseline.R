## =============================================================================
## generate_grunfeld_baseline.R
## Phase 6 v1.5 M1 - Grunfeld Cluster SE R baseline 生成
##
## 对标:
##   - R sandwich 3.1-3 vcovCL (Liang-Zeger 1986, Cameron-Gelbach-Miller 2011)
##   - R plm 2.6-4 Grunfeld data (10 firms x 20 years, 1955-1954..1974)
##   - 排幻觉点 E6: 小样本调整 G/(G-1)*(N-1)/(N-K), R vcovCL type="HC1" 默认
##
## 输出: tests/fixtures/econometrics/grunfeld_cluster_baseline.json
##   - 系数 (4 个: intercept, inv, value, capital)
##   - 单向聚类 (firm): V_cluster_firm (4x4)
##   - 单向聚类 (year): V_cluster_year (4x4)
##   - 双向聚类 (firm + year): V_cluster_twoway (4x4)
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))

suppressMessages(library(plm))
suppressMessages(library(sandwich))

cat("=== Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("sandwich version:", as.character(packageVersion("sandwich")), "\n")
cat("plm version:", as.character(packageVersion("plm")), "\n")

## ============================================================================
## STEP 1: 加载 Grunfeld 数据 (plm 包)
## ============================================================================
data("Grunfeld", package = "plm")
cat("\n=== Grunfeld data ===\n")
cat("N obs:", nrow(Grunfeld), "\n")
cat("Firms:", length(unique(Grunfeld$firm)), "\n")
cat("Years:", length(unique(Grunfeld$year)), "\n")
cat("Year range:", min(Grunfeld$year), "-", max(Grunfeld$year), "\n")
cat("Head:\n")
print(head(Grunfeld, 3))

## ============================================================================
## STEP 2: OLS 回归 inv ~ value + capital (含截距)
## 模型: inv_i = beta0 + beta1*value + beta2*capital + epsilon
## ============================================================================
cat("\n=== OLS regression: inv ~ value + capital ===\n")
fm <- lm(inv ~ value + capital, data = Grunfeld)
coefs <- coef(fm)
cat("Coefficients:\n")
print(coefs)
cat("N:", length(fitted(fm)), "\n")
cat("K:", length(coefs), "\n")

## ============================================================================
## STEP 3: 残差 + XtX_inv
## ============================================================================
residuals_ols <- residuals(fm)
X <- model.matrix(fm)  ## 含截距列
XtX <- t(X) %*% X
XtX_inv <- solve(XtX)
cat("\nX dimensions:", dim(X), "\n")
cat("XtX_inv dimensions:", dim(XtX_inv), "\n")

## ============================================================================
## STEP 4: 单向聚类 (firm) - sandwich::vcovCL type="HC1"
## 排幻觉点 E6: G/(G-1)*(N-1)/(N-K) 小样本调整
## ============================================================================
cat("\n=== One-way cluster (firm) vcovCL type=HC1 ===\n")
V_firm <- vcovCL(fm, cluster = ~firm, type = "HC1")
cat("V_firm:\n")
print(V_firm)
cat("Diagonal (variances):", diag(V_firm), "\n")
cat("SE:", sqrt(diag(V_firm)), "\n")

## ============================================================================
## STEP 5: 单向聚类 (year) - sandwich::vcovCL type="HC1"
## ============================================================================
cat("\n=== One-way cluster (year) vcovCL type=HC1 ===\n")
V_year <- vcovCL(fm, cluster = ~year, type = "HC1")
cat("V_year:\n")
print(V_year)
cat("SE:", sqrt(diag(V_year)), "\n")

## ============================================================================
## STEP 6: 双向聚类 (firm + year) - Cameron-Gelbach-Miller 2011
## V_twoway = V(firm) + V(year) - V(firm ∩ year)
## ============================================================================
cat("\n=== Two-way cluster (firm + year) vcovCL type=HC1 ===\n")
V_twoway <- vcovCL(fm, cluster = ~firm + year, type = "HC1")
cat("V_twoway:\n")
print(V_twoway)
cat("SE:", sqrt(diag(V_twoway)), "\n")

## ============================================================================
## STEP 7: 聚类数信息
## ============================================================================
n_firms <- length(unique(Grunfeld$firm))
n_years <- length(unique(Grunfeld$year))
n_obs <- nrow(Grunfeld)
n_coef <- length(coefs)
## 双向聚类交集的 cluster 数 (firm, year) 组合数
n_firm_year <- nrow(unique(Grunfeld[, c("firm", "year")]))

cat("\n=== Cluster counts ===\n")
cat("G_firm:", n_firms, "\n")
cat("G_year:", n_years, "\n")
cat("G_firm_year:", n_firm_year, "\n")
cat("N:", n_obs, "\n")
cat("K:", n_coef, "\n")

## ============================================================================
## STEP 8: 小样本调整系数 (审计核对)
## ============================================================================
adj_firm <- (n_firms / (n_firms - 1)) * ((n_obs - 1) / (n_obs - n_coef))
adj_year <- (n_years / (n_years - 1)) * ((n_obs - 1) / (n_obs - n_coef))
cat("\n=== Small-sample adjustment ===\n")
cat("adj_firm:", adj_firm, "\n")
cat("adj_year:", adj_year, "\n")

## ============================================================================
## STEP 9: 输出 JSON
## ============================================================================
## 将矩阵转为 list (按行展开)
mat_to_vec <- function(M) {
  as.numeric(t(M))
}

baseline <- list(
  metadata = list(
    generator = "generate_grunfeld_baseline.R",
    r_version = R.version.string,
    sandwich_version = as.character(packageVersion("sandwich")),
    plm_version = as.character(packageVersion("plm")),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    spec_reference = "docs/phases/phase6/PHASE6_ECONOMETRICS_SPEC.md",
    tolerance = 1e-8,
    note = "Cluster SE baseline: Liang-Zeger 1986 + Cameron-Gelbach-Miller 2011"
  ),
  data = list(
    name = "Grunfeld",
    source = "plm::Grunfeld",
    n_obs = n_obs,
    n_firms = n_firms,
    n_years = n_years,
    year_min = min(Grunfeld$year),
    year_max = max(Grunfeld$year)
  ),
  model = list(
    formula = "inv ~ value + capital",
    k = n_coef,
    coef_names = names(coefs),
    coefficients = as.numeric(coefs)
  ),
  cluster_info = list(
    G_firm = n_firms,
    G_year = n_years,
    G_firm_year = n_firm_year,
    N = n_obs,
    K = n_coef,
    adj_firm = adj_firm,
    adj_year = adj_year
  ),
  ## vcovCL cluster=~firm type=HC1
  V_firm = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_firm),
    diag = as.numeric(diag(V_firm)),
    se = as.numeric(sqrt(diag(V_firm)))
  ),
  ## vcovCL cluster=~year type=HC1
  V_year = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_year),
    diag = as.numeric(diag(V_year)),
    se = as.numeric(sqrt(diag(V_year)))
  ),
  ## vcovCL cluster=~firm+year type=HC1 (twoway)
  V_twoway = list(
    dim = c(n_coef, n_coef),
    values = mat_to_vec(V_twoway),
    diag = as.numeric(diag(V_twoway)),
    se = as.numeric(sqrt(diag(V_twoway)))
  )
)

## 写入 JSON
## 注: digits=17 保证 double 全精度 (R jsonlite 默认 digits=4 严重精度丢失)
out_file <- "tests/fixtures/econometrics/grunfeld_cluster_baseline.json"
json_str <- jsonlite::toJSON(baseline, auto_unbox = TRUE, pretty = TRUE,
                             digits = 17)
writeLines(json_str, out_file)
cat("\n=== Baseline written to:", out_file, "===\n")
cat("File size:", file.size(out_file), "bytes\n")

## ============================================================================
## STEP 10: 自验证 - 确认 JSON 可读回
## ============================================================================
cat("\n=== Self-verification ===\n")
loaded <- jsonlite::fromJSON(out_file)
cat("Loaded coef:", loaded$model$coefficients, "\n")
cat("Loaded V_firm diag:", loaded$V_firm$diag, "\n")
cat("Loaded V_twoway diag:", loaded$V_twoway$diag, "\n")
cat("All matches:", all.equal(loaded$model$coefficients, as.numeric(coefs)),
    all.equal(loaded$V_firm$diag, as.numeric(diag(V_firm))),
    all.equal(loaded$V_twoway$diag, as.numeric(diag(V_twoway))), "\n")

cat("\n=== Done ===\n")

## ============================================================================
## STEP 11: 生成 C++ 数据头文件 (grunfeld_data.hpp)
## 避免 C++ 测试中硬编码 200 行数据, 直接 constexpr 数组
## ============================================================================
cpp_file <- "tests/fixtures/econometrics/grunfeld_data.hpp"
cpp_lines <- c(
  "// AUTO-GENERATED by generate_grunfeld_baseline.R - DO NOT EDIT",
  "// Source: plm::Grunfeld (10 firms x 20 years, 1935-1954)",
  "// Model: inv ~ value + capital (with intercept)",
  "// Used by: test_cluster_se_integration.cpp",
  "#pragma once",
  "#include <cpphub/core/types.hpp>",
  "",
  "namespace cpphub {",
  "inline namespace v1 {",
  "namespace econometrics {",
  "namespace grunfeld {",
  "",
  paste0("constexpr Size N = ", n_obs, ";"),
  paste0("constexpr Size K = ", n_coef, ";  // intercept + value + capital"),
  paste0("constexpr Index N_FIRMS = ", n_firms, ";"),
  paste0("constexpr Index N_YEARS = ", n_years, ";"),
  "",
  "/// X matrix (N x K), row-major: [intercept, value, capital]",
  "/// Generated from plm::Grunfeld, OLS model.matrix(lm(inv ~ value + capital))",
  "constexpr Real X_DATA[N * K] = {"
)

## X matrix rows
for (i in 1:n_obs) {
  row_str <- paste(sprintf("%.17g", X[i, ]), collapse = ", ")
  cpp_lines <- c(cpp_lines, paste0("  ", row_str, if (i < n_obs) "," else ""))
}
cpp_lines <- c(cpp_lines, "};", "")

## y vector
cpp_lines <- c(cpp_lines, "/// y vector (N), Grunfeld$inv")
cpp_lines <- c(cpp_lines, paste0("constexpr Real Y_DATA[N] = {"))
for (i in 1:n_obs) {
  cpp_lines <- c(cpp_lines, paste0("  ", sprintf("%.17g", Grunfeld$inv[i]),
                                   if (i < n_obs) "," else ""))
}
cpp_lines <- c(cpp_lines, "};", "")

## firm cluster ID (0-based)
cpp_lines <- c(cpp_lines, "/// firm cluster ID (0-based, N)")
cpp_lines <- c(cpp_lines, paste0("constexpr Index FIRM_ID[N] = {"))
for (i in 1:n_obs) {
  cpp_lines <- c(cpp_lines, paste0("  ", Grunfeld$firm[i] - 1L,
                                   if (i < n_obs) "," else ""))
}
cpp_lines <- c(cpp_lines, "};", "")

## year cluster ID (0-based, relative to min year)
cpp_lines <- c(cpp_lines, "/// year cluster ID (0-based, N)")
cpp_lines <- c(cpp_lines, paste0("constexpr Index YEAR_ID[N] = {"))
for (i in 1:n_obs) {
  cpp_lines <- c(cpp_lines, paste0("  ", Grunfeld$year[i] - min(Grunfeld$year),
                                   if (i < n_obs) "," else ""))
}
cpp_lines <- c(cpp_lines, "};", "")

cpp_lines <- c(cpp_lines, "", "}  // namespace grunfeld",
               "}  // namespace econometrics",
               "}  // namespace v1",
               "}  // namespace cpphub",
               "")

writeLines(cpp_lines, cpp_file)
cat("\n=== C++ data header written to:", cpp_file, "===\n")
cat("File size:", file.size(cpp_file), "bytes\n")
cat("Lines:", length(cpp_lines), "\n")
