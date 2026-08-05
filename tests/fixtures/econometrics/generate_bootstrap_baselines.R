## =============================================================================
## generate_bootstrap_baselines.R
## Phase 6 v1.5 M4 - Bootstrap CI R baseline 生成 (法学院数据)
##
## 对标:
##   - Efron-Tibshirani 1993, An Introduction to the Bootstrap
##   - Table 3.1 (law school data, 15 schools, LSAT vs GPA)
##   - Pearson 相关系数 r 的 Bootstrap CI
##
## 数据: 法学院 15 学校 (Efron-Tibshirani 1993 Table 3.1)
##   LSAT = c(576, 635, 558, 578, 666, 580, 555, 661, 651, 605,
##            653, 575, 545, 572, 594)
##   GPA  = c(3.39, 3.30, 2.81, 3.03, 3.44, 3.07, 3.00, 3.43, 3.36, 3.13,
##            3.12, 2.74, 2.76, 2.88, 2.96)
##
## 输出: tests/fixtures/econometrics/bootstrap_baseline.json
##   - law_school_pearson_r           样本 Pearson 相关系数
##   - law_school_bootstrap_ci_lower  百分位 CI 下界 (2.5%)
##   - law_school_bootstrap_ci_upper  百分位 CI 上界 (97.5%)
##   - law_school_n_obs               15
##   - law_school_B                   999
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))

suppressMessages(library(sandwich))
## boot 包提供 BCa CI (如可用)
boot_available <- requireNamespace("boot", quietly = TRUE)
if (boot_available) {
  suppressMessages(library(boot))
}

cat("=== Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("sandwich version:", as.character(packageVersion("sandwich")), "\n")
cat("boot package available:", boot_available, "\n")

## ============================================================================
## STEP 1: 法学院数据 (Efron-Tibshirani 1993, Table 3.1)
## 15 所美国法学院的 LSAT (法学院入学考试) 与 GPA 平均分
## ============================================================================
cat("\n=== Law school data (Efron-Tibshirani 1993, Table 3.1) ===\n")
LSAT <- c(576, 635, 558, 578, 666, 580, 555, 661, 651, 605,
          653, 575, 545, 572, 594)
GPA  <- c(3.39, 3.30, 2.81, 3.03, 3.44, 3.07, 3.00, 3.43, 3.36, 3.13,
          3.12, 2.74, 2.76, 2.88, 2.96)
n_obs <- length(LSAT)
cat("N obs:", n_obs, "\n")
cat("LSAT range:", min(LSAT), "-", max(LSAT), "\n")
cat("GPA range:", min(GPA), "-", max(GPA), "\n")
cat("Head (LSAT, GPA):\n")
print(head(data.frame(LSAT = LSAT, GPA = GPA), 5))

## ============================================================================
## STEP 2: 样本 Pearson 相关系数
## Efron-Tibshirani 1993 关注的相关系数 r (而非回归系数 beta)
## ============================================================================
cat("\n=== Sample Pearson correlation ===\n")
pearson_r <- cor(LSAT, GPA, method = "pearson")
cat("Pearson r:", pearson_r, "\n")

## ============================================================================
## STEP 3: 配对 Bootstrap (set.seed(42), B=999)
## 每次从 1:15 有放回抽样 15 个索引,计算 r*
## ============================================================================
cat("\n=== Paired bootstrap (B=999, seed=42) ===\n")
set.seed(42)
B <- 999
r_star <- numeric(B)
for (b in 1:B) {
  idx <- sample(1:n_obs, n_obs, replace = TRUE)
  r_star[b] <- cor(LSAT[idx], GPA[idx], method = "pearson")
}
cat("Bootstrap r* mean:", mean(r_star), "\n")
cat("Bootstrap r* sd:", sd(r_star), "\n")
cat("Bootstrap r* range:", min(r_star), "-", max(r_star), "\n")

## ============================================================================
## STEP 4: 百分位 CI (2.5%, 97.5%)
## ============================================================================
cat("\n=== Percentile CI ===\n")
ci_lower <- quantile(r_star, 0.025)
ci_upper <- quantile(r_star, 0.975)
cat("2.5%:", ci_lower, "\n")
cat("97.5%:", ci_upper, "\n")
cat("CI width:", ci_upper - ci_lower, "\n")

## ============================================================================
## STEP 5: BCa CI (如 boot 包可用)
## Efron 1987 BCa (bias-corrected and accelerated) - 比 percentil 更精确
## ============================================================================
cat("\n=== BCa CI (Efron 1987) ===\n")
bca_lower <- NA_real_
bca_upper <- NA_real_
if (boot_available) {
  ## boot::boot 需要一个返回统计量的函数
  boot_stat <- function(data, indices) {
    cor(data[indices, "LSAT"], data[indices, "GPA"], method = "pearson")
  }
  law_data <- data.frame(LSAT = LSAT, GPA = GPA)
  set.seed(42)
  boot_obj <- boot(data = law_data, statistic = boot_stat, R = B)
  cat("boot t0:", boot_obj$t0, "\n")
  cat("boot se:", sd(boot_obj$t), "\n")
  ## BCa CI
  bca_ci <- boot.ci(boot_obj, type = "bca", conf = 0.95)
  if (!is.null(bca_ci) && !is.null(bca_ci$bca)) {
    bca_lower <- bca_ci$bca[4]
    bca_upper <- bca_ci$bca[5]
    cat("BCa 2.5%:", bca_lower, "\n")
    cat("BCa 97.5%:", bca_upper, "\n")
  } else {
    cat("BCa CI not available (possibly constant statistic)\n")
  }
  ## 也输出 percentile CI 通过 boot.ci
  pct_ci <- boot.ci(boot_obj, type = "perc", conf = 0.95)
  if (!is.null(pct_ci) && !is.null(pct_ci$percent)) {
    cat("boot.ci percentile 2.5%:", pct_ci$percent[4], "\n")
    cat("boot.ci percentile 97.5%:", pct_ci$percent[5], "\n")
  }
} else {
  cat("boot package not available, skipping BCa CI\n")
}

## ============================================================================
## STEP 6: 输出 JSON
## ============================================================================
baseline <- list(
  metadata = list(
    generator = "generate_bootstrap_baselines.R",
    r_version = R.version.string,
    sandwich_version = as.character(packageVersion("sandwich")),
    boot_available = boot_available,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    spec_reference = "docs/phases/phase6/PHASE6_ECONOMETRICS_SPEC.md",
    tolerance = 1e-8,
    note = "Bootstrap CI baseline: Efron-Tibshirani 1993 Table 3.1, law school data"
  ),
  data = list(
    name = "Law school (Efron-Tibshirani 1993)",
    source = "Efron-Tibshirani 1993, Table 3.1",
    n_obs = n_obs,
    variables = c("LSAT", "GPA"),
    description = "15 US law schools, LSAT vs GPA"
  ),
  ## 样本统计量
  law_school_pearson_r = as.numeric(pearson_r),
  law_school_n_obs = n_obs,
  law_school_B = B,
  law_school_seed = 42,
  ## Bootstrap 百分位 CI
  law_school_bootstrap_ci_lower = as.numeric(ci_lower),
  law_school_bootstrap_ci_upper = as.numeric(ci_upper),
  law_school_bootstrap_r_mean = as.numeric(mean(r_star)),
  law_school_bootstrap_r_sd = as.numeric(sd(r_star)),
  ## Bootstrap r* 分布 (供 C++ 侧分布形状对比)
  law_school_bootstrap_r_star = as.numeric(r_star),
  ## BCa CI (如可用)
  law_school_bca_ci_lower = as.numeric(bca_lower),
  law_school_bca_ci_upper = as.numeric(bca_upper),
  ## 原始数据 (供 C++ 测试直接加载)
  law_school_LSAT = as.numeric(LSAT),
  law_school_GPA = as.numeric(GPA)
)

## 写入 JSON
## 注: digits=17 保证 double 全精度 (R jsonlite 默认 digits=4 严重精度丢失)
out_file <- "tests/fixtures/econometrics/bootstrap_baseline.json"
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
cat("Loaded Pearson r:", loaded$law_school_pearson_r, "\n")
cat("Loaded CI lower:", loaded$law_school_bootstrap_ci_lower, "\n")
cat("Loaded CI upper:", loaded$law_school_bootstrap_ci_upper, "\n")
cat("Loaded N obs:", loaded$law_school_n_obs, "\n")
cat("Loaded B:", loaded$law_school_B, "\n")
cat("All matches:",
    all.equal(loaded$law_school_pearson_r, as.numeric(pearson_r)),
    all.equal(loaded$law_school_bootstrap_ci_lower, as.numeric(ci_lower)),
    all.equal(loaded$law_school_bootstrap_ci_upper, as.numeric(ci_upper)), "\n")

cat("\n=== Bootstrap baselines generated ===\n")
