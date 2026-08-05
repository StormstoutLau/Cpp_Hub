## =============================================================================
## generate_arellano_bond_baselines.R
## Phase 6 v1.5 M3 - Arellano-Bond 动态面板 GMM R baseline 生成
##
## 对标:
##   - R plm::pgmm (Arellano-Bond 1991, Blundell-Bond 1998)
##   - abdata (plm包, 140 UK firms 1976-1984) 或 EmplUK (备选)
##   - dynformula 构造动态面板公式
##
## 排幻觉点:
##   E11: R plm::pgmm工具变量矩阵构造, C++严格按Arellano-Bond 1991
##        R自动生成IV (lagged levels from lag 2 onwards), 可能含额外列
##        C++必须显式构造工具变量矩阵 W = [y_{t-2}, y_{t-3}, ..., y_1]
##        差分变换 (transformation="d") 消除个体效应
##
## 输出: tests/fixtures/econometrics/arellano_bond_baseline.json
##   - ab_coefficients        (alpha, beta_wage, beta_capital等)
##   - ab_ar1_statistic       AR(1)检验统计量
##   - ab_ar1_pvalue          AR(1)检验p值
##   - ab_ar2_statistic       AR(2)检验统计量
##   - ab_ar2_pvalue          AR(2)检验p值
##   - ab_sargan_statistic    Sargan/Hansen J检验统计量
##   - ab_sargan_pvalue       Sargan/Hansen J检验p值
##   - ab_n_obs, ab_n_entities, ab_n_periods
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))

suppressMessages(library(plm))

cat("=== Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("plm version:", as.character(packageVersion("plm")), "\n")

## ============================================================================
## STEP 1: 加载面板数据 (优先 abdata, 备选 EmplUK)
## abdata: Arellano-Bond 1991 原始数据 (140 UK firms, 1976-1984)
##   注: abdata 在 plm 2.6+ 可能不再内置 (改为 EmplUK 作为标准示例)
## EmplUK: plm 自带英国就业数据 (140 firms, 1976-1984, 1031 obs)
##
## 排幻觉点 E11a: data() 对不存在的数据集发 warning 而非 error,
##   tryCatch 无法捕获, 必须用 exists() 检查对象是否真正加载
## ============================================================================
cat("\n=== Panel data loading ===\n")

abdata_ok <- FALSE
tryCatch({
  data("abdata", package = "plm")
  abdata_ok <- exists("abdata", envir = .GlobalEnv)
}, error = function(e) {
  cat("abdata load error:", conditionMessage(e), "\n")
}, warning = function(w) {
  cat("abdata load warning:", conditionMessage(w), "\n")
})

if (abdata_ok) {
  panel_data <- abdata
  cat("Using abdata (Arellano-Bond 1991 original data)\n")
} else {
  ## 排幻觉点 E11a: 回退到 EmplUK (plm 标准示例数据)
  data("EmplUK", package = "plm")
  if (!exists("EmplUK", envir = .GlobalEnv)) {
    stop("Neither abdata nor EmplUK available in plm package")
  }
  panel_data <- EmplUK
  cat("abdata not available, using EmplUK (plm standard example)\n")
}

cat("N obs:", nrow(panel_data), "\n")
cat("Columns:", paste(colnames(panel_data), collapse = ", "), "\n")
cat("Head:\n")
print(head(panel_data, 3))

## 转为 pdata.frame (plm 面板数据格式)
## index: 第一个为个体 ID, 第二个为时间
idx_cols <- if (all(c("firm", "year") %in% colnames(panel_data))) {
  c("firm", "year")
} else {
  cat("WARNING: 'firm'/'year' columns not found, using first two columns as index\n")
  colnames(panel_data)[1:2]
}
cat("Index columns:", paste(idx_cols, collapse = ", "), "\n")

pdata <- pdata.frame(panel_data, index = idx_cols, drop.index = FALSE)

n_entities <- length(unique(pdata[[idx_cols[1]]]))
n_periods <- length(unique(pdata[[idx_cols[2]]]))
n_obs <- nrow(pdata)
cat("N entities:", n_entities, "\n")
cat("N periods:", n_periods, "\n")
cat("N obs:", n_obs, "\n")

## ============================================================================
## STEP 2: Arellano-Bond 动态面板 GMM 估计
## 模型: log(emp) 的动态面板
##   emp_t = alpha1 * emp_{t-1} + alpha2 * emp_{t-2}
##         + beta_w0 * wage_t + beta_w1 * wage_{t-1}
##         + beta_c0 * capital_t + beta_c1 * capital_{t-1} + u
##
## 排幻觉点 E11b: plm 2.6+ 弃用 dynformula(), 改用 multi-part formula:
##   y ~ lag(y, 1:2) + lag(x1, 0:1) + lag(x2, 0:1) | lag(y, 2:99)
##   - 左侧: 被解释变量
##   - 中间: 解释变量 (含滞后阶数)
##   - 右侧 | 后: GMM 工具变量 (lag 2 起的水平值, Arellano-Bond 1991)
##
## transformation = "d": 差分变换 (消除个体固定效应)
## effect = "twoways": 双向效应 (个体 + 时间)
## model = "twosteps": 两步 GMM (Arellano-Bond 1991 推荐)
## ============================================================================
cat("\n=== Arellano-Bond GMM estimation ===\n")
cat("Formula: log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + lag(log(capital), 0:1) | lag(log(emp), 2:99)\n")
cat("transformation='d', effect='twoways', model='twosteps'\n")

## plm 2.6+ multi-part formula (替代弃用的 dynformula)
## 等价于 dynformula(log(emp) ~ log(wage) + log(capital), list(2, 1, 1))
form <- log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) +
  lag(log(capital), 0:1) | lag(log(emp), 2:99)
cat("Formula:\n")
print(form)

## pgmm 估计
fm_ab <- tryCatch({
  pgmm(form,
        data = pdata,
        effect = "twoways",
        model = "twosteps",
        transformation = "d")
}, error = function(e) {
  cat("ERROR in pgmm (twosteps):", conditionMessage(e), "\n")
  cat("Trying onestep model...\n")
  tryCatch({
    pgmm(form,
          data = pdata,
          effect = "twoways",
          model = "onestep",
          transformation = "d")
  }, error = function(e2) {
    cat("ERROR in pgmm (onestep):", conditionMessage(e2), "\n")
    NULL
  })
})

## ============================================================================
## STEP 3: 提取估计结果
## ============================================================================
if (!is.null(fm_ab)) {
  cat("\n=== Arellano-Bond results ===\n")
  coefs_ab <- coef(fm_ab)
  cat("Coefficients:\n")
  print(coefs_ab)
  cat("Model:", attr(fm_ab, "model"), "\n")
  cat("Transformation:", attr(fm_ab, "transformation"), "\n")

  ## ============================================================================
  ## STEP 4: AR(1)/AR(2) 检验 + Sargan/Hansen J 检验
  ## AR(1): 一阶自相关 (差分残差应有一阶自相关, 因 MA(1) 结构)
  ## AR(2): 二阶自相关 (不应有, 若有则工具变量无效)
  ## Sargan/Hansen J: 过度识别检验 (矩条件有效性)
  ##
  ## 排幻觉点 E11c: plm 2.6+ summary.pgmm 结构变更:
  ##   - 旧版 (<2.6): s$AR 为 2x2 矩阵, s$sargan 为 c(stat, pval) 向量
  ##   - 新版 (>=2.6): s$m1/s$m2 为 htest 对象 (AR1/AR2),
  ##                   s$sargan 为 htest 对象, s$wald.coef/s$wald.td 为 htest
  ## ============================================================================
  cat("\n=== Specification tests ===\n")
  s_ab <- summary(fm_ab)

  ## AR 检验提取 (兼容 plm 2.6+ htest 和旧版矩阵)
  ar1_stat <- NA
  ar1_pval <- NA
  ar2_stat <- NA
  ar2_pval <- NA
  tryCatch({
    if (!is.null(s_ab$m1)) {
      ## plm 2.6+ htest 对象
      ar1_stat <- as.numeric(s_ab$m1$statistic)
      ar1_pval <- as.numeric(s_ab$m1$p.value)
    } else if (!is.null(s_ab$AR)) {
      ## 旧版 2x2 矩阵
      ar1_stat <- as.numeric(s_ab$AR[1, 1])
      ar1_pval <- as.numeric(s_ab$AR[1, 2])
    }
    if (!is.null(s_ab$m2)) {
      ar2_stat <- as.numeric(s_ab$m2$statistic)
      ar2_pval <- as.numeric(s_ab$m2$p.value)
    } else if (!is.null(s_ab$AR)) {
      ar2_stat <- as.numeric(s_ab$AR[2, 1])
      ar2_pval <- as.numeric(s_ab$AR[2, 2])
    }
  }, error = function(e) {
    cat("AR test extraction error:", conditionMessage(e), "\n")
  })
  cat("AR(1) statistic:", ar1_stat, ", p-value:", ar1_pval, "\n")
  cat("AR(2) statistic:", ar2_stat, ", p-value:", ar2_pval, "\n")

  ## Sargan/Hansen J 检验 (plm 2.6+ htest 对象)
  sargan_stat <- NA
  sargan_pval <- NA
  tryCatch({
    if (!is.null(s_ab$sargan)) {
      if (is.list(s_ab$sargan)) {
        ## plm 2.6+ htest 对象
        sargan_stat <- as.numeric(s_ab$sargan$statistic)
        sargan_pval <- as.numeric(s_ab$sargan$p.value)
      } else {
        ## 旧版向量
        sargan_stat <- as.numeric(s_ab$sargan[1])
        sargan_pval <- as.numeric(s_ab$sargan[length(s_ab$sargan)])
      }
    }
  }, error = function(e) {
    cat("Sargan test extraction error:", conditionMessage(e), "\n")
  })
  cat("Sargan statistic:", sargan_stat, ", p-value:", sargan_pval, "\n")

  ## Wald 检验 (plm 2.6+: s$wald.coef 系数联合, s$wald.td 时间虚拟变量)
  wald_coef_stat <- NA
  wald_coef_pval <- NA
  tryCatch({
    if (!is.null(s_ab$wald.coef)) {
      wald_coef_stat <- as.numeric(s_ab$wald.coef$statistic)
      wald_coef_pval <- as.numeric(s_ab$wald.coef$p.value)
      cat("Wald (coef) statistic:", wald_coef_stat,
          ", p-value:", wald_coef_pval, "\n")
    } else if (!is.null(s_ab$Wald)) {
      cat("Wald (legacy):\n")
      print(s_ab$Wald)
    }
  }, error = function(e) {
    cat("Wald test not available:", conditionMessage(e), "\n")
  })

} else {
  cat("ERROR: pgmm estimation failed, using NA for all results\n")
  coefs_ab <- rep(NA, 4)
  ar1_stat <- ar1_pval <- ar2_stat <- ar2_pval <- NA
  sargan_stat <- sargan_pval <- NA
}

## ============================================================================
## STEP 5: 输出 JSON
## ============================================================================
baseline <- list(
  metadata = list(
    generator = "generate_arellano_bond_baselines.R",
    r_version = R.version.string,
    plm_version = as.character(packageVersion("plm")),
    data_source = if (abdata_ok) "abdata" else "EmplUK",
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    spec_reference = "docs/phases/phase6/PHASE6_ECONOMETRICS_SPEC.md",
    tolerance = 1e-8,
    note = "Arellano-Bond 1991 dynamic panel GMM baseline",
    hallucination_notes = list(
      E11 = paste0("R plm::pgmm auto-generates IV matrix ",
                   "(lagged levels from lag 2); C++ must explicitly ",
                   "construct per Arellano-Bond 1991; diff transformation ",
                   "eliminates individual effects")
    )
  ),
  data = list(
    name = if (abdata_ok) "abdata" else "EmplUK",
    source = paste0("plm::", if (abdata_ok) "abdata" else "EmplUK"),
    n_obs = n_obs,
    n_entities = n_entities,
    n_periods = n_periods,
    formula = "dynformula(log(emp) ~ log(wage) + log(capital), list(2,1,1))",
    transformation = "d",
    effect = "twoways",
    model = "twosteps",
    reference = "Arellano-Bond 1991"
  ),
  ab_coefficients = as.numeric(coefs_ab),
  ab_coef_names = names(coefs_ab),
  ab_ar1_statistic = ar1_stat,
  ab_ar1_pvalue = ar1_pval,
  ab_ar2_statistic = ar2_stat,
  ab_ar2_pvalue = ar2_pval,
  ab_sargan_statistic = sargan_stat,
  ab_sargan_pvalue = sargan_pval,
  ab_n_obs = n_obs,
  ab_n_entities = n_entities,
  ab_n_periods = n_periods
)

## 写入 JSON
out_file <- "tests/fixtures/econometrics/arellano_bond_baseline.json"
json_str <- jsonlite::toJSON(baseline, auto_unbox = TRUE, pretty = TRUE,
                             digits = 17)
writeLines(json_str, out_file)
cat("\n=== Baseline written to:", out_file, "===\n")
cat("File size:", file.size(out_file), "bytes\n")

## ============================================================================
## STEP 6: 自验证 - 确认 JSON 可读回
## ============================================================================
cat("\n=== Self-verification ===\n")
loaded <- jsonlite::fromJSON(out_file)
cat("Loaded coef:", loaded$ab_coefficients, "\n")
cat("Loaded AR1 stat:", loaded$ab_ar1_statistic, "\n")
cat("Loaded AR2 stat:", loaded$ab_ar2_statistic, "\n")
cat("Loaded Sargan stat:", loaded$ab_sargan_statistic, "\n")
cat("Loaded n_obs:", loaded$ab_n_obs, "\n")
if (!is.null(fm_ab)) {
  cat("Coefficients match:",
      all.equal(loaded$ab_coefficients, as.numeric(coefs_ab)), "\n")
}

cat("\n=== Arellano-Bond baselines generated ===\n")
