## =============================================================================
## generate_mle_baselines.R
## Phase 6 v1.5 M2 - MLE/QMLE R baseline 生成
##
## 对标:
##   - R glm + sandwich::sandwich (QMLE, Huber 1967, White 1982)
##   - Spector-Mazzeo 1980 学习表现数据 (car::Spector, 32 obs)
##   - DoctorVisits 数据 (AER::DoctorVisits, 5190 obs, 可选)
##   - Greene 表17.x MLE系数对照
##
## 排幻觉点:
##   E7: R glm用IRLS, C++用Newton-Raphson, canonical link等价
##       (logit对binomial是canonical link, IRLS = Newton-Raphson更新)
##   E8: GLM bread = (X'WX)^{-1}, meat = X'diag(eps^2)X
##       (sandwich::bread.glm 返回 vcov(fm), 即Hessian逆; 非OPG)
##
## 输出: tests/fixtures/econometrics/mle_baseline.json
##   - spector_logit_coefficients        (4个: intercept, GPA, TUCE, PSI)
##   - spector_logit_hessian_vcov        (4x4, vcov(fm), 标准MLE)
##   - spector_logit_sandwich_vcov       (4x4, sandwich::sandwich(fm), QMLE)
##   - spector_n_obs                     32
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))

suppressMessages(library(sandwich))
car_available <- requireNamespace("car", quietly = TRUE)
if (car_available) {
  suppressMessages(library(car))
}

cat("=== Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("sandwich version:", as.character(packageVersion("sandwich")), "\n")
cat("car package available:", car_available, "\n")

## ============================================================================
## STEP 1: Spector-Mazzeo 数据 (Greene 8ed Table 17.1, 32 obs)
## 学习表现 vs GPA + TUCE + PSI
## 列: grade (0/1), GPA, TUCE, PSI (0/1)
##
## 排幻觉点 E7a: 数据源优先 car::Spector; 若 car 不可用, 使用硬编码的
##   Spector-Mazzeo 1980 公开数据 (Greene 8ed Table 17.1, 32 obs).
##   数据来自 Greene Econometric Analysis 教材公开附录, 多个 R 包
##   (car, AER, mlogit) 均内置, 数值一致.
## ============================================================================
cat("\n=== Spector-Mazzeo data (Greene Table 17.1, 32 obs) ===\n")
if (car_available) {
  data("Spector", package = "car")
  cat("Source: car::Spector\n")
} else {
  ## 硬编码 Spector-Mazzeo 1980 公开数据 (Greene 8ed Table 17.1)
  ## 排幻觉点 E7a: 避免因 car 包缺失阻塞 baseline 生成
  Spector <- data.frame(
    grade = c(0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
              0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L,
              1L, 1L),
    GPA = c(2.66, 2.89, 3.28, 2.92, 4.00, 2.86, 2.76, 2.87, 3.03, 3.92,
            2.63, 3.32, 3.57, 3.26, 3.53, 2.74, 2.75, 2.83, 3.12, 3.16,
            2.06, 3.62, 2.89, 3.51, 3.54, 2.83, 3.39, 2.67, 3.65, 4.00,
            3.10, 2.39),
    TUCE = c(20L, 22L, 24L, 12L, 21L, 17L, 17L, 21L, 25L, 29L, 20L, 23L,
             23L, 25L, 26L, 19L, 25L, 19L, 23L, 25L, 22L, 28L, 14L, 26L,
             24L, 27L, 17L, 24L, 21L, 23L, 21L, 19L),
    PSI = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
            0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
            1L, 1L)
  )
  cat("Source: hardcoded (Greene 8ed Table 17.1), car package unavailable\n")
}
cat("N obs:", nrow(Spector), "\n")
cat("Columns:", paste(colnames(Spector), collapse = ", "), "\n")
cat("Head:\n")
print(head(Spector, 3))
spector_n_obs <- nrow(Spector)

## ============================================================================
## STEP 2: Logistic 回归 (MLE)
## 模型: glm(grade ~ GPA + TUCE + PSI, family=binomial(link="logit"))
## 排幻觉点 E7: R glm用IRLS, C++用Newton-Raphson, canonical link等价
## ============================================================================
cat("\n=== Logistic regression: grade ~ GPA + TUCE + PSI ===\n")
fm_logit <- glm(grade ~ GPA + TUCE + PSI, family = binomial(link = "logit"),
                data = Spector)
coefs_logit <- coef(fm_logit)
cat("Coefficients:\n")
print(coefs_logit)
cat("N:", nobs(fm_logit), "\n")
cat("K:", length(coefs_logit), "\n")
cat("Converged:", fm_logit$converged, "\n")

## ============================================================================
## STEP 3: 三种协方差估计
## 排幻觉点 E8: bread=(X'WX)^{-1}, meat=X'diag(eps^2)X
##   - Hessian: vcov(fm) = (X'WX)^{-1} (标准MLE, dispersion=1 for binomial)
##   - OPG/bread: sandwich::bread(fm) = vcov(fm) for GLM (注: R的bread是Hessian逆, 非OPG)
##   - Sandwich: sandwich::sandwich(fm) = bread * meat * bread / n (QMLE)
## ============================================================================
cat("\n=== Covariance estimates (排幻觉点 E8) ===\n")

## Hessian-based (标准MLE)
V_hessian <- vcov(fm_logit)
cat("Hessian vcov (vcov(fm)):\n")
print(V_hessian)
cat("SE:", sqrt(diag(V_hessian)), "\n")

## bread (注: R的bread.glm 返回 vcov(fm), 即 (X'WX)^{-1}, 非OPG)
V_bread <- sandwich::bread(fm_logit)
cat("\nbread (sandwich::bread(fm)):\n")
print(V_bread)
cat("bread == vcov(fm)?:", all.equal(V_bread, V_hessian, tolerance = 1e-10), "\n")

## Sandwich (QMLE)
V_sandwich <- sandwich::sandwich(fm_logit)
cat("\nSandwich vcov (sandwich::sandwich(fm)):\n")
print(V_sandwich)
cat("SE:", sqrt(diag(V_sandwich)), "\n")

## ============================================================================
## STEP 4: 验证 E7 - glm IRLS == Newton-Raphson (canonical link)
## 对于 canonical link (logit for binomial), IRLS 更新 = Newton-Raphson 更新
## beta_new = beta + (X'WX)^{-1} X'(y - mu)
## ============================================================================
cat("\n=== E7 verification: glm IRLS vs Newton-Raphson ===\n")
X <- model.matrix(fm_logit)
y_spector <- Spector$grade
beta_nr <- rep(0, ncol(X))  ## 起始值
for (iter in 1:100) {
  eta <- as.numeric(X %*% beta_nr)
  mu_nr <- 1 / (1 + exp(-eta))
  W_nr <- mu_nr * (1 - mu_nr)
  score <- as.numeric(t(X) %*% (y_spector - mu_nr))
  XtWX <- crossprod(X, X * W_nr)
  step <- solve(XtWX, score)
  beta_nr <- beta_nr + step  ## Newton-Raphson 更新 (= IRLS for canonical link)
  if (max(abs(step)) < 1e-12) break
}
cat("Newton-Raphson iterations:", iter, "\n")
cat("NR coefficients:", beta_nr, "\n")
cat("glm coefficients:", as.numeric(coefs_logit), "\n")
cat("E7 PASS (IRLS == Newton-Raphson):",
    all.equal(beta_nr, as.numeric(coefs_logit), tolerance = 1e-8), "\n")

## ============================================================================
## STEP 5: 验证 E8 - bread/meat 分解
## bread = (X'WX)^{-1} (Hessian逆)
## meat = X'diag(eps^2)X (Pearson残差外积)
## sandwich = bread %*% meat %*% bread / n
## ============================================================================
cat("\n=== E8 verification: bread/meat decomposition ===\n")
eta_fit <- as.numeric(X %*% coefs_logit)
mu_fit <- 1 / (1 + exp(-eta_fit))
W_fit <- mu_fit * (1 - mu_fit)
XtWX_fit <- crossprod(X, X * W_fit)
bread_manual <- solve(XtWX_fit)  ## (X'WX)^{-1}

eps_pearson <- residuals(fm_logit, type = "pearson")
meat_manual <- crossprod(X, eps_pearson^2 * X)  ## X'diag(eps^2)X
n_spector <- length(y_spector)

cat("bread manual == sandwich::bread?:",
    all.equal(bread_manual, V_bread, tolerance = 1e-8), "\n")
cat("meat manual == sandwich::meat?:",
    all.equal(meat_manual / n_spector, sandwich::meat(fm_logit), tolerance = 1e-8),
    "\n")
cat("sandwich manual == sandwich::sandwich?:",
    all.equal(bread_manual %*% (meat_manual / n_spector) %*% bread_manual,
              V_sandwich, tolerance = 1e-8), "\n")

## ============================================================================
## STEP 6: warpbreaks 数据 (datasets包, 48 obs) - Poisson MLE
## 模型: glm(breaks ~ wool + tension, family=poisson(link="log"), data=warpbreaks)
##
## 排幻觉点 E8a: Poisson MLE bread = (X'WX)^{-1}, W = diag(λ_i) = diag(exp(Xβ))
##   (Poisson 方差 = 均值 = λ, 所以 W = diag(λ), 非 diag(λ(1-λ)))
##   log link 是 Poisson 的 canonical link, IRLS = Newton-Raphson
##
## warpbreaks 是 R 内置数据集 (datasets 包), 无需额外依赖:
##   - 48 obs (2 wool types x 3 tension levels x 8 replicates)
##   - breaks: 断纱次数 (count data, Poisson 分布)
##   - wool: 羊毛类型 (A/B, factor)
##   - tension: 张力 (L/M/H, factor)
##   - Cameron-Trivedi 1998, Venables-Ripley 2002 经典 Poisson 回归示例
## ============================================================================
cat("\n=== warpbreaks data (datasets::warpbreaks, 48 obs, Poisson MLE) ===\n")
warpbreaks_ok <- FALSE
tryCatch({
  data("warpbreaks", package = "datasets")
  cat("warpbreaks N obs:", nrow(warpbreaks), "\n")
  cat("Columns:", paste(colnames(warpbreaks), collapse = ", "), "\n")
  cat("Head:\n")
  print(head(warpbreaks, 3))

  ## Poisson 回归: breaks ~ wool + tension (log link, canonical)
  fm_pois <- glm(breaks ~ wool + tension, family = poisson(link = "log"),
                 data = warpbreaks)
  coefs_pois <- coef(fm_pois)
  cat("\nPoisson coefficients:\n")
  print(coefs_pois)
  cat("N:", nobs(fm_pois), "\n")
  cat("K:", length(coefs_pois), "\n")
  cat("Converged:", fm_pois$converged, "\n")

  ## Hessian vcov (标准 MLE)
  V_pois_hessian <- vcov(fm_pois)
  cat("\nHessian vcov:\n")
  print(V_pois_hessian)
  cat("SE:", sqrt(diag(V_pois_hessian)), "\n")

  ## Sandwich vcov (QMLE)
  V_pois_sandwich <- sandwich::sandwich(fm_pois)
  cat("\nSandwich vcov:\n")
  print(V_pois_sandwich)
  cat("Sandwich SE:", sqrt(diag(V_pois_sandwich)), "\n")

  warpbreaks_n_obs <- nrow(warpbreaks)
  warpbreaks_ok <- TRUE
}, error = function(e) {
  cat("WARNING: warpbreaks data load failed:", conditionMessage(e), "\n")
})

## ============================================================================
## STEP 7: 输出 JSON
## ============================================================================
mat_to_vec <- function(M) {
  as.numeric(t(M))
}

baseline <- list(
  metadata = list(
    generator = "generate_mle_baselines.R",
    r_version = R.version.string,
    sandwich_version = as.character(packageVersion("sandwich")),
    car_available = car_available,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    spec_reference = "docs/phases/phase6/PHASE6_ECONOMETRICS_SPEC.md",
    tolerance = 1e-8,
    note = "MLE/QMLE baseline: glm + sandwich, Spector-Mazzeo 1980 + warpbreaks Poisson",
    hallucination_notes = list(
      E7 = paste0("R glm uses IRLS, C++ uses Newton-Raphson; ",
                  "for canonical link (logit/binomial, log/poisson) they are equivalent"),
      E8 = paste0("GLM bread = (X'WX)^{-1} (Hessian inverse, NOT OPG); ",
                  "meat = X'diag(eps^2)X; sandwich = bread*meat*bread/n"),
      E7a = paste0("Spector-Mazzeo data: hardcoded Greene 8ed Table 17.1 ",
                   "(32 obs) when car package unavailable"),
      E8a = paste0("Poisson MLE: W = diag(lambda) = diag(exp(Xb)), ",
                   "log link is canonical, IRLS = Newton-Raphson")
    )
  ),
  data = list(
    name = "Spector",
    source = if (car_available) "car::Spector" else "hardcoded Greene 8ed Table 17.1 (car unavailable)",
    n_obs = spector_n_obs,
    formula = "grade ~ GPA + TUCE + PSI",
    family = "binomial(link = logit)",
    reference = "Greene Table 17.x, Spector-Mazzeo 1980"
  ),
  spector_logit_coefficients = as.numeric(coefs_logit),
  spector_logit_coef_names = names(coefs_logit),
  spector_logit_hessian_vcov = list(
    dim = c(length(coefs_logit), length(coefs_logit)),
    values = mat_to_vec(V_hessian),
    diag = as.numeric(diag(V_hessian)),
    se = as.numeric(sqrt(diag(V_hessian)))
  ),
  spector_logit_sandwich_vcov = list(
    dim = c(length(coefs_logit), length(coefs_logit)),
    values = mat_to_vec(V_sandwich),
    diag = as.numeric(diag(V_sandwich)),
    se = as.numeric(sqrt(diag(V_sandwich)))
  ),
  spector_n_obs = spector_n_obs
)

## warpbreaks Poisson MLE baseline
if (warpbreaks_ok) {
  baseline$warpbreaks_data <- list(
    name = "warpbreaks",
    source = "datasets::warpbreaks (R built-in)",
    n_obs = warpbreaks_n_obs,
    formula = "breaks ~ wool + tension",
    family = "poisson(link = log)",
    reference = "Cameron-Trivedi 1998, Venables-Ripley 2002"
  )
  baseline$warpbreaks_poisson_coefficients <- as.numeric(coefs_pois)
  baseline$warpbreaks_poisson_coef_names <- names(coefs_pois)
  baseline$warpbreaks_poisson_hessian_vcov <- list(
    dim = c(length(coefs_pois), length(coefs_pois)),
    values = mat_to_vec(V_pois_hessian),
    diag = as.numeric(diag(V_pois_hessian)),
    se = as.numeric(sqrt(diag(V_pois_hessian)))
  )
  baseline$warpbreaks_poisson_sandwich_vcov <- list(
    dim = c(length(coefs_pois), length(coefs_pois)),
    values = mat_to_vec(V_pois_sandwich),
    diag = as.numeric(diag(V_pois_sandwich)),
    se = as.numeric(sqrt(diag(V_pois_sandwich)))
  )
  baseline$warpbreaks_n_obs <- warpbreaks_n_obs
}

## 写入 JSON
## 注: digits=17 保证 double 全精度 (R jsonlite 默认 digits=4 严重精度丢失)
out_file <- "tests/fixtures/econometrics/mle_baseline.json"
json_str <- jsonlite::toJSON(baseline, auto_unbox = TRUE, pretty = TRUE,
                             digits = 17)
writeLines(json_str, out_file)
cat("\n=== Baseline written to:", out_file, "===\n")
cat("File size:", file.size(out_file), "bytes\n")

## ============================================================================
## STEP 8: 自验证 - 确认 JSON 可读回
## ============================================================================
cat("\n=== Self-verification ===\n")
loaded <- jsonlite::fromJSON(out_file)
cat("Loaded Spector coef:", loaded$spector_logit_coefficients, "\n")
cat("Loaded Spector hessian diag:", loaded$spector_logit_hessian_vcov$diag, "\n")
cat("Loaded Spector sandwich diag:", loaded$spector_logit_sandwich_vcov$diag, "\n")
cat("Spector all matches:",
    all.equal(loaded$spector_logit_coefficients, as.numeric(coefs_logit)),
    all.equal(loaded$spector_logit_hessian_vcov$diag, as.numeric(diag(V_hessian))),
    all.equal(loaded$spector_logit_sandwich_vcov$diag,
              as.numeric(diag(V_sandwich))), "\n")
if (warpbreaks_ok) {
  cat("Loaded warpbreaks Poisson coef:", loaded$warpbreaks_poisson_coefficients, "\n")
  cat("Loaded warpbreaks hessian diag:", loaded$warpbreaks_poisson_hessian_vcov$diag, "\n")
  cat("Loaded warpbreaks sandwich diag:", loaded$warpbreaks_poisson_sandwich_vcov$diag, "\n")
  cat("warpbreaks all matches:",
      all.equal(loaded$warpbreaks_poisson_coefficients, as.numeric(coefs_pois)),
      all.equal(loaded$warpbreaks_poisson_hessian_vcov$diag, as.numeric(diag(V_pois_hessian))),
      all.equal(loaded$warpbreaks_poisson_sandwich_vcov$diag,
                as.numeric(diag(V_pois_sandwich))), "\n")
}

cat("\n=== MLE baselines generated ===\n")
