## =============================================================================
## generate_r_baselines.R
## Phase 5 v1.4.0 - HFE R baseline generator
##
## SOURCE: spec docs/phases/phase5/PHASE5_HFE_SPEC.md §7.1
## R pkg:  highfrequency 1.0.3 (Boudt, Kleen, Sjoerup 2022, JSS doi:10.18637/jss.v104.i08)
## Output: tests/fixtures/hfe/baselines.json (CI gate single source of truth)
##
## CRITICAL (2026-08-02 实测): Rscript 非交互模式不自动加载用户库路径,
##          必须显式 .libPaths() 否则 library(highfrequency) 失败.
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages({
  library(highfrequency)
  library(jsonlite)
})

## ---- metadata ----
meta <- list(
  generator      = "generate_r_baselines.R",
  r_version      = R.version.string,
  hf_version     = as.character(packageVersion("highfrequency")),
  generated_at   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  spec_reference = "docs/phases/phase5/PHASE5_HFE_SPEC.md",
  tolerance_default = 1e-12,
  tolerance_standard = 1e-10
)

## ---- Test Case 1: Constant price series (zero returns, edge case) ----
## make_returns of constant prices -> all zero -> RV=0, BPV=0, etc.
prices_const <- rep(100.0, 10)
ret_const   <- makeReturns(prices_const)         # length n-1 = 9, all 0
rv_const    <- rRVar(ret_const)
rvol_const  <- sqrt(rv_const)
rq_const    <- rQuar(ret_const)
bpv_const   <- rBPCov(ret_const, makeReturns = FALSE)  # ret_const already returns
rsv_const   <- rSV(ret_const)

## ---- Test Case 2: Simple known returns (hand-verifiable) ----
## returns = c(0.01, -0.02, 0.03, -0.01, 0.02)
ret_known <- c(0.01, -0.02, 0.03, -0.01, 0.02)
rv_known    <- rRVar(ret_known)
rvol_known  <- sqrt(rv_known)
rq_known    <- rQuar(ret_known)
bpv_known   <- rBPCov(ret_known, makeReturns = FALSE)
rsv_known   <- rSV(ret_known)

## ---- Test Case 3: GBM simulated (seed=42, n=100) ----
set.seed(42)
n_gbm <- 100
prices_gbm <- 100 * exp(cumsum(rnorm(n_gbm, 0, 0.01)))
ret_gbm    <- makeReturns(prices_gbm)
rv_gbm     <- rRVar(ret_gbm)
rvol_gbm   <- sqrt(rv_gbm)
rq_gbm     <- rQuar(ret_gbm)
bpv_gbm    <- rBPCov(ret_gbm, makeReturns = FALSE)
rsv_gbm    <- rSVar(ret_gbm)

## ---- Test Case 4: GBM + jump (for BNS rejection) ----
set.seed(123)
n_jmp <- 200
prices_jmp <- 100 * exp(cumsum(rnorm(n_jmp, 0, 0.005)))
## inject a jump at midpoint
prices_jmp[100] <- prices_jmp[99] * 1.10  # +10% jump
ret_jmp <- makeReturns(prices_jmp)
rv_jmp     <- rRVar(ret_jmp)
bpv_jmp    <- rBPCov(ret_jmp, makeReturns = FALSE)
rq_jmp     <- rQuar(ret_jmp)
rsv_jmp    <- rSV(ret_jmp)

## ---- Test Case 5: BNS jump test (no jump, should NOT reject) ----
bns_nojump <- BNSjumpTest(ret_gbm, IVestimator = "BV",
                          IQestimator = "TP", makeReturns = FALSE,
                          alpha = 0.975)

## ---- Test Case 6: BNS jump test (with jump, should reject) ----
bns_jump <- BNSjumpTest(ret_jmp, IVestimator = "BV",
                        IQestimator = "TP", makeReturns = FALSE,
                        alpha = 0.975)

## ---- Test Case 7: Multi-asset rCov (2 assets) ----
set.seed(7)
n_multi <- 50
p1 <- 100 * exp(cumsum(rnorm(n_multi, 0, 0.01)))
p2 <- 100 * exp(cumsum(rnorm(n_multi, 0, 0.012)))
ret_multi <- cbind(makeReturns(p1), makeReturns(p2))
rcov_multi <- rCov(ret_multi, makeReturns = FALSE)

## ---- Test Case 8: aggregatePrice (highfrequency sampleTData) ----
data(sampleTData)
agg_30s <- aggregatePrice(sampleTData, alignBy = "seconds", alignPeriod = 30)
n_agg   <- nrow(agg_30s)
## aggregatePrice returns data.table; PRICE column is numeric
prices_agg <- as.numeric(agg_30s$PRICE)
first_price_agg <- prices_agg[1]
last_price_agg  <- prices_agg[n_agg]

## ---- Test Case 9: makeReturns on aggregated data ----
ret_agg <- makeReturns(prices_agg)
rv_agg  <- rRVar(ret_agg)

## =============================================================================
## Assemble baseline list
## =============================================================================
baselines <- list(
  metadata = meta,
  case1_constant = list(
    description = "Constant prices -> zero returns (edge case)",
    prices = prices_const,
    n_returns = length(ret_const),
    rv  = rv_const,
    rvol = rvol_const,
    rq  = rq_const,
    bpv = bpv_const,
    rsv_pos = rsv_const$rSVarupside,
    rsv_neg = rsv_const$rSVardownside
  ),
  case2_known = list(
    description = "Hand-verifiable returns c(0.01,-0.02,0.03,-0.01,0.02)",
    returns = ret_known,
    n_returns = length(ret_known),
    rv  = rv_known,
    rvol = rvol_known,
    rq  = rq_known,
    bpv = bpv_known,
    rsv_pos = rsv_known$rSVarupside,
    rsv_neg = rsv_known$rSVardownside
  ),
  case3_gbm = list(
    description = "GBM simulated seed=42 n=100 sigma=0.01",
    prices = prices_gbm,
    n_returns = length(ret_gbm),
    rv  = rv_gbm,
    rvol = rvol_gbm,
    rq  = rq_gbm,
    bpv = bpv_gbm,
    rsv_pos = rsv_gbm$rSVarupside,
    rsv_neg = rsv_gbm$rSVardownside
  ),
  case4_jump = list(
    description = "GBM + 10% jump at midpoint, seed=123 n=200",
    prices = prices_jmp,
    n_returns = length(ret_jmp),
    rv  = rv_jmp,
    bpv = bpv_jmp,
    rq  = rq_jmp,
    rsv_pos = rsv_jmp$rSVarupside,
    rsv_neg = rsv_jmp$rSVardownside
  ),
  case5_bns_nojump = list(
    description = "BNS test on GBM (no jump), alpha=0.975",
    n_returns = length(ret_gbm),
    z_statistic = bns_nojump$ztest,
    p_value = bns_nojump$pvalue,
    critical_value = bns_nojump$critical.value
  ),
  case6_bns_jump = list(
    description = "BNS test on jump series, alpha=0.975",
    n_returns = length(ret_jmp),
    z_statistic = bns_jump$ztest,
    p_value = bns_jump$pvalue,
    critical_value = bns_jump$critical.value
  ),
  case7_multi_asset = list(
    description = "2-asset rCov, seed=7 n=50",
    n_returns = n_multi - 1,
    rcov_00 = rcov_multi[1, 1],
    rcov_01 = rcov_multi[1, 2],
    rcov_10 = rcov_multi[2, 1],
    rcov_11 = rcov_multi[2, 2]
  ),
  case8_aggregate = list(
    description = "aggregatePrice sampleTData 30s",
    n_obs = n_agg,
    first_price = first_price_agg,
    last_price = last_price_agg
  ),
  case9_agg_returns = list(
    description = "makeReturns on aggregated 30s prices",
    n_returns = length(ret_agg),
    rv = rv_agg
  )
)

## ---- Write JSON ----
out_path <- "tests/fixtures/hfe/baselines.json"
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(baselines, auto_unbox = TRUE, digits = 17, pretty = TRUE),
           out_path)

cat("\n=== Baselines written to:", out_path, "===\n")
cat("R version:", R.version.string, "\n")
cat("highfrequency version:", as.character(packageVersion("highfrequency")), "\n\n")

## ---- Console summary (for hardcoding into C++ test) ----
cat("=== CASE 1 (constant) ===\n")
cat(sprintf("  rv=%.17g  rvol=%.17g  rq=%.17g  bpv=%.17g  rsv+=%.17g  rsv-=%.17g\n",
            rv_const, rvol_const, rq_const, bpv_const,
            rsv_const$rSVarupside, rsv_const$rSVardownside))

cat("=== CASE 2 (known returns) ===\n")
cat(sprintf("  rv=%.17g  rvol=%.17g  rq=%.17g  bpv=%.17g  rsv+=%.17g  rsv-=%.17g\n",
            rv_known, rvol_known, rq_known, bpv_known,
            rsv_known$rSVarupside, rsv_known$rSVardownside))

cat("=== CASE 3 (GBM seed=42) ===\n")
cat(sprintf("  rv=%.17g  rvol=%.17g  rq=%.17g  bpv=%.17g  rsv+=%.17g  rsv-=%.17g\n",
            rv_gbm, rvol_gbm, rq_gbm, bpv_gbm,
            rsv_gbm$rSVarupside, rsv_gbm$rSVardownside))
cat("  prices[1:5]:", paste(sprintf("%.17g", prices_gbm[1:5]), collapse=", "), "\n")

cat("=== CASE 4 (GBM+jump seed=123) ===\n")
cat(sprintf("  rv=%.17g  bpv=%.17g  rq=%.17g  rsv+=%.17g  rsv-=%.17g\n",
            rv_jmp, bpv_jmp, rq_jmp,
            rsv_jmp$rSVarupside, rsv_jmp$rSVardownside))

cat("=== CASE 5 (BNS no-jump) ===\n")
cat(sprintf("  z=%.17g  p=%.17g  crit=%.17g\n",
            bns_nojump$ztest, bns_nojump$pvalue, bns_nojump$critical.value))

cat("=== CASE 6 (BNS jump) ===\n")
cat(sprintf("  z=%.17g  p=%.17g  crit=%.17g\n",
            bns_jump$ztest, bns_jump$pvalue, bns_jump$critical.value))

cat("=== CASE 7 (multi-asset rCov) ===\n")
cat(sprintf("  [0,0]=%.17g  [0,1]=%.17g  [1,0]=%.17g  [1,1]=%.17g\n",
            rcov_multi[1,1], rcov_multi[1,2],
            rcov_multi[2,1], rcov_multi[2,2]))
cat("  p1[1:3]:", paste(sprintf("%.17g", p1[1:3]), collapse=", "), "\n")
cat("  p2[1:3]:", paste(sprintf("%.17g", p2[1:3]), collapse=", "), "\n")

cat("=== CASE 8 (aggregatePrice) ===\n")
cat(sprintf("  n_obs=%d  first=%.17g  last=%.17g\n",
            n_agg, first_price_agg, last_price_agg))

cat("=== CASE 9 (agg returns RV) ===\n")
cat(sprintf("  n_ret=%d  rv=%.17g\n", length(ret_agg), rv_agg))

cat("\nDONE.\n")
