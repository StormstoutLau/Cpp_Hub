## verify_bns_internals.R
## 排查 BNSjumpTest 内部函数: tt / hatIV / hatIQ
## 目标: 确定 theta = tt("BV") 的值, 以及 hatIQ("TP") 的实际公式
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

ret_known <- c(0.01, -0.02, 0.03, -0.01, 0.02)
n <- length(ret_known)

## ---- tt 函数源码 ----
cat("=== tt function source ===\n")
print(getAnywhere("tt"))

## ---- hatIV 函数源码 ----
cat("\n=== hatIV function source ===\n")
print(getAnywhere("hatIV"))

## ---- hatIQ 函数源码 ----
cat("\n=== hatIQ function source ===\n")
print(getAnywhere("hatIQ"))

## ---- 直接调用 ----
cat("\n=== Direct calls on CASE2 (n=5) ===\n")
cat(sprintf("  tt('BV')=%.17g\n", highfrequency:::tt("BV")))
cat(sprintf("  tt('RV')=%.17g\n", highfrequency:::tt("RV")))
cat(sprintf("  hatIV(ret_known, 'BV')=%.17g\n", highfrequency:::hatIV(ret_known, "BV")))
cat(sprintf("  hatIV(ret_known, 'RV')=%.17g\n", highfrequency:::hatIV(ret_known, "RV")))
cat(sprintf("  hatIQ(ret_known, 'TP')=%.17g\n", highfrequency:::hatIQ(ret_known, "TP")))
cat(sprintf("  hatIQ(ret_known, 'RQ')=%.17g\n", highfrequency:::hatIQ(ret_known, "RQ")))
cat(sprintf("  rTPQuar(ret_known)=%.17g\n", rTPQuar(ret_known)))
cat(sprintf("  rQuar(ret_known)=%.17g\n", rQuar(ret_known)))

## ---- 反推 theta ----
## Z = sqrt(N) * (RV - BPV) / sqrt((theta-2) * product)
## CASE2: R BNS z=-0.18669369174861714
rv <- rRVar(ret_known)
bpv <- rBPCov(ret_known, makeReturns=FALSE)
tpq <- rTPQuar(ret_known)
rq <- rQuar(ret_known)
z_r <- -0.18669369174861714
N <- n

## 反推: (theta-2) * tpq = (sqrt(N)*(rv-bpv)/z_r)^2
vartheta_implied <- (sqrt(N) * (rv - bpv) / z_r)^2
theta_minus_2 <- vartheta_implied / tpq
cat(sprintf("\n=== Reverse-engineered theta (BV, TP) ===\n"))
cat(sprintf("  vartheta_implied=%.17g\n", vartheta_implied))
cat(sprintf("  theta-2=%.17g  theta=%.17g\n", theta_minus_2, theta_minus_2 + 2))
cat(sprintf("  pi/2=%.17g  pi^2/4=%.17g\n", pi/2, pi^2/4))

## ---- 验证: 用反推的 theta 重算 Z ----
z_verify <- sqrt(N) * (rv - bpv) / sqrt(theta_minus_2 * tpq)
cat(sprintf("  Z verify=%.17g  (should match R z=%.17g)\n", z_verify, z_r))

## ---- 用 CASE3 (GBM, n=100) 验证 ----
set.seed(42)
prices_gbm <- 100 * exp(cumsum(rnorm(100, 0, 0.01)))
ret_gbm <- makeReturns(prices_gbm)
rv3 <- rRVar(ret_gbm)
bpv3 <- rBPCov(ret_gbm, makeReturns=FALSE)
tpq3 <- rTPQuar(ret_gbm)
bns3 <- BNSjumpTest(ret_gbm, IVestimator="BV", IQestimator="TP",
                    makeReturns=FALSE, alpha=0.975)
cat(sprintf("\n=== CASE3 (GBM n=100) verification ===\n"))
cat(sprintf("  R BNS z=%.17g  p=%.17g\n", bns3$ztest, bns3$pvalue))
cat(sprintf("  RV=%.17g  BPV=%.17g  TPQ=%.17g\n", rv3, bpv3, tpq3))
z3_verify <- sqrt(100) * (rv3 - bpv3) / sqrt(theta_minus_2 * tpq3)
cat(sprintf("  Z verify with theta-2=%.17g: z=%.17g\n", theta_minus_2, z3_verify))
