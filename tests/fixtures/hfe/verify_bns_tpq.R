## verify_bns_tpq.R
## 排查 BNSjumpTest / rTPQuar 实际公式 (2026-08-02)
## 现象: C++ CASE5 z=0.315 vs R z=0.693; CASE6 z=7.647 vs R z=4.667
## RV/BPV 已对标 1e-10, 差异只能来自 vartheta (TPQ)
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

## 用 CASE2 已知收益率验证 (n=5, 手算可验证)
ret_known <- c(0.01, -0.02, 0.03, -0.01, 0.02)
n <- length(ret_known)
rv_known <- rRVar(ret_known)
bpv_known <- rBPCov(ret_known, makeReturns = FALSE)
rq_known <- rQuar(ret_known)

cat("=== CASE2 known returns (n=5) ===\n")
cat(sprintf("  RV=%.17g  BPV=%.17g  RQ=%.17g\n", rv_known, bpv_known, rq_known))

## rTPQuar 源码
cat("\n=== rTPQuar source ===\n")
print(getAnywhere("rTPQuar"))

## 计算 rTPQuar
tpq_known <- rTPQuar(ret_known)
cat(sprintf("\n  rTPQuar(ret_known)=%.17g\n", tpq_known))

## 手算 TPQ (按 C++ 实现):
## rMAX = 3 * sqrt(BPV/n)
## TPQ = (n/3) * sum_{|r_i|<=rMAX} r_i^4
rmax_cpp <- 3 * sqrt(bpv_known / n)
tpq_cpp <- (n / 3) * sum(ret_known[abs(ret_known) <= rmax_cpp]^4)
cat(sprintf("  C++ TPQ (rMAX=3*sqrt(BPV/n), n/3*sum|r^4|)=%.17g\n", tpq_cpp))
cat(sprintf("  rMAX_cpp=%.17g\n", rmax_cpp))

## BNSjumpTest 源码
cat("\n=== BNSjumpTest source ===\n")
bns_src <- getAnywhere("BNSjumpTest")
print(bns_src)

## 用 CASE2 跑 BNSjumpTest 看实际 z
cat("\n=== BNSjumpTest on CASE2 (n=5, no jump) ===\n")
bns_case2 <- BNSjumpTest(ret_known, IVestimator = "BV",
                         IQestimator = "TP", makeReturns = FALSE,
                         alpha = 0.975)
cat(sprintf("  z=%.17g  p=%.17g  crit=[%.17g, %.17g]\n",
            bns_case2$ztest, bns_case2$pvalue,
            bns_case2$critical.value[1], bns_case2$critical.value[2]))

## 手算 BNS Z (按 C++ 公式):
## Z = sqrt(n) * (RV - BPV) / sqrt(vartheta)
## vartheta (TPQ) = (pi^2/4) * TPQ
pi_sq_4 <- (pi^2) / 4
vartheta_cpp <- pi_sq_4 * tpq_cpp
z_cpp <- sqrt(n) * (rv_known - bpv_known) / sqrt(vartheta_cpp)
cat(sprintf("\n  C++ Z (sqrt(n)*(RV-BPV)/sqrt(pi^2/4 * TPQ))=%.17g\n", z_cpp))
cat(sprintf("  C++ vartheta=%.17g\n", vartheta_cpp))

## 用 R TPQ 重算
vartheta_r <- pi_sq_4 * tpq_known
z_r_tpq <- sqrt(n) * (rv_known - bpv_known) / sqrt(vartheta_r)
cat(sprintf("  Z with R rTPQuar=%.17g  (vartheta=%.17g)\n", z_r_tpq, vartheta_r))

## 检查 R BNSjumpTest 内部 vartheta 计算
## 可能 R 用不同的 IQ estimator 或系数
cat("\n=== Diff R rTPQuar vs C++ TPQ ===\n")
cat(sprintf("  R rTPQuar=%.17g  C++ TPQ=%.17g  diff=%.17g\n",
            tpq_known, tpq_cpp, tpq_known - tpq_cpp))
