# verify_po.R - Phillips-Ouliaris 基准生成: urca ca.po (Phase 7C M3)
#
# 对照点 (spec §5.3):
#   Pu/Pz 统计量 vs urca ca.po 1e-8 (CI12: 双实现, Pz 优先)
#   Pu 方向依赖 / Pz 方向无关 (对调列后 Pz 统计量不变)
#
# 语义 (一手源码 ca_po_source.txt, urca 转储):
#   zl = z[2:n,], zr = z[1:(n-1,)] (水平对滞后水平, 非差分!)
#   res = resid(lm(zl ~ zr [+det]))           → Bartlett LRV Ω:
#     lmax = trunc(4·(nobs/100)^0.25) [short] / trunc(12·(nobs/100)^0.25) [long]
#     Ω = res'res/nobs + (1/nobs)·Σ_l (1-l/(lmax+1))·(Γ_l + Γ_l')
#   Pu: resu = resid(z[,1] ~ z[,-1] [+det]); w112 = Ω11 - Ω21'Ω22⁻¹Ω21
#       stat = nobs·w112/mean(resu²)         (方向依赖, 列序敏感)
#   Pz: Mzz = zl'zl/nobs; stat = nobs·tr(Ω·Mzz⁻¹)  (方向无关)
#   CV 表内嵌 (m-1 × 3 pct × 3 demean), m = 列数
#
# 输出: po_urca_baselines.txt (stat + cval + lmax)

.libPaths(c("F:/R/win-library/4.6", .libPaths()))
library(urca)

HERE <- "F:/Cpp_Hub/tests/fixtures/timeseries"
y <- as.matrix(read.csv(file.path(HERE, "coint_smoke_data.csv")))

emit <- function(fh, tag, obj) {
  a <- as.numeric(as.matrix(obj))
  write(sprintf("%s\t%s", tag, paste(sprintf("%.17g", a), collapse = " ")),
        file = fh, append = TRUE)
}

pairs <- list(
  list("y1y2", y[, c(1, 2)]),   # 协整对
  list("y2y1", y[, c(2, 1)]),   # 列对调 (Pu 变, Pz 不变 — CI12 断言点)
  list("y1y3", y[, c(1, 3)])    # 无协整
)

out <- file.path(HERE, "po_urca_baselines.txt")
fh <- file(out, open = "wt")
for (pr in pairs) {
  for (tp in c("Pu", "Pz")) {
    for (dm in c("none", "constant", "trend")) {
      for (lg in c("short", "long")) {
        r <- ca.po(pr[[2]], demean = dm, lag = lg, type = tp)
        base <- sprintf("po_%s_%s_%s_%s", pr[[1]], tp, dm, lg)
        emit(fh, base, r@teststat)
        emit(fh, paste0(base, "_cval"), r@cval)   # [10%, 5%, 1%]
        emit(fh, paste0(base, "_lmax"), r@lag)
      }
    }
  }
}
close(fh)
cat("wrote", out, "\n")
