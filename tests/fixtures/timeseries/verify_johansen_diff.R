# verify_johansen_diff.R - urca 侧 Johansen 双库 diff 网格 (M3 前置, §6.2.1)
#
# 冻结结论见 docs/phases/phase7/JOHANSEN_DUAL_LIB_DIFF.md:
#   - 参数映射: urca ca.jo(K) ≙ statsmodels coint_johansen(k_ar_diff = K-1)
#   - 情形映射: ecdet="none" ≡ SM det_order=0; ecdet="const"/"trend" 无 SM 对应
#   - spec transitory/longrun 数学恒等 (差被 Z1 张成), 仅转储 transitory
#
# 网格: K ∈ {2, 3} × ecdet ∈ {none, const, trend} (transitory)
# 输出: johansen_urca_grid.txt (teststat/cval 统一翻转为 r=0..N-1 与 SM 一致)

.libPaths(c("F:/R/win-library/4.6", .libPaths()))
library(urca)

HERE <- "F:/Cpp_Hub/tests/fixtures/timeseries"
y <- as.matrix(read.csv(file.path(HERE, "coint_smoke_data.csv")))

emit <- function(fh, tag, obj) {
  a <- as.numeric(as.matrix(obj))
  write(sprintf("%s\t%s", tag, paste(sprintf("%.17g", a), collapse = " ")),
        file = fh, append = TRUE)
}

out <- file.path(HERE, "johansen_urca_grid.txt")
fh <- file(out, open = "wt")
for (K in c(2, 3)) {
  for (ec in c("none", "const", "trend")) {
    tag <- sprintf("ur_K%d_%s", K, ec)
    jtr <- ca.jo(y, type = "trace", ecdet = ec, K = K, spec = "transitory")
    jmx <- ca.jo(y, type = "eigen", ecdet = ec, K = K, spec = "transitory")
    # urca 原始序 [r<=N-1, ..., r=0] → 翻转为 [r=0, ..., r<=N-1] (SM 行序)
    emit(fh, paste0(tag, "_eig"), jtr@lambda)
    emit(fh, paste0(tag, "_lr1"), rev(jtr@teststat))
    emit(fh, paste0(tag, "_lr2"), rev(jmx@teststat))
    # cval 行序同 teststat (行1=r<=N-1) → 翻转行序: 行1=r=0
    emit(fh, paste0(tag, "_cvt"), jtr@cval[nrow(jtr@cval):1, ])
    emit(fh, paste0(tag, "_cvm"), jmx@cval[nrow(jmx@cval):1, ])
    emit(fh, paste0(tag, "_nobs"), nrow(jtr@Z0))
  }
}
close(fh)
cat("wrote", out, "\n")
