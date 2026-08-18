# -*- coding: utf-8 -*-
# dump_var_r_values.R - R 基准数值机器输出 (Phase 7C M2 .inc 生成用)
#
# 输出 F:/Cpp_Hub/tests/fixtures/timeseries/var_r_values.txt
# 格式: key<TAB>value1 value2 ... (%.17g, tab/空格分隔, 嵌套用行)
# 内容: vars 交叉 + Spillover 主基准 (g.fevd/G.spillover/net/VARselect/滚动)
.libPaths(c("F:/R/win-library/4.6", .libPaths()))
suppressMessages(library(vars))
suppressMessages(library(Spillover))
setwd("F:/Cpp_Hub/tests/fixtures/timeseries")

d <- read.csv("var_smoke_data.csv")
y <- as.matrix(d[, c("y1", "y2", "y3")])
out <- file("var_r_values.txt", "w")
emit <- function(key, vals) {
  cat(key, paste(sprintf("%.17g", vals), collapse = " "), "\n",
      file = out, sep = "\t")
}

v2 <- VAR(y, p = 2, type = "const")
coefs <- sapply(v2$varresult, function(m) coef(m))  # 回归元×方程
emit("vars_coef", as.numeric(coefs))
emit("vars_loglik", logLik(v2))

# Cholesky 下三角 (df 修正 Σ)
n <- v2$obs; K <- 3
sig_df <- crossprod(resid(v2)) / (n - 2 * K - 1)
emit("vars_chol_lower", as.numeric(t(chol(sig_df))))

# fevd orth H=10
fv <- fevd(v2, n.ahead = 10)
emit("vars_fevd_orth_h10", as.numeric(t(sapply(fv, function(m) m[10, ]))))

# VARselect (maxlags=4) — 轨迹 p=1..4 (无 p=0)
# criteria 行序 (R 实测定档): AIC(n), HQ(n), SC(n), FPE(n)
vs <- VARselect(y, lag.max = 4, type = "const")
emit("vars_varselect_aic", as.numeric(vs$criteria["AIC(n)", ]))
emit("vars_varselect_hq", as.numeric(vs$criteria["HQ(n)", ]))
emit("vars_varselect_sc", as.numeric(vs$criteria["SC(n)", ]))
emit("vars_varselect_fpe", as.numeric(vs$criteria["FPE(n)", ]))

# --- Spillover 主基准 ---
g10 <- g.fevd(v2, n.ahead = 10)
emit("gfevd_dy_h10", c(as.numeric(g10$y1[10, ]), as.numeric(g10$y2[10, ]),
                       as.numeric(g10$y3[10, ])))
g10n <- g.fevd(v2, n.ahead = 10, normalized = FALSE)
emit("gfevd_dy_h10_raw", c(as.numeric(g10n$y1[10, ]), as.numeric(g10n$y2[10, ]),
                           as.numeric(g10n$y3[10, ])))
sp <- G.spillover(v2, n.ahead = 10, standardized = TRUE)
emit("gspillover_table_h10", as.numeric(sp[1:3, 1:3]))  # θ̃·100/K 表
emit("gspillover_to", as.numeric(sp[4, 1:3]))
emit("gspillover_from", as.numeric(sp[1:3, 4]))
emit("gspillover_tci", as.numeric(sp[4, 4]))
nt <- net(sp)
emit("net", as.numeric(nt[, 3]))  # Net 列

# H=50 敏感性 (V10)
g50 <- g.fevd(v2, n.ahead = 50)
emit("gfevd_dy_h50", c(as.numeric(g50$y1[50, ]), as.numeric(g50$y2[50, ]),
                       as.numeric(g50$y3[50, ])))
sp50 <- G.spillover(v2, n.ahead = 50, standardized = TRUE)
emit("gspillover_tci_h50", as.numeric(sp50[4, 4]))

# 滚动 (w=150, step=1, 泛化) — 逐窗口 TCI
yz <- as.zoo(y)
dyr <- roll.spillover(data = yz, width = 150, n.ahead = 10,
                      index = "generalized")
emit("roll_tci_all", as.numeric(dyr))

close(out)
cat("written var_r_values.txt\n")
