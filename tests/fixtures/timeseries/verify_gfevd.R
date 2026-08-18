# -*- coding: utf-8 -*-
# verify_gfevd.R - GFEVD/DY 主基准: R Spillover 0.1.1 (Phase 7C M2)
#
# 对照点 (spec §5.2.3/§5.3):
#   GFEVD (DY 框架) vs Spillover::g.fevd 1e-8 (主基准)
#   TCI/TO/FROM/NET vs Spillover::G.spillover + net 1e-8 (主基准)
#   滚动 vs Spillover::roll.spillover/dynamic.spillover (方向对照)
#
# g.fevd 一手源码语义 (Spillover 0.1.1 R/g.fevd.R L48-104 实录):
#   Sigma = summary(x)$covres       # ⚠️ df 修正版 Σ (÷(T−Kp−k)), 非 MLE
#   D = diag(1/sqrt(diag(Sigma)))
#   num[i,j](H) = Σ_{l=0..H-1} [(Φ_l Σ D)[i,j]]^2     # 逐期平方后累加
#               = σ_jj^{-1}·Σ_l [(Φ_l Σ)[i,j]]^2      # DY 框架系数 (V8)
#   den[i](H)   = Σ_{l=0..H-1} (Φ_l Σ Φ_l')[i,i]
#   θ = num/diag(den);  normalized → 行归一 (V7)
#   ⚠️ 与 PS 1998 框架差异: 分子系数 σ_jj^{-1} (DY) vs σ_ii^{-1} (PS) —
#      g.fevd 是 DY 实现; PS 轨由 C++ 自实现 + 差异断言 (V8)
#
# 包装载留档: Spillover 0.1.1 已于 2026-01-15 从 CRAN 存档 (依赖 fastSOM
#   亦存档且需编译); 本地裁剪装载 — g.fevd/G.spillover/net/dynamic.spillover
#   四函数零 fastSOM 依赖 (源码逐行核查), 仅 O.spillover/roll.net 的
#   ortho 分支被 stop() 替换; M2 主基准语义不受影响
#
# 用法: Rscript verify_gfevd.R  (读 var_smoke_data.csv, 12 位全精度)
.libPaths(c("F:/R/win-library/4.6", .libPaths()))
suppressMessages(library(vars))
suppressMessages(library(Spillover))
setwd("F:/Cpp_Hub/tests/fixtures/timeseries")

d <- read.csv("var_smoke_data.csv")
y <- as.matrix(d[, c("y1", "y2", "y3")])
cat(sprintf("T=%d K=%d\n", nrow(y), ncol(y)))

v2 <- VAR(y, p = 2, type = "const")

# --- GFEVD DY 框架 H=10 (normalized) ---
g10 <- g.fevd(v2, n.ahead = 10)           # list per-eq, each H×K
cat("== g.fevd H=10 normalized (list y1: rows=H) last row ==\n")
cat("y1 row10:\n"); print(g10$y1[10, ], digits = 12)
cat("y2 row10:\n"); print(g10$y2[10, ], digits = 12)
cat("y3 row10:\n"); print(g10$y3[10, ], digits = 12)

# 未归一 (PS/DY 分子系数核查用)
g10n <- g.fevd(v2, n.ahead = 10, normalized = FALSE)
cat("== g.fevd H=10 UNnormalized last rows ==\n")
cat("y1:\n"); print(g10n$y1[10, ], digits = 12)
cat("y2:\n"); print(g10n$y2[10, ], digits = 12)
cat("y3:\n"); print(g10n$y3[10, ], digits = 12)

# --- G.spillover 静态表 (standardized=TRUE, H=10) ---
sp <- G.spillover(v2, n.ahead = 10, standardized = TRUE)
cat("== G.spillover table ==\n")
print(sp, digits = 12)

# net (from spillover table)
nt <- net(sp)
cat("== net ==\n"); print(nt, digits = 12)

# --- H=50 敏感性 (V10) ---
g50 <- g.fevd(v2, n.ahead = 50)
cat("== g.fevd H=50 last rows ==\n")
cat("y1:\n"); print(g50$y1[50, ], digits = 12)
sp50 <- G.spillover(v2, n.ahead = 50, standardized = TRUE)
cat("H50 TCI (spillover table [K+2, K+1]):\n")
print(sp50[(nrow(sp50)-1):nrow(sp50), ncol(sp50)], digits = 12)

# --- 滚动 (window=150, step=1) — 方向对照 ---
yz <- as.zoo(y)
dy_roll <- roll.spillover(data = yz, width = 150, n.ahead = 10,
                          index = "generalized")
cat("== roll.spillover (generalized, w=150) ==\n")
cat("class:", class(dy_roll), " length:", length(dy_roll), "\n")
cat("first 3 windows:\n"); print(head(dy_roll, 3), digits = 12)
cat("last 3 windows:\n"); print(tail(dy_roll, 3), digits = 12)
