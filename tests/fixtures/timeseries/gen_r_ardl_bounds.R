# gen_r_ardl_bounds.R - R ARDL bounds F/t 交叉基准生成 (Phase 8 M1 Step 4)
#
# R ARDL 0.2.5 语义 (probe_r_ardl.R 探明, 2026-08-22):
#   - case 1: 公式 y ~ x - 1 (无截距); case 2/3: y ~ x (截距);
#     case 4/5: y ~ x + time(y) (截距 + 确定性趋势, order 不含 trend)
#   - case-trend 不相容 → throw (parse_case 实录, 本库 Step 4 校验镜像)
#   - k = 独立变量数 (PSS 语义, null.value 实录; ≠ statsmodels 偏移)
#   - t 检验仅 case 1/3/5 (case 2/4 报错, PSS Table CII 边界)
#   - F/t 对 trend 列线性变换不变 (Wald 不变性) → time() 尺度无关
#
# 用法: & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" gen_r_ardl_bounds.R
# 产出: r_ardl_bounds.csv (name,F,t; t 无该 case 则 NA)
suppressMessages(library(ARDL))

lines <- readLines("ardl_smoke_data.csv")
hdr <- sub("^#columns=", "", tail(grep("^#columns=", lines, value = TRUE), 1))
d <- read.csv("ardl_smoke_data.csv", header = FALSE, comment.char = "#")
names(d) <- strsplit(hdr, ",")[[1]]

out <- function(name, model, case, has_t) {
  f <- bounds_f_test(model, case = case)$statistic[[1]]
  t <- if (has_t) bounds_t_test(model, case = case)$statistic[[1]] else NA
  lines <<- c(lines, sprintf("%s,%.17g,%.17g", name, f, t))
}
lines <- character(0)

# k1: p=2, q=[1] — 五 case
m_k1_c <- ardl(k1_y ~ k1_x1, data = d, order = c(2, 1))
m_k1_n <- ardl(k1_y ~ k1_x1 - 1, data = d, order = c(2, 1))
m_k1_ct <- ardl(k1_y ~ k1_x1 + time(k1_y), data = d, order = c(2, 1))
out("uecm_c1_k1", m_k1_n, 1, TRUE)
out("uecm_c2_k1", m_k1_c, 2, FALSE)
out("uecm_c3_k1", m_k1_c, 3, TRUE)
out("uecm_c4_k1", m_k1_ct, 4, FALSE)
out("uecm_c5_k1", m_k1_ct, 5, TRUE)

# k2: p=2, q=[1,2] — case 3/5 (k_PSS=2)
m_k2_c <- ardl(k2_y ~ k2_x1 + k2_x2, data = d, order = c(2, 1, 2))
m_k2_ct <- ardl(k2_y ~ k2_x1 + k2_x2 + time(k2_y), data = d, order = c(2, 1, 2))
out("uecm_c3_k2", m_k2_c, 3, TRUE)
out("uecm_c5_k2", m_k2_ct, 5, TRUE)

writeLines(lines, "r_ardl_bounds.csv")
cat("written: r_ardl_bounds.csv\n")
