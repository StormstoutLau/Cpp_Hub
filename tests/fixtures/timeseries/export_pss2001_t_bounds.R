# export_pss2001_t_bounds.R - PSS 2001 Table CII (t 界) 导出 (Phase 8 M1 Step 4)
#
# R ARDL 0.2.5 bounds_t_test 的 CV 来源实锤 (源码实录):
#   crit_val_bounds_pss2001$t$<case> — PSS 2001 Table CII 转录
#   (alpha 档 {0.1, 0.05, 0.025, 0.01}; ≠ statsmodels F 表档 {10,5,1,0.1}% —
#    PSS9 对偶面: F 表第 3 档 1% 无 2.5%, t 表第 4 档 1% 无 0.1%)
#
# 用法: & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" export_pss2001_t_bounds.R
# 产出: pss2001_t_bounds.csv (alpha, k, I0, I1, case)
suppressMessages(library(ARDL))
ns <- asNamespace("ARDL")

out <- NULL
for (cs in c("i", "iii", "v")) {
  e <- eval(parse(text = paste0("crit_val_bounds_pss2001$t$", cs)), envir = ns)
  cat(cs, ": class =", class(e), " nrow =", nrow(e),
      " names =", paste(names(e), collapse = ","), "\n")
  e$case <- cs
  out <- rbind(out, e)
}
write.csv(out, "pss2001_t_bounds.csv", row.names = FALSE)
cat("written: pss2001_t_bounds.csv  rows =", nrow(out), "\n")
cat("--- 抽查 (k=1, 全 alpha×case) ---\n")
print(subset(out, k == 1))
cat("--- F 表交叉印证 (case iii k=1 5%) ---\n")
f <- eval(parse(text = "crit_val_bounds_pss2001$f$iii"), envir = ns)
print(subset(f, k == 1 & alpha == 0.05))
