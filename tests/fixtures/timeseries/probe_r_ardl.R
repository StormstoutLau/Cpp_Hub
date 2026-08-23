# probe_r_ardl.R - R ARDL bounds_test/t_test 语义探针 (Phase 8 M1 Step 4)
# 数据: ardl_smoke_data.csv (verify_ardl.py 生成, k1 夹具)
# 用法: & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" probe_r_ardl.R

suppressMessages(library(ARDL))
# CSV 无 header 行 (列名在尾部 #columns= 注释, verify_ardl.py 约定)
lines <- readLines("ardl_smoke_data.csv")
hdr <- sub("^#columns=", "", tail(grep("^#columns=", lines, value = TRUE), 1))
d <- read.csv("ardl_smoke_data.csv", header = FALSE, comment.char = "#")
names(d) <- strsplit(hdr, ",")[[1]]
cat("cols:", paste(names(d), collapse = ","), "\n")
cat("nrow:", nrow(d), "\n\n")

# p=2, q1=1 (与 statsmodels c_p2q1 对拍点一致)
m <- ardl(k1_y ~ k1_x1, data = d, order = c(2, 1))
cat("=== ardl order c(2,1) summary ===\n")
cat("order:", paste(m$order, collapse = ","), "\n")
cat("nobs:", length(m$coefficients) + m$df.residual, " K:", length(m$coefficients), "\n\n")

for (cs in 2:3) {
  cat("=== bounds_f_test case", cs, "===\n")
  bt <- bounds_f_test(m, case = cs)
  print(bt)
  cat("\n")
}

for (cs in c(1, 3)) {
  cat("=== bounds_t_test case", cs, "(含截距模型) ===\n")
  tryCatch({
    bt <- bounds_t_test(m, case = cs)
    print(bt)
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
  cat("\n")
}

# case 1: 无截距模型
m_n <- ardl(k1_y ~ k1_x1 - 1, data = d, order = c(2, 1))
cat("=== bounds_f_test case 1 (无截距模型) ===\n")
print(bounds_f_test(m_n, case = 1))
cat("\n=== bounds_t_test case 1 (无截距模型) ===\n")
print(bounds_t_test(m_n, case = 1))
cat("\n")

# case 4/5: 截距+趋势模型 (time() 确定性趋势, order 不含)
m_ct <- tryCatch(
  ardl(k1_y ~ k1_x1 + time(k1_y), data = d, order = c(2, 1)),
  error = function(e) {
    cat("time() failed:", conditionMessage(e), "\n")
    NULL
  })
if (!is.null(m_ct)) {
cat("=== bounds_f_test case 4 (ct 模型) ===\n")
tryCatch(print(bounds_f_test(m_ct, case = 4)), error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
cat("\n=== bounds_f_test case 5 (ct 模型) ===\n")
tryCatch(print(bounds_f_test(m_ct, case = 5)), error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
cat("=== bounds_t_test case 5 (ct 模型) ===\n")
tryCatch(print(bounds_t_test(m_ct, case = 5)), error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
cat("\n")
}

# 结构探查
cat("=== str(bounds_f_test case 3) ===\n")
bt <- bounds_f_test(m, case = 3)
str(bt, max.level = 2)
