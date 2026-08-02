userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

ret_known <- c(0.01, -0.02, 0.03, -0.01, 0.02)
n <- length(ret_known)
sum_r4 <- sum(ret_known^4)

cat("n =", n, "\n")
cat("sum r^4 =", sprintf("%.17g", sum_r4), "\n")
cat("(n/3)*sum r^4 =", sprintf("%.17g", (n/3)*sum_r4), "\n")
cat("rQuar(ret_known) =", sprintf("%.17g", rQuar(ret_known)), "\n")
cat("rQuar source:\n")
print(getAnywhere(rQuar))
