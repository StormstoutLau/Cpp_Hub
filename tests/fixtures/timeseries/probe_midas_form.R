.libPaths(c('F:/R/win-library/4.6', file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6"), .libPaths()))
suppressMessages(library(midasr))
set.seed(42)
n_lf <- 100; m <- 3
x_hf <- as.numeric(arima.sim(list(ar = 0.5), n = n_lf * m,
                             rand.gen = function(n) rnorm(n, 0, 0.8)))
w_true <- nealmon(c(5, -0.5), 4, NULL)
x_lag <- mls(x_hf, 0:3, m)
y_lf <- 2 + as.vector(x_lag %*% w_true) + rnorm(n_lf, 0, 0.3)
ctrl <- list(reltol = 1e-12, maxit = 10000)

cat("--- Form A: weight in formula, plain start ---\n")
r <- try({
  mrA <- midas_r(y_lf ~ mls(x_hf, 0:3, m, nealmon),
                 start = list(x_hf = c(5, -0.5)), control = ctrl)
  cat("A coef:", paste(sprintf("%.10g", coef(mrA)), collapse = ", "), "\n")
  cat("A conv:", mrA$convergence, " SSR:", sprintf("%.10g", deviance(mrA)), "\n")
})
if (inherits(r, "try-error")) cat("A FAILED\n")

cat("--- Form B: weight call in start ---\n")
r <- try({
  mrB <- midas_r(y_lf ~ mls(x_hf, 0:3, m),
                 start = list(x_hf = nealmon(c(5, -0.5), 4)), control = ctrl)
  cat("B coef:", paste(sprintf("%.10g", coef(mrB)), collapse = ", "), "\n")
})
if (inherits(r, "try-error")) cat("B FAILED\n")
