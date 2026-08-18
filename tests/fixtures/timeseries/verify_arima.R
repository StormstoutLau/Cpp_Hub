# verify_arima.R - ARIMA 基准生成: R stats::arima CSS/CSS-ML (Phase 7C M1)
#
# 对照点 (spec §1.3, AR6 配对):
#   ARIMA CSS vs R stats::arima method="CSS"    1e-10 (参数级实测校准)
#   ARIMA CSS-ML vs R method="CSS-ML"           1e-8
#   d=1 drift vs forecast::Arima (AR5)          1e-8
#
# 语义对齐 (一手实测留痕):
#   - d=0 夹具: include.mean=FALSE (C++ CSS 无均值项, spec 决策 3;
#     数据本身零均值生成 + demean 语义留给 innovations 对照)
#   - AR2 复核点: R CSS 的 n.cond 实测 = max(p, q)+? — 本脚本对 CSS 残差
#     长度/数值打印, 若与 C++ (n_cond−d = p 起) 不一致 → 以 R 实测为准回写
#   - MA 正号 (AR1): R coef 里 ma1 即 (1+θB) 的 θ
#   - loglik: R 纯 CSS 报 "part log-likelihood" (AR3 不与 ML 混比, 只对
#     CSS-ML/ML 的 loglik 数值对照)
#
# 用法: Rscript verify_arima.R  (读 arima_smoke_data.csv, 输出 12 位基准)
.libPaths(c('F:/R/win-library/4.6', file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6"), .libPaths()))

here <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", here[grepl("^--file=", here)])))
csv <- file.path(here, "arima_smoke_data.csv")
# keys 从原始行读 (read.csv 会把 #columns 行按逗号拆列, 2026-08-18 实测)
lns <- readLines(csv)
keys <- strsplit(sub("^#columns=", "", lns[grepl("^#columns", lns)]), ",")[[1]]
raw <- read.csv(csv, header = FALSE, stringsAsFactors = FALSE)
dat <- raw[!grepl("^#", raw[[1]]), ]
cat("columns:", paste(keys, collapse = ", "), "\n")
arma11 <- as.numeric(dat[[1]])
arma22 <- as.numeric(dat[[2]])
arma21 <- as.numeric(dat[[3]])
arma12 <- as.numeric(dat[[4]])
arima111d_level <- as.numeric(dat[[5]])
stopifnot(!any(is.na(arima111d_level)))

fmt <- function(x) format(x, digits = 15, scientific = TRUE)

cat("=== arma11 (T=250) stats::arima CSS / CSS-ML, include.mean=FALSE ===\n")
a <- list()
for (m in c("CSS", "CSS-ML")) {
  fit <- arima(arma11, order = c(1, 0, 1), method = m, include.mean = FALSE)
  a[[m]] <- fit
  cat(m, ":")
  cat(" coef=", paste(fmt(fit$coef), collapse = ", "))
  cat(" sigma2=", fmt(fit$sigma2))
  cat(" loglik=", fmt(as.numeric(logLik(fit))))
  cat(" aic=", fmt(fit$aic))
  cat(" n.cond=", fit$n.cond, "\n")
}
cat("resid head (CSS):", paste(fmt(head(a$CSS$residuals, 3)), collapse = ", "), "\n")
cat("resid tail (CSS):", paste(fmt(tail(a$CSS$residuals, 2)), collapse = ", "), "\n")
cat("nresid length CSS:", length(a$CSS$residuals), "\n")

cat("\n=== arma22 (T=300) CSS / CSS-ML ===\n")
for (m in c("CSS", "CSS-ML")) {
  fit <- arima(arma22, order = c(2, 0, 2), method = m, include.mean = FALSE)
  cat(m, ": coef=", paste(fmt(fit$coef), collapse = ", "),
      " sigma2=", fmt(fit$sigma2),
      " loglik=", fmt(as.numeric(logLik(fit))),
      " n.cond=", fit$n.cond, "\n")
}

cat("\n=== arma21 (T=300, p=2 q=1) CSS — AR2 裁决: n.cond 实测 ===\n")
for (m in c("CSS", "CSS-ML")) {
  fit <- arima(arma21, order = c(2, 0, 1), method = m, include.mean = FALSE)
  cat(m, ": coef=", paste(fmt(fit$coef), collapse = ", "),
      " sigma2=", fmt(fit$sigma2),
      " loglik=", fmt(as.numeric(logLik(fit))),
      " n.cond=", fit$n.cond, "\n")
}

cat("\n=== arma12 (T=300, p=1 q=2) CSS — AR2 定案: max(p,q)=2 vs p=1 ===\n")
for (m in c("CSS", "CSS-ML")) {
  fit <- arima(arma12, order = c(1, 0, 2), method = m, include.mean = FALSE)
  cat(m, ": coef=", paste(fmt(fit$coef), collapse = ", "),
      " sigma2=", fmt(fit$sigma2),
      " loglik=", fmt(as.numeric(logLik(fit))),
      " n.cond=", fit$n.cond, "\n")
}

cat("\n=== arima111d: level 数据, d=1 ===\n")
# 语义裁决 (2026-08-18 实测): stats::arima 对 d≥1 默认 include.mean=FALSE
#   → 漂移未建模 → AR 多项式近单位根吸收漂移 (伪 φ=0.998, 退化路径;
#   对照 C++ arima_fit d=1 include.mean=false 行为)
# include.mean=TRUE → 差分序列均值吸收漂移 (与 statsmodels demean 对齐)
for (m in c("CSS", "CSS-ML")) {
  fit <- arima(arima111d_level, order = c(1, 1, 1), method = m)
  cat(m, " mean=FALSE: coef=", paste(fmt(fit$coef), collapse = ", "),
      " sigma2=", fmt(fit$sigma2),
      " loglik=", fmt(as.numeric(logLik(fit))),
      " n.cond=", fit$n.cond, "\n")
}
for (m in c("CSS", "CSS-ML")) {
  fit <- arima(arima111d_level, order = c(1, 1, 1), method = m,
               include.mean = TRUE)
  cat(m, " mean=TRUE: coef=", paste(fmt(fit$coef), collapse = ", "),
      " sigma2=", fmt(fit$sigma2),
      " loglik=", fmt(as.numeric(logLik(fit))),
      " n.cond=", fit$n.cond, "\n")
}
# forecast::Arima with drift (AR5): CSS-ML 语义
library(forecast)
fa <- Arima(arima111d_level, order = c(1, 1, 1), include.drift = TRUE,
            method = "CSS-ML")
cat("forecast::Arima drift: coef=", paste(fmt(fa$coef), collapse = ", "),
    " sigma2=", fmt(fa$sigma2),
    " loglik=", fmt(as.numeric(logLik(fa))),
    " drift=", fmt(fa$coef["drift"]), "\n")
fa_css <- Arima(arima111d_level, order = c(1, 1, 1), include.drift = TRUE,
                method = "CSS")
cat("forecast::Arima drift CSS: coef=", paste(fmt(fa_css$coef), collapse = ", "),
    " sigma2=", fmt(fa_css$sigma2), "\n")
