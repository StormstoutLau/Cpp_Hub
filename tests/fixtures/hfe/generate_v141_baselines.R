## =============================================================================
## generate_v141_baselines.R
## Phase 5 v1.4.1 - Realized Kernel R baseline 生成
##
## 对标: R highfrequency 1.0.3 rKernelCov (BNS 2008 ECTA)
## 输出: 硬编码到 tests/unit/hfecon/test_realized_kernel.cpp
## =============================================================================
userLib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(userLib, .libPaths()))
suppressMessages(library(highfrequency))

cat("=== highfrequency version:", as.character(packageDescription("highfrequency")$Version), "===\n\n")

## ============================================================================
## CASE A: 核函数反向验证 (通过 rKernelCov 反推 KK() 值)
## 构造 r = [1, 1, 0, ..., 0] (n=30), 使得:
##   gamma_0 = 2, gamma_1 = 1, gamma_h = 0 (h>=2)
##   RK(H=1, no DOF) = gamma_0 + 2*k(0)*gamma_1 = 2 + 2*1 = 4  (所有核 k(0)=1)
##   RK(H=h, no DOF) = 2 + 2*k((h-1)/h)*1 = 2 + 2*k((h-1)/h)
##   => k((h-1)/h) = (RK - 2) / 2
## ============================================================================
cat("=== CASE A: Kernel reverse engineering ===\n")
n <- 30
ret_a <- c(1, 1, rep(0, n - 2))
cat(sprintf("ret = [1, 1, 0, ..., 0], n = %d\n", n))
cat(sprintf("Expected: gamma_0 = 2, gamma_1 = 1, gamma_h = 0 (h>=2)\n\n"))

kernels <- listAvailableKernels()
cat("Available kernels:", paste(kernels, collapse=", "), "\n\n")

## 反推 k(x) at x = 0, 1/2, 2/3, 3/4, 4/5, 5/6, 6/7, 7/8, 8/9, 9/10
x_targets <- c(0.0, 1/2, 2/3, 3/4, 4/5, 5/6, 6/7, 7/8, 8/9, 9/10)
cat("Reverse-engineered k(x) values:\n")
cat(sprintf("%-25s %s\n", "Kernel", paste(sprintf("%10.6f", x_targets), collapse="  ")))
for (kn in kernels) {
  kvals <- numeric(length(x_targets))
  for (i in seq_along(x_targets)) {
    x <- x_targets[i]
    if (x == 0) {
      H <- 1
    } else {
      H <- round(x / (1 - x))  ## x = (H-1)/H => H = x/(1-x) => round
      if (H < 1) H <- 1
    }
    rk <- rKernelCov(rData = ret_a, kernelType = kn, kernelParam = H,
                     kernelDOFadj = FALSE)
    rk_val <- as.numeric(rk)
    kvals[i] <- (rk_val - 2) / 2
  }
  cat(sprintf("%-25s %s\n", kn, paste(sprintf("%10.6f", kvals), collapse="  ")))
}

## ============================================================================
## CASE B: Realized Kernel 主估计量 R baseline
## 使用已知收益率序列, 多种核 x 多种 bandwidth x DOF on/off
## ============================================================================

cat("\n\n=== CASE B: rKernelCov R baselines ===\n")

## B.1: 简单已知收益率 r = [0.01, -0.02, 0.03, -0.01, 0.02] (n=5)
ret_b1 <- c(0.01, -0.02, 0.03, -0.01, 0.02)
cat(sprintf("\n--- B.1: ret = [0.01, -0.02, 0.03, -0.01, 0.02], n=%d ---\n", length(ret_b1)))
cat(sprintf("RV (gamma_0) = %.18f\n", sum(ret_b1^2)))

for (kn in c("rectangular", "Bartlett", "Second", "Epanechnikov", "Cubic",
             "Fifth", "Sixth", "Seventh", "Eighth", "Parzen",
             "TukeyHanning", "ModifiedTukeyHanning")) {
  for (H in c(1, 2, 3)) {
    for (dof in c(TRUE, FALSE)) {
      rk <- rKernelCov(rData = ret_b1, kernelType = kn, kernelParam = H,
                       kernelDOFadj = dof)
      rk_val <- as.numeric(rk)
      cat(sprintf("rKernelCov(kernel=%-20s H=%d DOF=%-5s) = %.18f\n",
                  paste0('"', kn, '"'), H, as.character(dof), rk_val))
    }
  }
}

## B.2: GBM 收益率 (n=100, seed=42, sigma=0.01) — 使用 CASE3 价格
## 从 baselines.json 读取 CASE3 价格
cat("\n--- B.2: GBM n=100 seed=42 sigma=0.01 (CASE3 prices) ---\n")
case3_prices <- c(
  101.38039917595452, 100.80951930678135, 101.17625276642408, 101.81858984699473,
  102.23104330040793, 102.12260864846189, 103.67793931303507, 103.57984520742802,
  105.69176743555116, 105.62550457608084, 107.01281132027738, 109.48800645723944,
  107.97788160279309, 107.67727062775589, 107.5338095049555,  108.21985032207174,
  107.91266902871402, 105.08375782085513, 102.5502637483704,  103.91301863018465,
  103.59486924448341, 101.76586358015204, 101.59106069971737, 102.83258659625074,
  104.80004776280842, 104.34988550719638, 104.08177023796344, 102.26272240106542,
  102.73431454136686, 102.07891967627012, 102.54489858699878, 103.27022650783616,
  104.34473178668154, 103.71127977897947, 104.23629964143039, 102.46183081603857,
  101.66120415375512, 100.79983118852508, 98.395454060551174, 98.431003483923988,
  98.633978964810709, 98.278495919828828, 99.02643907679554,  98.309417634148389,
  96.973429388244014, 97.394057491444727, 96.60700611354865,  98.012231123074883,
  97.5902719916715,   98.232222712900949, 98.548966622169928, 97.779520992766194,
  99.332462747317521, 99.973127669179803, 100.06290448109721, 100.3400121863053,
  101.02393092902517, 101.11472441741427, 98.13311336912875,  98.413076475124512,
  98.052332359422962, 98.234123562975157, 98.807338936726623, 100.2001064478679,
  99.474002678242329, 100.7781691944837,  101.11719977608917, 102.17277970661081,
  103.11785780390004, 103.86389771145825, 102.78610482778274, 102.69344754230983,
  103.33576022608229, 102.35511239852755, 101.8010046441415,  102.39418642982434,
  103.18378568270319, 103.66342999507516, 102.74925863684794, 101.62543305377729,
  103.17441434013229, 103.44086674398767, 103.53239054954682, 103.40729910498969,
  102.17962169825519, 102.80687524182999, 102.58388274191515, 102.39657502742747,
  103.3567636733232,  104.20962124264086, 105.67048531652755, 105.16850611660557,
  105.85469687979499, 107.33754273540536, 106.15184672834231, 105.2420209750828,
  104.05767083080239, 102.55027156853281, 102.63232670461113, 103.30491982877388
)
ret_b2 <- makeReturns(case3_prices)
cat(sprintf("n = %d, RV = %.18f\n\n", length(ret_b2), sum(ret_b2^2)))

for (kn in c("rectangular", "Bartlett", "Parzen", "TukeyHanning")) {
  for (H in c(1, 3, 5, 10)) {
    rk <- rKernelCov(rData = ret_b2, kernelType = kn, kernelParam = H,
                     kernelDOFadj = TRUE)
    rk_val <- as.numeric(rk)
    cat(sprintf("rKernelCov(kernel=%-20s H=%-2d DOF=TRUE ) = %.18f\n",
                paste0('"', kn, '"'), H, rk_val))
  }
}

## B.3: 噪声稳健性测试 — 纯噪声序列, RK 应显著小于 RV
cat("\n--- B.3: Pure noise series (n=1000, noise sd=0.001) ---\n")
set.seed(999)
noise <- rnorm(1000, mean = 0, sd = 0.001)
cat(sprintf("n = %d, RV = %.18f\n", length(noise), sum(noise^2)))
for (kn in c("rectangular", "Bartlett", "Parzen")) {
  for (H in c(5, 10, 20)) {
    rk <- rKernelCov(rData = noise, kernelType = kn, kernelParam = H,
                     kernelDOFadj = TRUE)
    rk_val <- as.numeric(rk)
    cat(sprintf("rKernelCov(kernel=%-20s H=%-2d DOF=TRUE ) = %.18f  (RK/RV = %.6f)\n",
                paste0('"', kn, '"'), H, rk_val, rk_val / sum(noise^2)))
  }
}

## B.4: gamma_1 提取 (用于验证 gamma_1 字段)
cat("\n--- B.4: gamma_1 verification ---\n")
cat(sprintf("ret_b1 gamma_1 = sum(r[1:n-1]*r[2:n]) = %.18f\n",
            sum(ret_b1[1:(length(ret_b1)-1)] * ret_b2[2:length(ret_b1)])))
## rKernelCov with H=1, DOF=FALSE, Rectangular:
## RK = gamma_0 + 2*k(0)*gamma_1 = gamma_0 + 2*gamma_1
## => gamma_1 = (RK - gamma_0) / 2
rk_h1_rect_nodof <- as.numeric(rKernelCov(rData = ret_b1, kernelType = "rectangular",
                                            kernelParam = 1, kernelDOFadj = FALSE))
gamma_0_b1 <- sum(ret_b1^2)
gamma_1_b1 <- (rk_h1_rect_nodof - gamma_0_b1) / 2
cat(sprintf("ret_b1: gamma_0 = %.18f, RK(H=1,rect,nodof) = %.18f => gamma_1 = %.18f\n",
            gamma_0_b1, rk_h1_rect_nodof, gamma_1_b1))

## 手算 gamma_1 for ret_b1 = [0.01, -0.02, 0.03, -0.01, 0.02]
## gamma_1 = sum_{i=0}^{3} r[i]*r[i+1] = 0.01*(-0.02) + (-0.02)*0.03 + 0.03*(-0.01) + (-0.01)*0.02
##         = -0.0002 - 0.0006 - 0.0003 - 0.0002 = -0.0013
cat(sprintf("ret_b1: hand-computed gamma_1 = %.18f\n", -0.0013))

cat("\nDONE.\n")
