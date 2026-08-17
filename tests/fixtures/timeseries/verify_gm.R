# verify_gm.R - rugarch archm 对照基准 (Phase 7C M0, spec §2.3 Step5 次锚)
#
# 用途: 同一 gm_smoke_data.csv (verify_gm.py 生成) 上 rugarch::ugarchspec
#       archm=TRUE, archpow=1|2 拟合, 输出全精度参数/llf — C++ 对照容差 1e-4
#       (solver 敏感: rugarch hybrid/nloptr vs arch SLSQP 落点差, spec 风险 #6)
# 注: rugarch 无 log 形 (GM1: archpow 仅 1=σ/2=σ²) — log 形唯一锚 = arch
#
# 三步法 (spec §2.3 Step 5.3, fix() 互验 — 隔离似然差 vs 优化器落点差):
#   5.1 C++ vs arch (主锚, verify_gm.py, 1e-4 实测)
#   5.2 C++ vs rugarch (本脚本, 1e-4)
#   5.3 rugarch 参数 → arch fix 似然 (Python 侧 verify_gm_fix.py, 可选)
#
# 运行: Rscript verify_gm.R  (工作目录 = 本脚本所在目录)
.libPaths(c('F:/R/win-library/4.6',
            file.path(Sys.getenv('USERPROFILE'), 'R', 'win-library', '4.6'),
            .libPaths()))
library(rugarch)

y <- as.numeric(read.csv('gm_smoke_data.csv', header = TRUE)$r)

for (pw in c(1, 2)) {
  spec <- ugarchspec(
    variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), archm = TRUE, archpow = pw),
    distribution.model = 'norm')
  fit <- ugarchfit(spec, data = y, solver = 'hybrid')
  cf <- coef(fit)
  # archpow=1: mu, archm, omega, alpha1, beta1; archpow=2 同结构
  cat(sprintf('archpow=%d  mu=%.12f  archm=%.12f  omega=%.12f  alpha1=%.12f  beta1=%.12f\n',
              pw, cf['mu'], cf['archm'], cf['omega'], cf['alpha1'], cf['beta1']))
  cat(sprintf('           llf=%.10f  (与 arch/Cpp 落点差容差 1e-4, spec 风险 #6)\n',
              likelihood(fit)))
}
