# verify_za.R - urca ur.za 对照基准生成 (Phase 7C M0, Model A/B/C 主对照)
#
# 用途: 同一合成数据 (za_smoke_data.csv, 由 verify_za.py 生成) 上
#       urca ur.za 固定 lag=1 基准, 供 test_zivot_andrews.cpp 断言 (spec §3.2.2, 1e-8)
# 映射约定 (2026-08-17 实测): urca bpoint = C++ break_index + 1
# 环境: R 4.6.1; urca 1.3-4 位于 F:/R/win-library/4.6 (沙箱限制, 双库并列加载)
.libPaths(c('F:/R/win-library/4.6',
            file.path(Sys.getenv('USERPROFILE'), 'R', 'win-library', '4.6'),
            .libPaths()))
library(urca)

y <- as.numeric(read.csv('za_smoke_data.csv', header = TRUE)$y)
# 注: 在本脚本所在目录运行 (Rscript verify_za.R), CSV 由 verify_za.py 生成

for (m in c('intercept', 'trend', 'both')) {
  z <- ur.za(y, model = m, lag = 1)
  cat(sprintf('%-10s stat=%.12f bpoint=%d cval=[%s]\n', m, z@teststat, z@bpoint,
              paste(sprintf('%.2f', z@cval), collapse = ' ')))
}
