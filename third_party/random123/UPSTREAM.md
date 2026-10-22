# Random123（vendored）

- 上游: https://github.com/DEShawResearch/random123（RESEARCH §7.1 E1 锚点；raw philox.h 即取自该仓）
- 固定 commit: `9545ff6413f258be2f04c1d319d99aaef7521150`（main 分支，2026-10-22 取，经 GitHub API 页面级确认）
- 归档来源: `https://codeload.github.com/DEShawResearch/random123/zip/9545ff6413f258be2f04c1d319d99aaef7521150`（430 644 字节）
- 拷贝范围: `include/Random123/**` 原样拷贝（27 个文件，2026-10-22 逐文件 SHA-256 比对 0 diff），**无任何修改**；另拷贝 `LICENSE`（与上游逐字节相同）
- 许可证: BSD-3-Clause（D. E. Shaw Research 2010-2012，见 LICENSE），与本项目 Apache-2.0 兼容
- 依据: RISK-027 / DESIGN §6.1（不用 submodule：避免下游 wheel 构建追加 recursive 负担）
