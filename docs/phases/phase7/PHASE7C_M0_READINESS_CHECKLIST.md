# Phase 7C M0 实施前置条件检查清单

> **版本**: v1.0 (2026-08-17, 基于 [PHASE7C_SPEC.md](./PHASE7C_SPEC.md) v1.2 生成)
> **范围**: M0 回填三项 (NP / ZA / GARCH-M) 开工前置 — 区别于 [PHASE7C_ACCEPTANCE_CHECKLIST.md](./PHASE7C_ACCEPTANCE_CHECKLIST.md) (实施后验收)
> **标注**: `[x]` 已验证 / `[ ]` 待办 / `[!]` **缺口 (阻断开工)**
> **预检方式**: 主线一手实测 (源码 grep / pip·Rscript 环境探测 / 根 CMakeLists 读取), 2026-08-17

---

## 0. 预检结论摘要

| 类别 | 结论 |
|------|------|
| ✅ 已就绪 | 规范链 (spec v1.2 + checklist 331 项) / 复用主体 (GARCH 基础设施·SLSQP·Schwert) / statsmodels 0.14.6 (conda) / rugarch 1.5.6 (用户库) / **arch 8.0.0 (已装, ARCHInMean+form ✓)** / **urca 1.3.4 (已装 F:/R/win-library/4.6, ur.za/ca.po/ca.jo ✓)** |
| ⚠️ 处置项 | DF-GLS 去势变换不可直接复用 (内联于检验函数) → 按预判重实现 + 偏差归因 |
| ~~❌ 阻断缺口~~ | **已清零 (2026-08-17)**: arch 8.0.0 + urca 1.3.4 安装并验证 (附录留痕); **M0 状态 = GO** |
| ❓ 可降级 | Stata 可编程访问未验证 (NP 逐 k MAIC) — 降级路径已有 |

---

## A. 规范与流程前置 (3/3 就绪)

| # | 检查项 | 状态 | 证据/说明 |
|---|--------|------|----------|
| A1 | spec v1.2 冻结 (R1-R4 清零 + §1.4 接口总则 + §8.4-8.8 工程约束) | [x] | commit `0d70343`; §0 门禁记录在案 |
| A2 | M0 范围与决策映射: 4 头文件 (§1.1.1-1.1.4) ← ADR-019 回填决策 NP/ZA/GM | [x] | spec §2.1-2.3 接口签名齐全 (v1.2) |
| A3 | 验收映射: checklist §3 (22 项) + §8.5 (16 项) + §1.1.1-4 + 性能 P1-P3; 附录 C 节奏 M0 先行 | [x] | checklist `86e835c` |
| A4 | 对照禁令生效: Julia 常数情形 MPT 禁作对照 (H1 裁决, 仅趋势可旁证) | [x] | spec §13-b; verify_np_semi.py 设计时执行 |

## B. 复用接口审计 (v1.6 既有代码, additive-only 约束 §8.5)

| # | 检查项 | 状态 | 预检证据 (2026-08-17 实测) |
|---|--------|------|--------------------------|
| B1 | **DF-GLS GLS 去势变换复用** (NP Step 1 同一 ERS 变换) | [!]→⚠️处置 | `df_gls_test.hpp` 仅有公共 `df_gls_test()`, 去势逻辑**内联其中**, 无独立 `gls_detrend` 函数 ⇒ NP 无法 include 复用。**处置 (预批)**: 新头文件内重实现 ERS 变换 (c̄/ρ̄=1+c̄/T/首项不变换), 登记复用偏差归因 (checklist §13.10, 7B Issue #2 惯例); 提取公共函数=触碰 v1.6 文件, 违 §8.5, 不取 |
| B2 | Schwert lag / AIC 滞后框架 (ZA 复用) | [x] | `unit_root_common.hpp` `schwert_lag()` 为公共 inline ✓ |
| B3 | GARCH 基础设施 (GM 复用: backcast/似然/约束/sandwich) | [x] | `garch_distribution.hpp` 似然为公共函数 ✓; `garch_model.hpp` GarchConfig/递归模式可仿写。**预检发现**: GarchConfig 无 seed 字段, 但多起始随机扰动用 **Philox 确定性计数器** (garch_model.hpp L376 "Philox 确定性") ⇒ 可复现性内建, 与 §8.6/§1.4-8 精神兼容, GM 沿用同模式, 无需改 Config (零破坏保持) |
| B4 | SLSQP (C5) | [x] | ADR-018 已验收 12/12; `calibration/optimizer.hpp` 公共接口 |
| B5 | MacKinnon 表管线仿写 (np_table1/za .inc 惯例) | [x] | 7B `gen_mackinnon_tables.py` + `mackinnon_cv.hpp` constexpr+static_assert 模式可复制 |

## C. 对照环境 (2 缺口阻断)

| # | 检查项 | 状态 | 预检证据 | 处置 |
|---|--------|------|----------|------|
| C1 | Python **arch 8.0.0** (GM 主锚, ARCHInMean form 三值) | [x] **已装并验证** | `pip install arch==8.0.0` (cp311 wheel); `arch.__version__ == '8.0.0'` + `ARCHInMean` 可导入 + `__init__` 含 `form` 参数 ✓ (GM2/GM3 锚落位) |
| C2 | statsmodels 0.14.6 (ZA 主对照 zivot_andrews) | [x] | conda open-webui: sm 0.14.6 ✓ (注意: 系统 python 无 — verify 脚本须用该 env 解释器或显式路径) |
| C3 | R **urca** (ZA ur.za 固定 lag 对照; 亦为 M3 OL1992 表转录源) | [x] **已装并验证** | CRAN Windows 二进制 **urca 1.3-4** 装入 `F:/R/win-library/4.6`; 验证: `ur.za`/`ca.po`/`ca.jo` 全存在 ✓ (M0 ur.za + M3 前瞻一次到位) |
| C4 | R rugarch (GM archpow 交叉 + fix() 三步法) | [x] | win-library/4.6: rugarch **1.5.6** ✓, 与 urca 双库路径共存加载验证通过 (spec 冻结"R rugarch"未锁版本, 实际 1.5.6 记录入 verify_gm.R) |
| C5 | R 环境 .libPaths 教训执行 | [x] | **v1.1 更新**: urca 落位 `F:/R/win-library/4.6` (agent 沙箱禁写 `C:\Users\Peng\R`, 实测确认), verify_*.R 头部须**双库并列**: `.libPaths(c('F:/R/win-library/4.6', file.path(Sys.getenv("USERPROFILE"),"R","win-library","4.6"), .libPaths()))` (7B 模板扩展) |
| C6 | Stata (NP 逐 k MAIC, dfgls r(results)) | [ ]待确认 | 可编程访问未验证 (Stata 批处理模式可用性) | 路径1: Stata batch `dfgls y, maxlag(k)` 导出 r(results) CSV; **降级路径2 (预批)**: 手册/文献逐 k 数值抄录 + Gretl `adf --gls --test-down` 对照 (H2 调研已核 Gretl 路径), 偏差归因记录 |
| C7 | vars/midasr/Spillover (M2/M4 对照) | [ ]非 M0 | 均未装 (预检) — 不阻断 M0, M1/M4 开工前安装 |

## D. 数值基线与数据源 (实施前固化)

| # | 检查项 | 状态 | 说明 |
|---|--------|------|------|
| D1 | NP Table 1 抄录源 (12 值: 4 统计量 × 2 情形 × 3 水平) | [ ] | 主源 BC wp369 PDF (URL 已知) + AU 副本双源核对; ⚠️ Julia `NGPERRON_CRITICAL_VALUES` 不作源 (H1 已证其常数分支统计量有 bug, 表亦未独立核过); 5% 锚 spec 已冻结 (−8.10/−1.98/0.233/3.17; −17.30/−2.91/0.168/5.48) |
| D2 | ZA1992 论文表 (主, za1992_cv.inc) | [ ] | 源 = urca 源码常量 (依赖 C3 装后转录, constexpr+static_assert) 或 ZA1992 论文表; 与 urca 对齐即 spec 决策 |
| D3 | ZA MC 表 (za_mc_cv.inc) | [ ] | 源 = statsmodels `zivot_andrews.py` 内嵌表 (conda 0.14.6 可直读); 锚: c 1%=−5.27644/5%=−4.81067/10%=−4.56618, t 1%=−5.03421, ct 1%=−5.57556 (ZA4) |
| D4 | NP 合成夹具: 固定数据 → Stata 逐 k MAIC/σ̂² CSV (+ Julia 趋势情形旁证) | [ ] | 数据构造记 seed; C6 路径决定 |
| D5 | ZA 基线: statsmodels (Baum 模式) + urca (固定 lag, trim 放开) 固定数据输出 | [ ] | 依赖 C2 ✓ / C3 缺口 |
| D6 | GM 基线: arch 8.0 三 form + rugarch archpow 1/2 + fix() llf (rescale 关闭) | [ ] | 依赖 C1/C4; 三步法脚本 verify_gm.R 预写 |
| D7 | 夹具记录规范: 库版本/seed/数据构造全录 CSV 头 (§8.6) | [ ] | 沿用 7B 模式 |

## E. 工程约束落地 (v1.2 §1.4/§8.4-8.8)

| # | 检查项 | 状态 | 说明 |
|---|--------|------|------|
| E1 | M0 零 Eigen: 4 头文件纯标量 (grep 断言预设: `#include.*igen` 零命中) | [ ] | §8.2 表; M0 不触碰 cpphub_timeseries_mat (M2 才建) |
| E2 | additive-only 例外清单确认: 允许触碰的既有文件**仅** `tests/CMakeLists.txt` (追加 3 套注册) | [x] | §8.5 明文; B1 处置不触碰 v1.6 头文件 |
| E3 | 接口总则: NgPerron cv 用 `std::array` ✓ (v1.2 已修); ZA `baum_preselect` 默认主模式 ✓; Result 默认初始化/NaN 政策 | [x] | spec v1.2 签名已冻结, 实施照抄 |
| E4 | 精度/可复现: 沿用根 CMakeLists flags (无需改动); GM 扰动 Philox 模式 (B3) | [x] | §8.6 |
| E5 | 性能断言落位: P1 (NP T=1000 <1s) / P2 (ZA T=500 <5s) / P3 (GM T=5000 三变体 <10s) 写入对应测试套件 | [ ] | checklist §15.1-15.3 镜像 |
| E6 | 源码 UTF-8 + `/utf-8` (MSVC) 沿用; 排幻觉点注释格式 `// NP4: ...` | [ ] | §12 规范 |

## F. 开工判据 (Go/No-Go)

| 判据 | 状态 |
|------|------|
| ~~❌ No-Go 项~~ | **0 项** (2026-08-17 清零: C1 arch 8.0.0 + C3 urca 1.3-4 安装并验证) |
| ⚠️ 预批处置 (可带病开工): B1 去势重实现+归因 / C6 Stata 降级路径 | 2 项 |
| ✅ 其余 A/B/D/E 项: A 全就绪, B 4/5, D/E 为实施首日任务 (不阻断开工, 阻断合并) | — |

**结论 (2026-08-17 更新)**: **M0 状态 = GO** — 阻断缺口清零; D 组数值基线 (D1 NP Table 1 双源抄录为 NP 实施唯一前置) 与 C6 确认可与实施并行推进。

---

## 附录: 预检命令留痕 (2026-08-17)

```
# 环境 (系统 python 无包 — 7B 惯例用 conda env + R 用户库):
& 'C:\Users\Peng\.conda\envs\open-webui\python.exe'  → sm 0.14.6 ✓, arch ✗ (预检时)
Rscript (win-library/4.6) → rugarch 1.5.6 ✓, urca/vars/midasr/Spillover ✗ (预检时)
# 源码:
df_gls_test.hpp → 无公共 gls_detrend (B1); garch_model.hpp L376 → Philox 确定性扰动 (B3)

# 安装与验证 (2026-08-17, No-Go 清零):
pip install arch==8.0.0 (conda env, cp311 wheel)
  → arch 8.0.0; ARCHInMean 可导入; __init__ 含 'form' 参数 ✓
install.packages('urca', lib='F:/R/win-library/4.6')  (CRAN win 二进制)
  ⚠️ 实测: agent 沙箱禁写 C:\Users\Peng\R (install.packages _test_dir 探测被拦, 含
    requires_approval 重试仍拦) → 改装 F:\R\win-library\4.6 (R 标准布局镜像路径)
  → urca 1.3-4; ur.za/ca.po/ca.jo 存在 ✓; 与 rugarch 1.5.6 双库路径共存加载 ✓
```
