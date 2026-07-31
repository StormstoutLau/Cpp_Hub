# D. CVA with Wrong-Way Risk (WWR) — 任务规范

## 执行站: B 站 (scott-lau-GTR-Pro.local)
## 模型: opencode/deepseek-v4-flash-free --auto

## 背景
当前 `include/cpphub/risk/xva.hpp` 的 CVA 假设违约与暴露独立:
  CVA_independent = -(1-R_c) * ∫ EPE_disc(t) * dPD_c(t)

现实场景中, 交易对手违约往往发生在暴露上升时 (Wrong-Way Risk, WWR),
例如: 出售看跌期权 → 标的下跌 → 交易对手违约 → 暴露同时上升.
需建模违约-暴露相关性, 计算 WWR-adjusted CVA.

## 文献
- Pykhtin & Zhu (2007) "A Guide to Modeling Counterparty Credit Risk" (Gaussian Copula WWR)
- Hull & White (2012) "CVA and Wrong Way Risk" (近似公式)
- BCBS (2015) "Margin requirements for non-centrally cleared derivatives"
- Gregory (2015) "XVA" Ch.15-17

## 需新增文件
- `include/cpphub/risk/wrong_way_risk.hpp` — WWR CVA 模型
- `tests/unit/risk/test_wrong_way_risk.cpp` — 单元测试
- `tests/CMakeLists.txt` — 注册测试 (用 cpphub_add_test 宏)

## 接口规范 (wrong_way_risk.hpp)

```cpp
#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include "cpphub/risk/xva.hpp"
#include <vector>

namespace cpphub {
inline namespace v1 {

// WWR 配置
struct WWRConfig {
    Real rho_ww = 0.0;  // 暴露-违约相关性 ∈ [-0.99, 0.99]
                        // rho_ww > 0: WWR (违约时暴露高)
                        // rho_ww < 0: RWR (Right-Way Risk, 违约时暴露低)
    // Gaussian copula 参数化:
    // 暴露驱动 M_V ~ N(0,1), 违约驱动 M_D ~ N(0,1), corr(M_V, M_D) = rho_ww
    // V(t) = f(M_V, t), default if M_D > Φ^{-1}(1-PD(t))
};

// ============ 方法 1: Hull-White WWR 近似 (半解析) ============
// Hull-White (2012): 用 WWR 调整因子近似
// CVA_wwr ≈ CVA_independent * adjustment_factor(rho_ww, PD, sigma_V)
// adjustment_factor = exp(rho_ww * sigma_V * sqrt(T) * Φ^{-1}(PD) + 0.5*(rho_ww*sigma_V*sqrt(T))^2)
// 其中 sigma_V 为暴露 V(t) 的波动率
//
// 输入:
//   profile: 独立假设下的暴露轮廓 (已计算)
//   pd_counterparty: 交易对手 PD 曲线
//   cfg: WWR 配置 (rho_ww)
//   sigma_V: 暴露 V(t) 的波动率 (年化, 可从 MC 路径估计)
// 输出: WWR 调整后的 XVAResult (CVA 放大, DVA/FVA 不变)
XVAResult compute_cva_wwr_hw(
    const ExposureProfile& profile,
    const PDCurve& pd_counterparty,
    const PDCurve& pd_self,
    const XVAConfig& xva_cfg,
    const WWRConfig& wwr_cfg,
    Real sigma_V,
    Real risk_free_price = 0.0);

// ============ 方法 2: Gaussian Copula WWR (Pykhtin-Zhu 2007, 半解析) ============
// 条件期望: E[max(V(t),0) | M_D = m]
// V(t) 假设对数正态: ln V(t) ~ N(mu_V(t), sigma_V(t)^2)
// 条件化: ln V(t) | M_D=m ~ N(mu_V + rho*sigma_V*m, sigma_V^2*(1-rho^2))
// E[max(V,0)|M_D=m] = V0 * exp(mu_cond + 0.5*sigma_cond^2) * Phi(d1_cond)
// 其中 d1_cond = (mu_cond + sigma_cond^2 - ln K_adj) / sigma_cond
//
// WWR CVA = -(1-R) * ∫_0^T E[max(V(t),0) * 1_{τ∈(t,t+dt)}] dt
//         = -(1-R) * ∫_0^T [∫ E[max(V,0)|M_D=m] * f_{M_D}(m) dm] * dt
// 数值实现: 在每个 t_i, 对 M_D 做高斯求积 (Gauss-Hermite)
XVAResult compute_cva_wwr_copula(
    const ExposureProfile& profile,
    const PDCurve& pd_counterparty,
    const PDCurve& pd_self,
    const XVAConfig& xva_cfg,
    const WWRConfig& wwr_cfg,
    Real sigma_V,
    Real risk_free_price = 0.0);

// ============ 方法 3: MC WWR 模拟 ============
// 路径生成: 标的资产 GBM + 违约时间通过 copula 相关
// 每条路径:
//   1. 生成 M_V, M_D ~ N(0,1) with corr = rho_ww
//   2. V(t) 沿 GBM 路径 (用 M_V 驱动)
//   3. 违约时间 τ = inf{t : M_D > Φ^{-1}(1-PD(t))} (单因子模型)
//   4. 若 τ ≤ T, 损失 = (1-R) * max(V(τ), 0) * P_d(0, τ)
// CVA_wwr = -E[损失]
//
// 输入: value_fn (同 xva.hpp), gen (路径生成器), n_paths, seed
Real compute_cva_wwr_mc(
    const MultiAssetGBMPathGenerator& gen,
    std::function<Real(Real, const std::vector<Real>&)> value_fn,
    const std::vector<Real>& exposure_times,
    const PDCurve& pd_counterparty,
    const ZeroCurve& discount_curve,
    const XVAConfig& xva_cfg,
    const WWRConfig& wwr_cfg,
    Real sigma_V,
    Size n_paths,
    uint64_t seed);

// ============ 辅助: 从 MC 路径估计暴露波动率 ============
Real estimate_exposure_volatility(
    const std::vector<std::vector<Real>>& V_samples,
    const std::vector<Real>& times);

}  // namespace v1
}  // namespace cpphub
```

## 测试要求 (test_wrong_way_risk.cpp, 至少 12 测试)

1. **rho_ww=0 退化**: WWR CVA (HW 近似) == 独立 CVA (容差 1e-10)
2. **rho_ww>0 放大**: WWR CVA > 独立 CVA (WWR 增加损失)
3. **rho_ww<0 缩小**: RWR CVA < 独立 CVA (Right-Way Risk)
4. **rho_ww 单调性**: CVA 随 rho_ww 单调递增
5. **极端 rho_ww→1**: CVA 显著放大 (≥ 2x 独立 CVA)
6. **PD=0 边界**: WWR CVA = 0
7. **R=1 边界**: WWR CVA = 0 (无损失)
8. **Copula vs HW 近似一致性**: 两种方法在 rho_ww 小时 (≤0.3) 结果接近 (容差 10%)
9. **MC WWR vs Copula 一致性**: MC 结果与半解析 Copula 一致 (容差 5%, n_paths≥10000)
10. **MC rho_ww=0 退化为独立**: MC WWR 在 rho_ww=0 时 == 独立 MC CVA
11. **暴露波动率估计**: estimate_exposure_volatility 在已知 GBM 下匹配解析值
12. **数值稳定性**: 大 rho_ww (0.9) + 长期限 (10y) 不产生 NaN/Inf

## 实现要点
1. 复用 `xva.hpp` 的 ExposureProfile / XVAConfig / XVAResult / compute_xva
2. 复用 `credit_curve.hpp` 的 PDCurve
3. 复用 `ois_curve.hpp` 的 ZeroCurve
4. 复用 `multi_asset_path_generator.hpp` 的 MultiAssetGBMPathGenerator
5. Gauss-Hermite 求积: 可用 n=20 节点 (查表或 Newton 迭代)
6. MC WWR 路径: 生成相关高斯 [M_V, M_D], M_V 驱动资产价格, M_D 决定违约时间
7. 违约时间映射: τ = inf{t : PD(t) > Φ(M_D)} (单因子 Gaussian copula)

## 验证标准
- 编译: g++ -std=c++17 -O2 (B 站 GCC 13.3.0)
- 测试: 全部通过
- 跨平台: 主站 MSVC 2022 编译通过
