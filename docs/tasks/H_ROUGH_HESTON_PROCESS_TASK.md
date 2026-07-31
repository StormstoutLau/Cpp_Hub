# H. Rough Heston 过程层 (分数阶 Euler 采样 + MC 定价) — 任务规范

## 执行站: B 站 (scott-lau-GTR-Pro.local)
## 模型: opencode/deepseek-v4-flash-free --auto

## 背景
Rough Heston (El Euch-Rosenbaum 2018) 方差过程服从 fractional CIR:
  D^α v(t) = κ(θ - v(t)) + σ√v(t) dW(t), α = H + 1/2 ∈ (0.5, 1)
价格过程: dS/S = (r-q)dt + √v(t) dZ(t), dZ = ρdW + √(1-ρ²)dW⊥

本任务实现过程层 (路径采样 + MC 定价), 与 A 站解析层 (CF) 交叉验证.

## 文献
- El Euch, Rosenbaum (2018) "The characteristic function of rough Heston model"
- Bennedsen, Lunde, Pakkanen (2017) "Hybrid scheme for Brownian semistationary processes" (采样参考)
- Smith (2017) "Option pricing under the rough Heston model" (Euler 方案)

## 复用资源
- `models/diffusion/process.hpp`: StochasticProcess 基类 (generate_path 接口)
- `models/diffusion/rough_bergomi.hpp`: RoughBergomiProcess (Volterra 核采样参考)
- `models/diffusion/rbergomi_hybrid_scheme.hpp`: RLFbmHybridSampler (近似核采样参考)
- `pricing/analytic/heston_cf.hpp`: Heston CF (H=0.5 退化基准)
- `pricing/fourier/cos_method.hpp`: COSEngine (MC vs COS 交叉验证, 依赖 A 站 CF)
- `core/rng.hpp`: Philox4x64, box_muller

## 需新增文件
1. `include/cpphub/models/diffusion/rough_heston.hpp` — Rough Heston 过程类
2. `tests/unit/models/test_rough_heston_process.cpp` — 过程层单元测试
3. `tests/CMakeLists.txt` — 注册测试

## 数学规范

### 1. 分数阶 CIR 积分形式 (Volterra 型)
v(t) = v(0) + (1/Γ(α)) ∫₀^t (t-s)^{α-1} [κ(θ - v(s)) ds + σ√v(s) dW(s)]

### 2. 分数阶 Euler 离散化 (Smith 2017)
网格 t_j = j·Δt, j=0,...,N, Δt = T/N.

确定性核积分:
  K_{j,i} = ∫_{tᵢ}^{tᵢ₊₁} (t_{j+1}-s)^{α-1} ds = [(t_{j+1}-tᵢ)^α - (t_{j+1}-tᵢ₊₁)^α] / α

方差更新 (Full Truncation 避免负方差):
  v_{j+1} = v(0) + (κθ/Γ(α)) · (t_{j+1})^α
            - (κ/Γ(α)) · Σᵢ₌₀ʲ K_{j,i} · vᵢ⁺
            + (σ/Γ(α)) · Σᵢ₌₀ʲ √(K_{j,i} · vᵢ⁺) · Zᵢ

其中 vᵢ⁺ = max(vᵢ, 0) (full truncation), Zᵢ ~ N(0,1) i.i.d.

注: 随机积分项近似为 √(K_{j,i}) · √vᵢ · Zᵢ (Itô 等距近似).

### 3. 价格过程 (Euler)
S_{j+1} = S_j · exp((r - q - 0.5·v_j⁺)·Δt + √(v_j⁺·Δt) · dZ_j)
dZ_j = ρ·dW_j + √(1-ρ²)·dW⊥_j

其中 dW_j 与方差过程的 Zᵢ 相关 (杠杆效应).
简化: 方差过程的 Z_j 即为 dW_j/√Δt, 价格过程用同一个 W_j.

### 4. H=0.5 (α=1) 退化
α=1 时 K_{j,i} = ∫_{tᵢ}^{tᵢ₊₁} 1 ds = Δt (常数),
v_{j+1} = v(0) + κθ·t_{j+1} - κ·Σᵢ vᵢ·Δt + σ·Σᵢ √(vᵢ·Δt)·Zᵢ
即标准 CIR Euler, 与 Heston 过程 (Euler scheme) 一致.

## 接口规范

```cpp
#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/models/diffusion/process.hpp"
#include "cpphub/core/rng.hpp"
#include <vector>

namespace cpphub {
inline namespace v1 {

struct RoughHestonParams {
    Real H;       // Hurst ∈ (0, 0.5), α = H + 0.5
    Real kappa;   // 均值回归速度
    Real theta;   // 长期方差
    Real sigma;   // vol of vol
    Real rho;     // 相关性 ∈ (-1, 1)
    Real v0;      // 初始方差
    Real S0;      // 标的现价
    Real r;       // 利率
    Real q;       // 股息率
};

// 预计算 Volterra 核系数 K_{j,i} (上三角, j≥i)
// K[j][i] = [(t_{j+1}-t_i)^α - (t_{j+1}-t_{i+1})^α] / α
// 返回 N×N 矩阵 (j=0..N-1, i=0..j)
std::vector<std::vector<Real>> rough_heston_kernel(Real T, Size n_steps, Real alpha);

// Rough Heston 过程 (fractional CIR + GBM 价格)
class RoughHestonProcess : public StochasticProcess {
public:
    RoughHestonProcess(const RoughHestonParams& params, Real T, Size n_steps);

    Size dimension() const override { return 2; }  // [S, v]
    Real spot() const override { return params_.S0; }

    // 生成一条路径: path[0]=S 路径, path[1]=v 路径
    // path[asset_idx][step_idx], step_idx=0..n_steps
    void generate_path(Real T, Size n_steps, std::span<Real> path, Philox4x64& rng) const override;

    // 便捷接口: 返回 [S_path, v_path], 各为 n_steps+1 长度
    std::vector<std::vector<Real>> generate_path(Philox4x64& rng, Real dt_scale = 1.0) const;

    const RoughHestonParams& params() const { return params_; }
    Real T() const { return T_; }
    Size n_steps() const { return n_steps_; }
    Real alpha() const { return params_.H + 0.5; }

private:
    RoughHestonParams params_;
    Real T_;
    Size n_steps_;
    std::vector<Real> t_grid_;
    std::vector<std::vector<Real>> K_;  // 预计算核
};

// ============ MC 定价 ============
struct RoughHestonMCResult {
    Real price;
    Real std_error;
    Real ci_lower;
    Real ci_upper;
    Size n_paths;
};

// 欧式 call/put MC 定价
RoughHestonMCResult rough_heston_price_european(
    const RoughHestonParams& params,
    Real T, Real K, bool is_call,
    Size n_steps, Size n_paths,
    uint64_t seed);

}  // namespace v1
}  // namespace cpphub
```

## 测试要求 (至少 10 测试)

1. **核矩阵正确性**: K[j][j] = [(Δt)^α - 0]/α = Δt^α/α (对角元), K[j][i]=0 for i>j
2. **核积分总和**: Σᵢ K_{N-1,i} = T^α/α (总积分 = ∫₀^T (T-s)^{α-1}ds = T^α/α)
3. **α=1 退化为常数核**: α=1 时 K_{j,i} = Δt (所有元素相等)
4. **方差非负**: Full truncation 保证 v_j ≥ 0 对所有 j (1000 路径验证)
5. **H=0.5 退化为 Heston Euler**: α=1 时 Rough Heston MC 价格与 Heston Euler MC 一致 (容差 2×std_error, n_paths≥10000)
6. **H=0.5 MC vs Heston CF**: α=1 时 MC 价格与 Heston COS 价格一致 (容差 5%, n_paths≥50000)
7. **杠杆效应**: rho<0 时 OTM put 的 IV > OTM call 的 IV (skew)
8. **路径确定性**: 相同 seed 生成相同路径 (可复现性)
9. **价格鞅性**: E[S_T] = S₀·exp((r-q)·T) (容差 1%, n_paths≥50000)
10. **数值稳定性**: H=0.05 (极端粗糙) + 长期限 (T=2y) 不产生 NaN/Inf

## 实现要点

### 1. 核预计算 (关键性能优化)
- K_{j,i} 只依赖 α, T, n_steps, 可在构造时预计算
- 存储为 N×N 下三角矩阵 (j≥i)
- 复杂度 O(N²) 一次, 路径采样 O(N²) per path
- N=200 时每路径约 0.04ms, 50000 路径约 2s

### 2. 方差更新 (Full Truncation)
```
for j = 0 to N-1:
    v_pos = max(v_j, 0)
    drift_sum = 0, diff_sum = 0
    for i = 0 to j:
        drift_sum += K[j][i] * max(v_i, 0)
        diff_sum += sqrt(K[j][i] * max(v_i, 0)) * Z_i
    v_{j+1} = v0 + (kappa*theta/Gamma(alpha)) * t_{j+1}^alpha
              - (kappa/Gamma(alpha)) * drift_sum
              + (sigma/Gamma(alpha)) * diff_sum
```

### 3. 价格过程 (相关布朗运动)
```
for j = 0 to N-1:
    dW = sqrt(dt) * Z_j       // 方差过程噪声
    dZ = rho * dW + sqrt(1-rho^2) * sqrt(dt) * Z_perp  // 价格过程噪声
    S_{j+1} = S_j * exp((r - q - 0.5*max(v_j,0))*dt + sqrt(max(v_j,0)) * dZ)
```
注: Z_j 与方差更新中的 Z_i 是同一个序列 (杠杆相关).

### 4. H=0.5 退化验证
- α=1 时核退化为常数 Δt, 方差更新退化为标准 CIR Euler
- 与 `models/diffusion/heston.hpp` 的 Euler scheme 对比
- 容差: MC 噪声, 用 2×std_error 或 5% 容差

### 5. MC vs COS 交叉验证 (需 A 站 CF)
- 等 A 站完成后, 用 make_rough_heston_cf + COSEngine 计算 call 价格
- MC 价格 (n_paths=50000) vs COS 价格, 容差 5%
- H<0.5 时验证 rough volatility 的 IV skew 比 Heston 更陡

## 验证标准
- 编译: g++ -std=c++17 -O2 (B 站 GCC 13.3.0)
- 测试: 全部通过
- 跨平台: 主站 MSVC 2022 编译通过
- 精度: H=0.5 退化容差 5%, 鞅性容差 1%

## 依赖
- 本任务可与 A 站 (G任务, CF) 并行开发
- H=0.5 退化测试用现有 Heston 基础设施, 不依赖 A 站
- MC vs COS 交叉验证 (测试 6) 依赖 A 站 CF, 若 A 站未完成可跳过此测试
