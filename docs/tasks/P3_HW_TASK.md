# P3: Hull-White 完整短期利率模型扩展任务

## 任务定位
**目标**: 在现有 `include/cpphub/models/ir/short_rate.hpp` 的 HullWhite 类基础上，扩展利率衍生品定价功能，完成 IR 闭环。
**工作站**: B 站 (`scott-lau-GTR-Pro.local`)
**模型**: `opencode/deepseek-v4-flash-free` (opencode Zen 免费模型)
**工作目录**: `/tmp/oc_codetest`

## 关键约束 (来自 project_memory)
1. B 站 opencode 路径为 `~/.opencode/bin/opencode` (1.18.9)，非交互 SSH 必须用全路径
2. B 站必须 `unset http_proxy https_proxy`，否则路由到死端口 7890 卡死
3. B 站无 mihomo 代理，opencode.ai 直连可达 (HTTP 200)
4. 头文件 `#include` 必须放在 namespace 外，禁止嵌套命名空间
5. 模型选择: `opencode/deepseek-v4-flash-free` (B 站实测 100% 通过率)
6. 使用 `--auto` 参数允许 agent 自主读写

## 现有代码状态

### `include/cpphub/models/ir/short_rate.hpp` (HullWhite 类已有功能)
- `HullWhiteParams{r0, kappa, sigma}` — 参数结构
- `zero_coupon_bond(T)` — 从市场曲线插值零息债价格
- `forward_rate(T)` — 数值微分计算远期利率
- `theta(t)` — 时变漂移项 (数值微分)
- `hw_zero_coupon_bond(T)` — HW 闭式零息债
- `bond_option(T_opt, T_bond, K, is_call)` — 债券期权 (Black 形式)
- `simulate_path(T, n_steps, path, rng)` — Euler 路径模拟

### 测试文件 `tests/unit/models/test_short_rate.cpp` (已有 21 测试)
- Vasicek: 12 测试 (验证/解析解/Call-Put parity/路径模拟)
- CIR: 6 测试
- HullWhite: 7 测试 (验证/曲线匹配/theta/Bond option parity/Simulate)
- G2: 7 测试

## 需要扩展的功能 (按优先级)

### 1. Jamshidian 分解 (Jamshidian 1989) — 核心
**目的**: 把 coupon-bearing bond option 分解为 zero-coupon bond options 之和。
**算法**:
- 给定 coupon bond: 现金流 `{(t_i, c_i)}`, i=1..N
- 找 r* 使: Σ c_i * A(t,T_i) * exp(-B(t,T_i) * r*) = K (strike)
- 用 Newton 或 bisection 求解 r*
- Call on coupon bond = Σ c_i * P(0,t_i) * [P(t,T_i) 的 call on r*]

**接口**:
```cpp
// 在 HullWhite 类中新增:
// Jamshidian 分解: 求 r* 使 Σ c_i * A(t,T_i) * exp(-B(t,T_i)*r*) = K
Real jamshidian_r_star(Real T_opt,
                        const std::vector<Real>& payment_times,
                        const std::vector<Real>& cashflows,
                        Real K) const;

// Coupon-bearing bond option via Jamshidian decomposition
// cashflows[i] 在 payment_times[i] 时支付
Real coupon_bond_option(Real T_opt,
                         const std::vector<Real>& payment_times,
                         const std::vector<Real>& cashflows,
                         Real K, bool is_call) const;
```

### 2. Swaption 定价 (Payer/Receiver) — 核心
**算法**: 
- Swaption = option on swap, 用 Jamshidian 分解为 ZCB options 之和
- Payer swaption = call on fixed leg (pay fixed, receive float)
- 在 T_opt 时刻, swap value = P_float - P_fixed = (1 - K*τ) - Σ K*τ*P(T_opt,T_i)
- 找 r* 使 swap value = 0, 然后分解

**接口**:
```cpp
// Payer swaption (right to pay fixed, receive floating)
// T_opt: option expiry, T_start: swap start, T_end: swap end
// n_periods: number of fixed-leg payments
// K: strike (swap rate)
// is_payer: true=payer, false=receiver
Real swaption(Real T_opt, Real T_start, Real T_end,
              Size n_periods, Real K, bool is_payer) const;
```

### 3. Caplet / Cap 定价 — 核心
**算法**: HW 下 caplet = put on ZCB
- caplet(t_start, t_end, K_cap) = max(L(t_start, t_end) - K_cap, 0) * τ * P(t_end)
- 在 HW 下, L(t_start, t_end) = (1/τ) * (1/P(t_start,t_end) - 1)
- 化简: caplet = (1 + τ*K_cap)^{-1} * max(1 - (1+τ*K_cap)*P(t_start,t_end), 0)
        = (1 + τ*K_cap)^{-1} * put(T_start, T_end, K'=1/(1+τ*K_cap))
- Cap = Σ caplets

**接口**:
```cpp
// Caplet on LIBOR: L(t_start, t_end) vs K_cap
Real caplet(Real t_start, Real t_end, Real K_cap) const;
// Cap = Σ caplets over schedule
Real cap(const std::vector<Real>& reset_times, Real K_cap) const;
// Floor / floorlet (对称)
Real floorlet(Real t_start, Real t_end, Real K_floor) const;
Real floor(const std::vector<Real>& reset_times, Real K_floor) const;
```

### 4. HW 解析 θ(t) (可选, 增强稳定性)
当前 `theta(t)` 用数值微分计算 f'(t) + κ*f(t) + σ²/(2κ)*(1-exp(-2κt))。可以新增**解析形式**假设远期曲线分段常数:
```cpp
// 解析 theta(t) 假设远期利率在 [T_i, T_{i+1}) 上为常数
Real theta_analytic(Real t) const;
```

## 测试要求 (新增 ~15 测试, 加在 test_short_rate.cpp 末尾)

### Jamshidian 分解 (4 测试)
1. `JamshidianRStarSolvesEquation` — 把 r* 代回 Σ c_i * P(T_i) 应等于 K
2. `JamshidianRStarMonotonicInK` — K 越高, r* 越低
3. `CouponBondOptionCallPutParity` — C - P = P(0,T_opt)*(Σ c_i*P(0,T_i) - K)
4. `CouponBondOptionVsZCB` — 单一 cashflow 时应等于 bond_option

### Swaption (4 测试)
5. `SwaptionCallPutParity` — Payer + Receiver = P(0,T_start) - P(0,T_end)*K*τ*Σ
6. `SwaptionPositive` — ATM swaption 价值 > 0
7. `SwaptionMonotonicInStrike` — Payer 随 K 单调递减
8. `SwaptionDecreasesWithExpiry` — 短 expiry 的 swaption 价值较低

### Cap/Floor (4 测试)
9. `CapletNonNegative` — caplet >= 0
10. `CapMonotonicInStrike` — Cap 随 K 单调递减
11. `CapFloorParity` — Cap - Floor = swap value
12. `CapletVsZCBPut` — 验证 caplet = (1+τK)^{-1} * put 公式

### 数值稳定性 (3 测试)
13. `SwaptionFlatCurveMatchesBS` — 平坦曲线 + 特殊参数下应与 BS 类似
14. `CapletAtZeroMaturity` — t_start=t_end 时 caplet = 0
15. `JamshidianConvergence` — Newton 迭代在 20 步内收敛

## 文献参考
- Jamshidian (1989) "An exact bond option formula", Journal of Finance 44(1), 205-209
- Brigo & Mercurio (2006) "Interest Rate Models - Theory and Practice" Ch.3.3 (HW), Ch.3.11 (Jamshidian)
- Hull (2018) "Options, Futures, and Other Derivatives" 10th ed, Ch.31 (利率衍生品)
- Andersen & Piterbarg (2010) "Interest Rate Modeling" Vol I, Ch.5-6

## 实施流程

1. **更新代码**: `cd /tmp/oc_codetest && git pull`
2. **opencode 启动命令** (B 站):
```bash
unset http_proxy https_proxy
cd /tmp/oc_codetest
~/.opencode/bin/opencode run -m opencode/deepseek-v4-flash-free --auto \
  "在 include/cpphub/models/ir/short_rate.hpp 的 HullWhite 类中实现 P3 任务: (1) Jamshidian 分解 jamshidian_r_star 和 coupon_bond_option; (2) swaption (payer/receiver); (3) caplet/cap/floorlet/floor。在 tests/unit/models/test_short_rate.cpp 末尾追加 15 个测试。所有现有测试必须仍然通过。编译: cmake --build build --config Release。运行测试: ctest --test-dir build --output-on-failure -R 'ShortRate|HullWhite|G2|Vasicek|CIR'。参考 docs/tasks/P3_HW_TASK.md 中的接口规范和测试要求。"
```
3. **验收**: 测试通过, 跨平台编译成功 (GCC + MSVC), 提交 commit

## 验收标准
- [ ] 新增 4 个核心方法: jamshidian_r_star, coupon_bond_option, swaption, caplet/cap/floorlet/floor
- [ ] 新增 ≥15 个测试, 全部通过
- [ ] 全量回归测试 ≥ 989/989 通过
- [ ] 跨平台编译 (GCC + MSVC) 成功
- [ ] commit message: `feat(v1.3 P3): Hull-White complete — Jamshidian/swaption/cap-floor`
