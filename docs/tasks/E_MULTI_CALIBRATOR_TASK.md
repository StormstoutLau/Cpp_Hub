# E. 多模型标定器 (Bates/VG/CEV) — 任务规范

## 执行站: A 站 (scott-lau-NEX.local)
## 模型: opencode/deepseek-v4-flash-free --auto

## 背景
当前 `include/cpphub/calibration/calibrator.hpp` 仅有 HestonCalibrator + SABRCalibrator.
SVI/SSVI 已有各自 calibrate() 方法 (svi.hpp/ssvi.hpp).
需新增 Bates/VG/CEV 三个标定器, 复用现有 DE + LM 优化框架.

## 文献
- Bates (1996) "Jumps and Stochastic Volatility" (Bates CF)
- Madan-Carr-Chang (1998) "The Variance Gamma Process" (VG CF)
- Cox (1975) "Notes on option pricing I: Constant elasticity of variance process" (CEV)
- Schoutens (2003) "Lévy Processes in Finance" (CF 标定参考)

## 复用资源
- `calibrator.hpp`: Calibrator 基类, HestonCalibrator/SABRCalibrator 模板
- `optimizer.hpp`: DifferentialEvolution + LevenbergMarquardt
- `objective.hpp`: ObjectiveFunction::make_iv_objective
- `bates_cf.hpp`: make_bates_cf (Bates 特征函数)
- `characteristic_functions.hpp`: make_vg_cf (VG 特征函数)
- `cev_analytic.hpp`: CEV 解析定价
- `cos_method.hpp` 或 `fourier_engine.hpp`: COS 方法定价 (用于 Bates/VG)

## 需修改文件
- `include/cpphub/calibration/calibrator.hpp` — 新增 3 个标定器类
- `tests/unit/calibration/test_calibration_framework.cpp` — 新增测试 (或独立文件 test_multi_calibrator.cpp)
- `tests/CMakeLists.txt` — 注册新测试 (若独立文件)

## 接口规范 (在 calibrator.hpp 新增)

### 1. BatesCalibrator
```cpp
// 参数向量 x = [v0, kappa, theta, sigma_v, rho, lambda, mu_J, sigma_J]
// 8 参数标定 (Heston 5 + Merton 跳跃 3)
// Objective: sum_i w_i (IV_model_i - IV_market_i)^2
// Constraints: v0>0, kappa>0, theta>0, sigma_v>0, |rho|<1,
//              lambda>0, mu_J∈R, sigma_J>0
struct BatesParams {
    Real v0, kappa, theta, sigma_v, rho;  // Heston 部分
    Real lambda, mu_J, sigma_J;            // Merton 跳跃部分
};

class BatesCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "BatesCalibrator"; }

    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 1.0},      // v0
                {1e-4, 10.0},     // kappa
                {1e-4, 1.0},       // theta
                {1e-4, 5.0},       // sigma_v
                {-0.99, 0.99},     // rho
                {1e-4, 5.0},       // lambda (跳跃强度)
                {-0.5, 0.5},       // mu_J (跳跃均值)
                {1e-4, 1.0}};      // sigma_J (跳跃波动率)
    }
    BatesParams extract_params(const std::vector<Real>& x) const {
        return BatesParams{x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]};
    }
    void set_market(Real S, Real r, Real q) { S_ = S; r_ = r; q_ = q; }

private:
    Real S_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
};
```

**定价方法**: 用 make_bates_cf + COS 方法 (cos_method.hpp) 计算 call 价格,
再用 bsm_implied_vol 反推 IV.

### 2. VGCalibrator
```cpp
// 参数向量 x = [sigma, nu, theta]
// 3 参数标定 (VG 无跳跃强度, 纯 Levy 过程)
// Objective: sum_i w_i (IV_model_i - IV_market_i)^2
struct VGParams {
    Real sigma, nu, theta;
};

class VGCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "VGCalibrator"; }

    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 2.0},      // sigma
                {1e-4, 5.0},       // nu
                {-1.0, 1.0}};      // theta
    }
    VGParams extract_params(const std::vector<Real>& x) const {
        return VGParams{x[0], x[1], x[2]};
    }
    void set_market(Real S, Real r, Real q) { S_ = S; r_ = r; q_ = q; }

private:
    Real S_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
};
```

**定价方法**: 用 make_vg_cf + COS 方法计算 call 价格, 再反推 IV.
注意 VG 的 Feller 条件: 1 - theta*nu - sigma^2*nu/2 > 0, 标定时需 clamp.

### 3. CEVCalibrator
```cpp
// 参数向量 x = [sigma, beta]
// 2 参数标定 (CEV 弹性参数 + 波动率)
// Objective: sum_i w_i (IV_model_i - IV_market_i)^2
struct CEVCalibParams {
    Real sigma, beta;
};

class CEVCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "CEVCalibrator"; }

    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 2.0},      // sigma
                {0.01, 0.99}};    // beta (0=正常, 1=对数正态, 中间为 CEV)
    }
    CEVCalibParams extract_params(const std::vector<Real>& x) const {
        return CEVCalibParams{x[0], x[1]};
    }
    void set_market(Real S, Real r, Real q) { S_ = S; r_ = r; q_ = q; }

private:
    Real S_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
};
```

**定价方法**: 用 cev_analytic.hpp 的解析定价 (非中心卡方分布), 再反推 IV.

## 测试要求 (至少 15 测试, 每标定器 5 测试)

### BatesCalibrator (5 测试)
1. **合成数据标定**: 用已知 Bates 参数生成 IV surface, 标定恢复参数 (容差 5%)
2. **Feller 条件检查**: 标定结果满足 2*kappa*theta > sigma_v^2 (或警告)
3. **跳跃参数恢复**: lambda/mu_J/sigma_J 在合成数据下可恢复 (容差 10%, 跳跃参数较难)
4. **DE + LM 优于纯 LM**: DE 初始化的标定目标函数值 ≤ 纯 LM
5. **边界处理**: 参数接近边界时不崩溃

### VGCalibrator (5 测试)
1. **合成数据标定**: 用已知 VG 参数生成 IV surface, 标定恢复参数 (容差 3%)
2. **Feller 条件**: 1 - theta*nu - sigma^2*nu/2 > 0
3. **theta=0 对称**: theta=0 时 IV smile 对称
4. **DE + LM 优于纯 LM**: DE 初始化更好
5. **边界处理**: nu→0 退化为 BS, 不崩溃

### CEVCalibrator (5 测试)
1. **合成数据标定**: 用已知 CEV 参数生成 IV, 标定恢复参数 (容差 2%)
2. **beta=1 退化为 BS**: beta→1 时 CEV IV → BS IV
3. **beta=0 正态**: beta→0 时 IV skew 为正
4. **DE + LM 优于纯 LM**: DE 初始化更好
5. **边界处理**: beta 接近 0/1 时不崩溃

## 实现要点
1. **复用 HestonCalibrator 模板**: DE 全局搜索 → LM 精炼 → 结果封装
2. **COS 方法定价**: 参考 fourier_engine.hpp 的 COSMethod 类, 或直接用 cos_method.hpp
3. **合成数据生成**: 用已知参数生成 IV surface 作为 "market quotes", 验证标定可恢复参数
4. **IV 反推**: 复用 bsm_implied_vol (已在 calibrator.hpp 中定义)
5. **参数 clamp**: 在 iv_fn lambda 中 clamp 参数到有效区间, 避免数值异常

## 验证标准
- 编译: g++ -std=c++17 -O2 (A 站 GCC 13.3.0)
- 测试: 全部通过
- 跨平台: 主站 MSVC 2022 编译通过
- 标定质量: 合成数据参数恢复误差 ≤ 5% (Bates 跳跃参数 ≤ 10%)
