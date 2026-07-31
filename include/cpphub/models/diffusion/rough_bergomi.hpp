#pragma once
// SOURCE: Bayer, Friz, Gatheral (2016) "Pricing under rough volatility"
//         Quantitative Finance, 16(6), 887-904.
// SOURCE: Bennedsen, Lunde, Pakkanen (2017) "Hybrid scheme for modeling
//         rough fractional Brownian motion".
//         Decisions in Economics and Finance, 40(1), 393-419.
// SOURCE: Gatheral, Jaisson, Rosenbaum (2018) "Volatility is rough".
//         Springer.
//
// 模块: Rough Bergomi (rBergomi) 模型 — 粗糙波动率模型
//
// ==================== rBergomi 模型数学 ====================
//
// 价格过程 (风险中性):
//   dS_t / S_t = sqrt(v_t) * dZ_t,  S_0 给定
//   Z 与 W (方差驱动 Brownian) 相关 ρ:  Z = ρ*W + sqrt(1-ρ²)*W^⊥
//
// 瞬时方差过程 (非马尔可夫, 路径依赖):
//   v_t = ξ_0(t) * exp( η * W̃^H_t - 0.5 * η² * t^{2H} )
//
// 其中:
//   - H ∈ (0, 1/2): 粗糙参数 (Hurst 指数, 越小越粗糙, 典型 0.05-0.15)
//   - η > 0: 波动率波动率 (vol of vol)
//   - ξ_0(t): 初始远期方差曲线 (t=0 时市场对 t 时刻方差的预期, 常数则为平期)
//   - ρ ∈ [-1, 1]: 价格-方差相关性 (典型负值, lever effect)
//   - W̃^H_t: Riemann-Liouville 分数 Brownian 运动 (RL-fBm):
//       W̃^H_t = sqrt(2H + 1) * ∫_0^t (t-s)^{H - 1/2} dW_s
//
// RL-fBm 与标准 fBm 的区别:
//   - 标准 fBm B^H_t 有平稳增量, 但 RL-fBm W̃^H_t = ∫ K(t,s) dW_s 不平稳
//   - rBergomi 用 RL-fBm (非平稳), 使 v_t 在 t=0 时确定, 演化更自然
//
// 关键性质:
//   - v_t 是 log-normal (因 W̃^H_t 是高斯)
//   - log(v_t) 的方差: Var(log v_t) = η² * t^{2H} (注意: 不是 t, 是 t^{2H})
//   - 当 H → 1/2 时, W̃^H → W (标准 Brownian), rBergomi → Bergomi (标准随机波动率)
//   - H 小 (如 0.1): 方差过程路径极粗糙, 对应实际市场波动率观察
//
// ==================== fBm 模拟: Hybrid Scheme ====================
//
// RL-fBm 的精确模拟需要 O(N²) 的 Cholesky 分解 (N = 路径步数).
// Hybrid Scheme (Bennedsen 2017) 将积分核分为两部分:
//
//   W̃^H_t_i = sqrt(2H+1) * [ ∫_0^{t_i} (t_i - s)^{H-1/2} dW_s ]
//                      ≈ ∑_{k=0}^{i-1} (t_i - t_k)^{H-1/2} ΔW_k  (Riemann 和)
//
// 对核 (t-s)^{H-1/2} 在 s 远离 t 时近似平坦 (H-1/2 > -1/2),
// 远端可用大步长近似; 近端 (s 接近 t) 核奇异, 需精确处理.
//
// Hybrid 划分 (参数 b 为分界点, 典型 b=1):
//   - 近端 (i-k ≤ b): 精确核 (t_i - t_k)^{H-1/2}, 用 Cholesky 模拟前 b+1 步
//   - 远端 (i-k > b): 用平顶近似 (t_i - t_k)^{H-1/2} ≈ (t_{i-b} - t_k)^{H-1/2}
//     即用第 (i-b) 步的核值, 配合 Euler 步进
//
// 实际实现 (本文件采用简化版 Cholesky 精确法, N ≤ 256 可接受):
//   完整 Cholesky: 构造协方差矩阵 Σ_{ij} = Cov(W̃^H_{t_i}, W̃^H_{t_j}),
//   Σ_{ij} = (2H+1) / (H + 1/2)² * ∫_0^{min(t_i,t_j)} (t_i-s)^{H-1/2} (t_j-s)^{H-1/2} ds
//   积分有闭式 (用 Beta 函数). 然后 L = chol(Σ), W̃^H = L * Z.
//
// 复杂度: Cholesky O(N³), 对 N=100 步路径约 1ms, 可接受.
// 对大 N (>500) 应切到 Hybrid Scheme, 但研究用途 N≤256 足够.

#include <cmath>
#include <complex>
#include <vector>
#include <span>
#include <stdexcept>
#include <algorithm>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/process.hpp"

namespace cpphub {
inline namespace v1 {

// ============ rBergomi 参数结构 ============
struct RoughBergomiParams {
    Real H;       // 粗糙参数 (Hurst 指数), (0, 0.5)
    Real eta;     // 波动率波动率 (vol of vol) > 0
    Real rho;     // 价格-方差相关性 [-1, 1]
    Real xi0;     // 初始远期方差 (常数曲线, > 0). 完整实现应支持 ξ_0(t) 函数
    Real S0;      // 初始价格 > 0
    Real r;       // 无风险利率
    Real q;       // 股息率
};

inline void validate_rough_bergomi_params(const RoughBergomiParams& p) {
    if (p.H <= 0.0 || p.H >= 0.5)
        throw std::invalid_argument("rBergomi: H must be in (0, 0.5)");
    if (p.eta <= 0.0)
        throw std::invalid_argument("rBergomi: eta must be positive");
    if (p.rho < -1.0 || p.rho > 1.0)
        throw std::invalid_argument("rBergomi: rho must be in [-1, 1]");
    if (p.xi0 <= 0.0)
        throw std::invalid_argument("rBergomi: xi0 must be positive");
    if (p.S0 <= 0.0)
        throw std::invalid_argument("rBergomi: S0 must be positive");
}

// ============ RL-fBm 协方差矩阵计算 ============
// Cov(W̃^H_{t_i}, W̃^H_{t_j}) = (2H+1) / (H+1/2)² * I(t_i, t_j)
// 其中 I(a, b) = ∫_0^{min(a,b)} (a-s)^{H-1/2} (b-s)^{H-1/2} ds
//
// 闭式解 (令 m = min(a,b), M = max(a,b), α = H+1/2):
//   I(a, b) = (m^{2H+1} / (2H+1)) * 2F1(-α, 1; 2H+2; m/M)
//   其中 2F1 为超几何函数.
//
// 简化: 对均匀网格 t_i = i*dt, 用数值积分 (Gauss-Legendre) 避免超几何函数.
// 对角线 (i=j): I(t,t) = t^{2H+1} / (2H+1)
//
// 更简洁的形式 (利用 RL-fBm 的核结构):
//   令 K(t,s) = (t-s)^{H-1/2} / Gamma(H+1/2)  (归一化核)
//   则 W̃^H_t = ∫_0^t K(t,s) dW_s,  Cov = ∫_0^{min(t_i,t_j)} K(t_i,s) K(t_j,s) ds
//
// 数值实现: 对均匀网格 dt, 协方差:
//   C[i][j] = dt * Σ_{k=0}^{min(i,j)-1} (t_i - t_k)^{H-1/2} (t_j - t_k)^{H-1/2} * dt
//           ≈ dt * ∫_0^{min(t_i,t_j)} (t_i-s)^{H-1/2} (t_j-s)^{H-1/2} ds
//   (用左 Riemann 和近似积分, 注意: 这是近似, 精确需 Gauss 积分)
//
// 为保证精度, 本实现用解析公式:
//   C[i][j] = (2H+1) * integral_{0}^{min(t_i,t_j)} (t_i-s)^{H-1/2} (t_j-s)^{H-1/2} ds
//   积分用代换 u = s/m, m=min(t_i,t_j):
//     = (2H+1) * m^{2H+1} * B(1, H+1/2, ...)  — 用超几何函数
//
// 简化实现: 直接用 Riemann 和 (对研究用途足够, 步长 dt 足够小时误差 O(dt)):
inline std::vector<std::vector<Real>> rbergomi_fbm_covariance(
    Real T, Size n_steps, Real H) {
    // 返回 n_steps x n_steps 协方差矩阵 (索引 0..n_steps-1 对应 t_1..t_n)
    // W̃^H_{t_i} = sqrt(2H+1) * ∫_0^{t_i} (t_i-s)^{H-1/2} dW_s
    // Cov(W̃^H_i, W̃^H_j) = (2H+1) * ∫_0^{min(t_i,t_j)} (t_i-s)^{H-1/2} (t_j-s)^{H-1/2} ds

    Real dt = T / static_cast<Real>(n_steps);
    Real alpha = H - 0.5;  // 核指数 (负值, 因 H < 0.5)
    Real coef = 2.0 * H + 1.0;

    std::vector<std::vector<Real>> C(n_steps, std::vector<Real>(n_steps, 0.0));

    for (Size i = 0; i < n_steps; ++i) {
        Real t_i = static_cast<Real>(i + 1) * dt;  // t_1 = dt, ..., t_n = T
        for (Size j = 0; j <= i; ++j) {
            Real t_j = static_cast<Real>(j + 1) * dt;
            Real m = std::min(t_i, t_j);  // = t_j (因 j <= i)
            // 积分 ∫_0^m (t_i - s)^{alpha} (t_j - s)^{alpha} ds
            // 用 Riemann 和 (左端点, k*dt 到 (k+1)*dt):
            Real integral = 0.0;
            Size n_int = std::min(i, j) + 1;  // 积分区间 [0, m] 分成 n_int 段
            for (Size k = 0; k < n_int; ++k) {
                Real s = static_cast<Real>(k) * dt;  // 左端点
                Real d_i = t_i - s;
                Real d_j = t_j - s;
                if (d_i > 0.0 && d_j > 0.0) {
                    integral += std::pow(d_i, alpha) * std::pow(d_j, alpha) * dt;
                }
            }
            C[i][j] = coef * integral;
            C[j][i] = C[i][j];  // 对称
        }
    }
    return C;
}

// ============ Cholesky 分解 (动态尺寸, std::vector 版本) ============
inline std::vector<std::vector<Real>> cholesky_dynamic(
    const std::vector<std::vector<Real>>& A) {
    Size n = A.size();
    if (n == 0) return {};
    for (auto& row : A) {
        if (row.size() != n) throw std::invalid_argument("cholesky: matrix must be square");
    }
    std::vector<std::vector<Real>> L(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) {
        for (Size j = 0; j <= i; ++j) {
            Real sum = A[i][j];
            for (Size k = 0; k < j; ++k) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 0.0) {
                    // 数值不稳定 (矩阵非正定, 通常因 Riemann 和近似误差)
                    // 用对角线微小正数修复
                    sum = std::max(sum, 1e-14);
                }
                L[i][j] = std::sqrt(sum);
            } else {
                if (std::abs(L[j][j]) < 1e-15) {
                    L[i][j] = 0.0;
                } else {
                    L[i][j] = sum / L[j][j];
                }
            }
        }
    }
    return L;
}

// ============ RL-fBm 采样器 (Cholesky 精确法) ============
// 给定均匀网格 t_1=dt, ..., t_n=T, 生成 W̃^H_{t_1}, ..., W̃^H_{t_n}
// 复杂度: O(n³) Cholesky (一次性) + O(n²) 每路径
class RLFbmSampler {
public:
    RLFbmSampler(Real T, Size n_steps, Real H)
        : T_(T), n_steps_(n_steps), H_(H) {
        if (T <= 0.0) throw std::invalid_argument("RLFbmSampler: T must be positive");
        if (n_steps == 0) throw std::invalid_argument("RLFbmSampler: n_steps must be positive");
        if (H <= 0.0 || H >= 0.5)
            throw std::invalid_argument("RLFbmSampler: H must be in (0, 0.5)");

        // 构造协方差矩阵并 Cholesky 分解 (一次性, 后续复用)
        auto C = rbergomi_fbm_covariance(T, n_steps, H);
        L_ = cholesky_dynamic(C);
    }

    // 采样一条 W̃^H 路径 (返回 n_steps 个值, 对应 t_1..t_n)
    // Z 为 n_steps 个独立 N(0,1) (已由调用方生成)
    std::vector<Real> sample(const std::vector<Real>& Z) const {
        if (Z.size() != n_steps_)
            throw std::invalid_argument("RLFbmSampler: Z size mismatch");
        std::vector<Real> W(n_steps_, 0.0);
        for (Size i = 0; i < n_steps_; ++i) {
            Real sum = 0.0;
            for (Size j = 0; j <= i; ++j) {
                sum += L_[i][j] * Z[j];
            }
            W[i] = sum;
        }
        return W;
    }

    Real T() const { return T_; }
    Size n_steps() const { return n_steps_; }
    Real H() const { return H_; }

    // 方差过程 v_t 的理论矩 (用于测试验证):
    // E[v_t] = ξ_0(t) (鞅条件, 因 exp(-0.5 η² t^{2H}) 抵消 W̃^H 的均值贡献)
    // Var(log v_t) = η² * t^{2H}
    // E[log v_t] = log(ξ_0(t)) - 0.5 η² t^{2H}
    static Real log_v_mean(Real t, Real xi0, Real eta, Real H) {
        return std::log(xi0) - 0.5 * eta * eta * std::pow(t, 2.0 * H);
    }
    static Real log_v_var(Real t, Real eta, Real H) {
        return eta * eta * std::pow(t, 2.0 * H);
    }

private:
    Real T_, H_;
    Size n_steps_;
    std::vector<std::vector<Real>> L_;  // Cholesky 下三角
};

// ============ Rough Bergomi 过程类 ============
class RoughBergomiProcess : public StochasticProcess {
public:
    explicit RoughBergomiProcess(RoughBergomiParams p)
        : params_(p) {
        validate_rough_bergomi_params(p);
    }

    Size dimension() const override { return 2; }  // (S, v) 或等价 (W, W̃^H)
    Real spot() const override { return params_.S0; }

    // rBergomi 无闭式特征函数 (与 rough Heston 不同), 返回 0 表示不可用
    Complex characteristic_function(Complex /*u*/, Real /*tau*/) const override {
        return Complex(0.0, 0.0);
    }

    // 生成一条路径 (价格过程 S_t)
    // path[0] = S0, path[i+1] = S_{t_{i+1}}
    void generate_path(Real T, Size n_steps,
                       std::span<Real> path, Philox4x64& rng) const override {
        if (path.size() < n_steps + 1)
            throw std::invalid_argument("rBergomi: path buffer too small");
        if (n_steps == 0) {
            path[0] = params_.S0;
            return;
        }

        Real dt = T / static_cast<Real>(n_steps);
        Real sqrt_dt = std::sqrt(dt);
        Real H = params_.H;
        Real eta = params_.eta;
        Real rho = params_.rho;
        Real xi0 = params_.xi0;
        Real sqrt_1m_rho2 = std::sqrt(1.0 - rho * rho);

        // 1. 生成方差驱动 Brownian 增量 dW (用于 W̃^H 和 v_t)
        std::vector<Real> dW(n_steps);
        // 同时生成正交 Brownian 增量 dW_perp (用于价格过程的独立部分)
        std::vector<Real> dW_perp(n_steps);
        for (Size i = 0; i < n_steps; ++i) {
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            auto [z1, z2] = box_muller(u1, u2);
            dW[i] = z1 * sqrt_dt;
            dW_perp[i] = z2 * sqrt_dt;
        }

        // 2. 构造 W̃^H 路径 (用 RLFbmSampler, 每次构造会重算 Cholesky — 低效)
        // 优化: 缓存 sampler (见 generate_path_with_sampler)
        RLFbmSampler sampler(T, n_steps, H);
        // W̃^H 采样需要 N(0,1) 序列, 用 dW / sqrt_dt 还原
        std::vector<Real> Z_for_fbm(n_steps);
        for (Size i = 0; i < n_steps; ++i) {
            Z_for_fbm[i] = dW[i] / sqrt_dt;  // 还原 N(0,1)
        }
        std::vector<Real> WH = sampler.sample(Z_for_fbm);
        // WH[i] = W̃^H_{t_{i+1}} (t_1=dt, ..., t_n=T)

        // 3. 构造价格路径
        // v_{t_i} = ξ_0 * exp(η * W̃^H_{t_i} - 0.5 * η² * t_i^{2H})
        // dS/S = sqrt(v_t) * dZ_t,  dZ_t = ρ*dW_t + sqrt(1-ρ²)*dW_perp_t
        path[0] = params_.S0;
        Real S = params_.S0;
        Real log_S = std::log(S);

        for (Size i = 0; i < n_steps; ++i) {
            Real t_i = static_cast<Real>(i + 1) * dt;  // 当前时间 t_{i+1}
            // v_{t_{i+1}} (用 WH[i], 对应 t_{i+1})
            Real t_to_2H = std::pow(t_i, 2.0 * H);
            Real log_v = std::log(xi0) + eta * WH[i] - 0.5 * eta * eta * t_to_2H;
            Real v = std::exp(log_v);

            // 价格增量: dS/S = sqrt(v_t) * dZ_t
            // 这里用 v_{t_{i+1}} 近似 (Euler), 更精确可用 v_{t_i} (上一步)
            // 但 v_0 = ξ_0 (t=0 时 W̃^H_0=0), 第一步用 ξ_0
            Real v_step;
            if (i == 0) {
                v_step = xi0;  // v_0
            } else {
                // 用上一步的 v (更稳定,避免前瞻)
                Real t_prev = static_cast<Real>(i) * dt;
                Real t_prev_2H = std::pow(t_prev, 2.0 * H);
                Real log_v_prev = std::log(xi0) + eta * WH[i - 1] - 0.5 * eta * eta * t_prev_2H;
                v_step = std::exp(log_v_prev);
            }

            Real dZ = rho * dW[i] + sqrt_1m_rho2 * dW_perp[i];
            // log-Euler: d(log S) = -0.5*v*dt + sqrt(v)*dZ
            log_S += -0.5 * v_step * dt + std::sqrt(v_step) * dZ;
            S = std::exp(log_S);
            path[i + 1] = S;
        }
    }

    // 优化版: 接收预构造的 sampler (避免每次重算 Cholesky)
    void generate_path_with_sampler(Real T, Size n_steps,
                                     std::span<Real> path, Philox4x64& rng,
                                     const RLFbmSampler& sampler) const {
        if (path.size() < n_steps + 1)
            throw std::invalid_argument("rBergomi: path buffer too small");
        if (n_steps == 0) {
            path[0] = params_.S0;
            return;
        }
        if (sampler.n_steps() != n_steps || std::abs(sampler.T() - T) > 1e-10)
            throw std::invalid_argument("rBergomi: sampler mismatch");

        Real dt = T / static_cast<Real>(n_steps);
        Real sqrt_dt = std::sqrt(dt);
        Real H = params_.H;
        Real eta = params_.eta;
        Real rho = params_.rho;
        Real xi0 = params_.xi0;
        Real sqrt_1m_rho2 = std::sqrt(1.0 - rho * rho);

        std::vector<Real> Z_for_fbm(n_steps);
        std::vector<Real> dW_perp(n_steps);
        for (Size i = 0; i < n_steps; ++i) {
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            auto [z1, z2] = box_muller(u1, u2);
            Z_for_fbm[i] = z1;       // 用于 fBm (N(0,1))
            dW_perp[i] = z2 * sqrt_dt;  // 正交 Brownian 增量
        }

        std::vector<Real> WH = sampler.sample(Z_for_fbm);

        path[0] = params_.S0;
        Real log_S = std::log(params_.S0);

        for (Size i = 0; i < n_steps; ++i) {
            Real t_i = static_cast<Real>(i + 1) * dt;
            Real t_to_2H = std::pow(t_i, 2.0 * H);
            Real log_v = std::log(xi0) + eta * WH[i] - 0.5 * eta * eta * t_to_2H;
            Real v = std::exp(log_v);
            (void)v;  // v 用于诊断, 实际定价用 v_step

            Real v_step;
            if (i == 0) {
                v_step = xi0;
            } else {
                Real t_prev = static_cast<Real>(i) * dt;
                Real t_prev_2H = std::pow(t_prev, 2.0 * H);
                Real log_v_prev = std::log(xi0) + eta * WH[i - 1] - 0.5 * eta * eta * t_prev_2H;
                v_step = std::exp(log_v_prev);
            }

            // dW_i = Z_for_fbm[i] * sqrt_dt (还原)
            Real dW_i = Z_for_fbm[i] * sqrt_dt;
            Real dZ = rho * dW_i + sqrt_1m_rho2 * dW_perp[i];
            log_S += -0.5 * v_step * dt + std::sqrt(v_step) * dZ;
            path[i + 1] = std::exp(log_S);
        }
    }

    // 生成方差路径 (用于测试和诊断)
    // v_path[i+1] = v_{t_{i+1}}
    void generate_variance_path(Real T, Size n_steps,
                                 std::span<Real> v_path, Philox4x64& rng,
                                 const RLFbmSampler& sampler) const {
        if (v_path.size() < n_steps + 1)
            throw std::invalid_argument("rBergomi: v_path buffer too small");
        if (n_steps == 0) {
            v_path[0] = params_.xi0;
            return;
        }

        Real dt = T / static_cast<Real>(n_steps);
        Real H = params_.H;
        Real eta = params_.eta;
        Real xi0 = params_.xi0;

        std::vector<Real> Z(n_steps);
        for (Size i = 0; i < n_steps; ++i) {
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            auto [z1, z2] = box_muller(u1, u2);
            Z[i] = z1;
            (void)z2;
        }

        std::vector<Real> WH = sampler.sample(Z);
        v_path[0] = xi0;  // v_0 = ξ_0
        for (Size i = 0; i < n_steps; ++i) {
            Real t_i = static_cast<Real>(i + 1) * dt;
            Real t_to_2H = std::pow(t_i, 2.0 * H);
            Real log_v = std::log(xi0) + eta * WH[i] - 0.5 * eta * eta * t_to_2H;
            v_path[i + 1] = std::exp(log_v);
        }
    }

    const RoughBergomiParams& params() const { return params_; }

private:
    RoughBergomiParams params_;
};

// ============ rBergomi MC 定价 (欧式期权) ============
// 因无闭式 CF, 用 MC 定价. 支持预构造 sampler 加速批量定价.
struct RoughBergomiMCResult {
    Real price;
    Real std_error;
    Real ci_lower;
    Real ci_upper;
    Size n_paths;
};

inline RoughBergomiMCResult rbergomi_price_european(
    const RoughBergomiParams& p, Real K, Real T, bool is_call,
    Size n_paths, uint64_t seed, Size n_steps = 100) {
    validate_rough_bergomi_params(p);
    if (T <= 0.0) throw std::invalid_argument("rBergomi MC: T must be positive");
    if (K <= 0.0) throw std::invalid_argument("rBergomi MC: K must be positive");
    if (n_paths == 0) throw std::invalid_argument("rBergomi MC: n_paths must be positive");

    RoughBergomiProcess process(p);
    RLFbmSampler sampler(T, n_steps, p.H);

    Real df = std::exp(-p.r * T);
    Real sum_payoff = 0.0, sum_payoff2 = 0.0;

    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(seed + j);
        std::vector<Real> path(n_steps + 1);
        process.generate_path_with_sampler(T, n_steps, path, rng, sampler);
        Real ST = path.back();
        Real payoff = is_call ? std::max(ST - K, 0.0) : std::max(K - ST, 0.0);
        sum_payoff += payoff;
        sum_payoff2 += payoff * payoff;
    }

    Real mean_payoff = sum_payoff / static_cast<Real>(n_paths);
    Real var_payoff = sum_payoff2 / static_cast<Real>(n_paths) - mean_payoff * mean_payoff;
    Real se = std::sqrt(std::max(var_payoff, 0.0) / static_cast<Real>(n_paths));
    Real price = df * mean_payoff;

    RoughBergomiMCResult result;
    result.price = price;
    result.std_error = df * se;
    result.ci_lower = price - 1.96 * result.std_error;
    result.ci_upper = price + 1.96 * result.std_error;
    result.n_paths = n_paths;
    return result;
}

}  // namespace v1
}  // namespace cpphub
