#pragma once
// SOURCE: Brown (1995) "Brownian bridge and the Peano curve" (原始构造思想)
// SOURCE: Caflisch, Morokoff, Owen (1997) "Valuation of mortgage-backed securities
//         using Brownian bridges to reduce effective dimension" J. Comput. Finance 1, 27-46.
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.5 (Variance Reduction)
// SOURCE: Acworth, Broadie, Glasserman (1998) "A comparison of some Monte Carlo and quasi
//         Monte Carlo techniques for option pricing" J. Comput. Finance 1, 105-115.
// 模块: Brownian Bridge 构造 (与 Sobol 低偏差序列配合, 降低有效维度)
//
// ==================== Brownian Bridge 数学 ====================
//
// 标准 Brownian motion W(t) 的样本路径构造:
//   方案 A (incremental): W(t_{i+1}) = W(t_i) + sqrt(t_{i+1}-t_i) * Z_i
//   方案 B (Brownian Bridge): 先生成 W(T), 再递归填充中点
//
// Brownian Bridge 构造 (等间隔 0 = t_0 < t_1 < ... < t_n = T):
//   1. 生成 W(T) = sqrt(T) * Z_0
//   2. 给定 W(t_left) 和 W(t_right), 中点:
//      W(t_mid) = 0.5*(W(t_left) + W(t_right)) + sqrt((t_mid-t_left)*(t_right-t_mid)/2) * Z
//   3. 递归直到所有节点被填充
//
// 任意时间网格 (非等间隔):
//   给定 0 = t_0 < t_1 < ... < t_n, 构造 W(t_0), W(t_1), ..., W(t_n)
//   - 第一步: W(t_n) ~ N(0, t_n), 即 W(t_n) = sqrt(t_n) * Z_0
//   - 后续: 选择下一个待生成节点 t_k (已生成 t_i, t_j, i < k < j)
//     条件分布: W(t_k) | W(t_i), W(t_j) ~ N(μ, σ²)
//       μ = ((t_j - t_k) * W(t_i) + (t_k - t_i) * W(t_j)) / (t_j - t_i)
//       σ² = (t_k - t_i) * (t_j - t_k) / (t_j - t_i)
//
// 优势: 最重要的维度 (Z_0) 控制 W(T), 这对欧式期权定价至关重要
//       与 Sobol 配合可将有效维度从 n 降到 ~log2(n), 显著提升收敛速度
//
// ==================== PCA 路径构造 (多资产相关 Brownian) ====================
//
// 多资产场景: W_i(t) (i=1..d) 为 d 个相关 Brownian motion, Corr(W_i, W_j) = ρ_ij
// 路径上 n+1 个时间点 × d 个资产 = (n+1)*d 维联合正态
// PCA 构造: 对协方差矩阵 Σ 做 eigenvalue 分解 Σ = V Λ V^T
//   - 按特征值降序排列: λ_1 ≥ λ_2 ≥ ... ≥ λ_{nd}
//   - 路径 X = V * sqrt(Λ) * Z (Z 为独立标准正态)
//   - 前 k 维 Sobol 控制 (1 - Σ_{i>k} λ_i / Σ λ_i) 的方差
// PCA 比 Brownian Bridge 更优 (理论最优), 但计算成本更高 (O(n^3 d^3))

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/pricing/monte_carlo/sobol.hpp"
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <functional>

namespace cpphub {
inline namespace v1 {

// ============ Brownian Bridge (单资产, 任意时间网格) ============
// 构造 W(t_0)=0, W(t_1), ..., W(t_n) 的样本路径
// 与 Sobol 序列配合: uniforms[0] → W(t_n), uniforms[1] → 中点, ...
class BrownianBridge {
public:
    // 等间隔网格: t_i = i * T / n_steps, i = 0, 1, ..., n_steps
    explicit BrownianBridge(Size n_steps, Real T = 1.0)
        : n_steps_(n_steps), T_(T) {
        if (n_steps == 0) throw std::invalid_argument("BrownianBridge: n_steps must be positive");
        if (T <= 0.0) throw std::invalid_argument("BrownianBridge: T must be positive");
        times_.resize(n_steps_ + 1);
        for (Size i = 0; i <= n_steps_; ++i) {
            times_[i] = T_ * static_cast<Real>(i) / static_cast<Real>(n_steps_);
        }
        build_schedule();
    }

    // 任意时间网格: times[0] = 0 < times[1] < ... < times[n]
    explicit BrownianBridge(const std::vector<Real>& times)
        : n_steps_(times.size() - 1), T_(times.back()), times_(times) {
        if (times.size() < 2) throw std::invalid_argument("BrownianBridge: need >= 2 time points");
        if (times[0] != 0.0) throw std::invalid_argument("BrownianBridge: times[0] must be 0");
        for (Size i = 1; i < times.size(); ++i) {
            if (times[i] <= times[i - 1]) {
                throw std::invalid_argument("BrownianBridge: times must be strictly increasing");
            }
        }
        build_schedule();
    }

    // 从 uniforms (dim = n_steps) 生成 Brownian 增量 dW_i = W(t_i) - W(t_{i-1}), i=1..n
    // 返回长度 n_steps, increments[0] = W(t_1) - W(t_0) = W(t_1)
    std::vector<Real> generate(const std::vector<Real>& uniforms) const {
        auto W = generate_path(uniforms);  // W[0..n-1] = W(t_1)..W(t_n), 长度 n_steps
        std::vector<Real> increments(n_steps_);
        // dW_1 = W(t_1) - W(t_0) = W(t_1) - 0 = W[0]
        increments[0] = W[0];
        for (Size i = 1; i < n_steps_; ++i) {
            increments[i] = W[i] - W[i - 1];
        }
        return increments;
    }

    // 从 uniforms (dim = n_steps) 生成 Brownian 路径 W(t_1), W(t_2), ..., W(t_n)
    // 注意: 不包含 W(t_0)=0, 返回长度 n_steps
    std::vector<Real> generate_path(const std::vector<Real>& uniforms) const {
        if (uniforms.size() < schedule_.size()) {
            throw std::invalid_argument("BrownianBridge: uniforms size < schedule size");
        }
        std::vector<Real> W(n_steps_ + 1, 0.0);  // W[0] = 0
        for (Size s = 0; s < schedule_.size(); ++s) {
            Real u = uniforms[s];
            // Sobol 可能产生 0 或 1, 需要 clip
            if (u <= 0.0) u = 1e-15;
            if (u >= 1.0) u = 1.0 - 1e-15;
            Real Z = inv_normal_cdf(u);
            const auto& step = schedule_[s];
            Real mean = (static_cast<Real>(step.right - step.index) * W[step.left]
                       + static_cast<Real>(step.index - step.left) * W[step.right])
                      / static_cast<Real>(step.right - step.left);
            W[step.index] = mean + step.sigma * Z;
        }
        // 返回 W(t_1)..W(t_n), 去掉 W(t_0)=0
        return std::vector<Real>(W.begin() + 1, W.end());
    }

    // 生成完整路径 (包含 W(t_0)=0), 长度 n_steps+1
    std::vector<Real> generate_full_path(const std::vector<Real>& uniforms) const {
        auto path = generate_path(uniforms);
        path.insert(path.begin(), 0.0);
        return path;
    }

    // 所需 Sobol 维度 (= n_steps_)
    Size dimension() const noexcept { return n_steps_; }
    Size n_steps() const noexcept { return n_steps_; }
    Real T() const noexcept { return T_; }
    const std::vector<Real>& times() const noexcept { return times_; }

private:
    Size n_steps_;
    Real T_;
    std::vector<Real> times_;  // 时间网格 (长度 n_steps+1, times_[0]=0)

    struct BridgeStep {
        Size index;   // 待生成的节点索引 (0..n_steps)
        Size left;    // 已生成的左节点
        Size right;   // 已生成的右节点
        Real sigma;   // 条件标准差
    };

    std::vector<BridgeStep> schedule_;

    void build_schedule() {
        schedule_.clear();
        // 第一步: 生成 W(t_n) ~ N(0, T)
        schedule_.push_back({n_steps_, 0, n_steps_, std::sqrt(times_[n_steps_])});

        // 递归填充中点 (bisection)
        std::function<void(Size, Size)> bridge = [&](Size left, Size right) {
            if (right - left <= 1) return;
            Size mid = (left + right) / 2;
            Real t_left = times_[left];
            Real t_right = times_[right];
            Real t_mid = times_[mid];
            // σ² = (t_mid - t_left) * (t_right - t_mid) / (t_right - t_left)
            Real sigma_sq = (t_mid - t_left) * (t_right - t_mid) / (t_right - t_left);
            Real sigma = std::sqrt(sigma_sq);
            schedule_.push_back({mid, left, right, sigma});
            bridge(left, mid);
            bridge(mid, right);
        };
        bridge(0, n_steps_);
    }
};

// ============ 多资产 Brownian Bridge (相关 Brownian) ============
// 给定 d 个资产的相关矩阵 R = LL^T, 生成 d 条相关 Brownian 路径
// 策略: 每个时间步用 d 维 Sobol 生成 d 个独立 Brownian Bridge, 再用 L 相关化
// 维度: n_steps * n_assets (按时间分组, 每组 n_assets 维)
class MultiAssetBrownianBridge {
public:
    // 等间隔时间网格 + 相关矩阵
    MultiAssetBrownianBridge(Size n_steps, Real T, Size n_assets,
                               const std::vector<std::vector<Real>>& correlation)
        : n_steps_(n_steps), T_(T), n_assets_(n_assets) {
        if (n_steps == 0) throw std::invalid_argument("MultiAssetBB: n_steps must be positive");
        if (T <= 0.0) throw std::invalid_argument("MultiAssetBB: T must be positive");
        if (n_assets == 0) throw std::invalid_argument("MultiAssetBB: n_assets must be positive");
        bb_ = std::make_shared<BrownianBridge>(n_steps, T);
        L_ = cholesky_semi_definite(correlation);
        if (L_.size() != n_assets) {
            throw std::invalid_argument("MultiAssetBB: correlation size mismatch");
        }
    }

    // 从 uniforms (dim = n_steps * n_assets) 生成 d 条相关 Brownian 路径
    // 返回 paths[asset][step], step 0..n_steps-1 (不含 W(0)=0)
    // uniforms 布局: uniforms[step * n_assets + asset] 用于第 step 步第 asset 个 Brownian
    std::vector<std::vector<Real>> generate_paths(const std::vector<Real>& uniforms) const {
        const Size expected_dim = n_steps_ * n_assets_;
        if (uniforms.size() < expected_dim) {
            throw std::invalid_argument("MultiAssetBB: uniforms size < n_steps * n_assets");
        }

        // 步骤 1: 对每个资产, 用 Brownian Bridge 生成独立路径
        // 需要从 uniforms 中按 (step, asset) 布局提取每资产的 n_steps 维子序列
        // 但 Brownian Bridge 期望按 BB 顺序 (W(T) 优先, 然后中点), 所以需要重排
        // 简化: 对每个 asset, 提取 uniforms[asset::n_assets] 作为该 asset 的 Sobol 序列
        // (即第 asset 个 asset 的第 k 个 Sobol 点 = uniforms[k * n_assets + asset])
        std::vector<std::vector<Real>> independent_paths(n_assets_);
        for (Size a = 0; a < n_assets_; ++a) {
            std::vector<Real> asset_uniforms(n_steps_);
            for (Size s = 0; s < n_steps_; ++s) {
                asset_uniforms[s] = uniforms[s * n_assets_ + a];
            }
            independent_paths[a] = bb_->generate_path(asset_uniforms);
        }

        // 步骤 2: 每个时间步应用 Cholesky 相关化
        // W_correlated[a][s] = Σ_b L[a][b] * W_independent[b][s]
        std::vector<std::vector<Real>> correlated_paths(n_assets_,
            std::vector<Real>(n_steps_, 0.0));
        for (Size s = 0; s < n_steps_; ++s) {
            for (Size a = 0; a < n_assets_; ++a) {
                Real sum = 0.0;
                for (Size b = 0; b <= a; ++b) {
                    sum += L_[a][b] * independent_paths[b][s];
                }
                correlated_paths[a][s] = sum;
            }
        }
        return correlated_paths;
    }

    // 生成 Brownian 增量 (而非路径), 返回 paths[asset][step], step 0..n_steps-1
    std::vector<std::vector<Real>> generate_increments(const std::vector<Real>& uniforms) const {
        auto paths = generate_paths(uniforms);
        // W(t_0)=0, W(t_1), ..., W(t_n) → 增量 dW_i = W(t_i) - W(t_{i-1})
        // paths[asset] 长度 n_steps (W(t_1)..W(t_n)), 增量 = paths[a][s] - (s==0 ? 0 : paths[a][s-1])
        std::vector<std::vector<Real>> incs(n_assets_, std::vector<Real>(n_steps_, 0.0));
        for (Size a = 0; a < n_assets_; ++a) {
            for (Size s = 0; s < n_steps_; ++s) {
                Real prev = (s == 0) ? 0.0 : paths[a][s - 1];
                incs[a][s] = paths[a][s] - prev;
            }
        }
        return incs;
    }

    Size dimension() const noexcept { return n_steps_ * n_assets_; }
    Size n_steps() const noexcept { return n_steps_; }
    Size n_assets() const noexcept { return n_assets_; }
    Real T() const noexcept { return T_; }
    const std::vector<std::vector<Real>>& cholesky_L() const noexcept { return L_; }

private:
    Size n_steps_;
    Real T_;
    Size n_assets_;
    std::shared_ptr<BrownianBridge> bb_;
    std::vector<std::vector<Real>> L_;  // Cholesky 下三角

    // 半正定 Cholesky (允许秩亏)
    static std::vector<std::vector<Real>> cholesky_semi_definite(
            const std::vector<std::vector<Real>>& A) {
        const Size n = A.size();
        if (n == 0) throw std::invalid_argument("Cholesky: empty matrix");
        for (Size i = 0; i < n; ++i) {
            if (A[i].size() != n) throw std::invalid_argument("Cholesky: non-square matrix");
        }
        std::vector<std::vector<Real>> L(n, std::vector<Real>(n, 0.0));
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j <= i; ++j) {
                Real s = A[i][j];
                for (Size k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
                if (i == j) {
                    if (s < -1e-12)
                        throw std::invalid_argument("Cholesky: matrix not positive semidefinite");
                    L[i][j] = std::sqrt(std::max(s, 0.0));
                } else {
                    L[i][j] = (L[j][j] > 1e-15) ? s / L[j][j] : 0.0;
                }
            }
        }
        return L;
    }
};

// ============ PCA 路径构造 (理论最优, 多资产) ============
// 对协方差矩阵 Σ 做 eigenvalue 分解, 按特征值降序排列
// 路径 X = V * sqrt(Λ) * Z, 前 k 维 Sobol 控制最大方差
// 适用于非等间隔时间网格和复杂相关结构
// 注: 完整 PCA 需要 O((n*d)^3), 对大 n*d 不实用, 仅提供接口
//      实际使用中 Brownian Bridge 已足够, PCA 仅作高级选项
class BrownianBridgePCA {
public:
    // 等间隔网格 + 相关矩阵 (n_assets 个资产, n_steps 步)
    BrownianBridgePCA(Size n_steps, Real T, Size n_assets,
                       const std::vector<std::vector<Real>>& correlation)
        : n_steps_(n_steps), T_(T), n_assets_(n_assets) {
        if (n_steps == 0 || n_assets == 0 || T <= 0.0) {
            throw std::invalid_argument("BrownianBridgePCA: invalid parameters");
        }
        build_covariance_matrix(correlation);
        decompose_covariance();
    }

    // 从 uniforms (dim = n_steps * n_assets) 生成 d 条相关 Brownian 路径
    // uniforms 按特征值降序对应: uniforms[0] → 最大特征值方向
    std::vector<std::vector<Real>> generate_paths(const std::vector<Real>& uniforms) const {
        const Size dim = n_steps_ * n_assets_;
        if (uniforms.size() < dim) {
            throw std::invalid_argument("BrownianBridgePCA: uniforms size < dimension");
        }
        // Z = inv_normal_cdf(uniforms), 按特征值降序
        std::vector<Real> Z(dim);
        for (Size i = 0; i < dim; ++i) {
            Real u = uniforms[i];
            if (u <= 0.0) u = 1e-15;
            if (u >= 1.0) u = 1.0 - 1e-15;
            Z[i] = inv_normal_cdf(u);
        }
        // X = V * sqrt(Λ) * Z = (V * sqrt(Λ)) * Z
        // eigenvectors_[i][j] = V[j][i] (第 i 个特征向量)
        // X[j] = Σ_i V[j][i] * sqrt(λ_i) * Z[i]
        std::vector<Real> X(dim, 0.0);
        for (Size i = 0; i < dim; ++i) {
            Real sqrt_lambda = std::sqrt(std::max(eigenvalues_[i], 0.0));
            Real coeff = sqrt_lambda * Z[i];
            for (Size j = 0; j < dim; ++j) {
                X[j] += eigenvectors_[i][j] * coeff;
            }
        }
        // 解布局 X 为 [asset][step] (step 0..n_steps-1, 不含 W(0)=0)
        // 但 Brownian 路径需要 W(0)=0, 而 PCA 直接生成 W(t_1)..W(t_n)
        std::vector<std::vector<Real>> paths(n_assets_,
            std::vector<Real>(n_steps_, 0.0));
        for (Size a = 0; a < n_assets_; ++a) {
            for (Size s = 0; s < n_steps_; ++s) {
                paths[a][s] = X[a * n_steps_ + s];
            }
        }
        return paths;
    }

    Size dimension() const noexcept { return n_steps_ * n_assets_; }
    const std::vector<Real>& eigenvalues() const noexcept { return eigenvalues_; }
    // 累计方差解释比例 (前 k 维)
    Real cumulative_variance_ratio(Size k) const {
        if (k == 0 || k > eigenvalues_.size()) return 0.0;
        Real total = 0.0, partial = 0.0;
        for (Size i = 0; i < eigenvalues_.size(); ++i) {
            total += eigenvalues_[i];
            if (i < k) partial += eigenvalues_[i];
        }
        return (total > 0.0) ? partial / total : 0.0;
    }

private:
    Size n_steps_;
    Real T_;
    Size n_assets_;
    std::vector<Real> eigenvalues_;                    // 降序排列
    std::vector<std::vector<Real>> eigenvectors_;      // eigenvectors_[i] = 第 i 个特征向量

    // 构造协方差矩阵 Σ (size = n_steps*n_assets)
    // Σ[(a,s)][(b,t)] = Cov(W_a(t_s), W_b(t_t)) = ρ_ab * min(t_s, t_t)
    void build_covariance_matrix(const std::vector<std::vector<Real>>& correlation) {
        const Size dim = n_steps_ * n_assets_;
        if (correlation.size() != n_assets_) {
            throw std::invalid_argument("PCA: correlation size mismatch");
        }
        for (Size a = 0; a < n_assets_; ++a) {
            if (correlation[a].size() != n_assets_) {
                throw std::invalid_argument("PCA: correlation non-square");
            }
        }
        cov_.assign(dim, std::vector<Real>(dim, 0.0));
        Real dt = T_ / static_cast<Real>(n_steps_);
        for (Size a = 0; a < n_assets_; ++a) {
            for (Size s = 0; s < n_steps_; ++s) {
                Real t_s = (s + 1) * dt;  // W(t_1)..W(t_n), t_s = (s+1)*dt
                Size idx_as = a * n_steps_ + s;
                for (Size b = 0; b < n_assets_; ++b) {
                    for (Size t = 0; t < n_steps_; ++t) {
                        Real t_t = (t + 1) * dt;
                        Size idx_bt = b * n_steps_ + t;
                        cov_[idx_as][idx_bt] = correlation[a][b] * std::min(t_s, t_t);
                    }
                }
            }
        }
    }

    // Jacobi eigenvalue 分解 (对称矩阵)
    // 返回 eigenvalues (降序) 和 eigenvectors (按特征值降序排列)
    void decompose_covariance() {
        const Size n = cov_.size();
        eigenvalues_.resize(n);
        std::vector<std::vector<Real>> V(n, std::vector<Real>(n, 0.0));
        for (Size i = 0; i < n; ++i) V[i][i] = 1.0;

        std::vector<std::vector<Real>> A = cov_;  // 工作副本
        const Size max_sweeps = 100;
        const Real tol = 1e-12;

        for (Size sweep = 0; sweep < max_sweeps; ++sweep) {
            // 计算非对角元平方和
            Real off = 0.0;
            for (Size i = 0; i < n; ++i) {
                for (Size j = i + 1; j < n; ++j) {
                    off += A[i][j] * A[i][j];
                }
            }
            if (off < tol) break;

            for (Size p = 0; p < n; ++p) {
                for (Size q = p + 1; q < n; ++q) {
                    if (std::abs(A[p][q]) < 1e-20) continue;
                    Real app = A[p][p], aqq = A[q][q], apq = A[p][q];
                    Real phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
                    Real c = std::cos(phi), s = std::sin(phi);
                    // 旋转
                    for (Size i = 0; i < n; ++i) {
                        Real aip = A[i][p], aiq = A[i][q];
                        A[i][p] = c * aip - s * aiq;
                        A[i][q] = s * aip + c * aiq;
                    }
                    for (Size j = 0; j < n; ++j) {
                        Real apj = A[p][j], aqj = A[q][j];
                        A[p][j] = c * apj - s * aqj;
                        A[q][j] = s * apj + c * aqj;
                    }
                    for (Size i = 0; i < n; ++i) {
                        Real vip = V[i][p], viq = V[i][q];
                        V[i][p] = c * vip - s * viq;
                        V[i][q] = s * vip + c * viq;
                    }
                }
            }
        }
        for (Size i = 0; i < n; ++i) eigenvalues_[i] = A[i][i];

        // 按特征值降序排列
        std::vector<Size> order(n);
        for (Size i = 0; i < n; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
            [&](Size a, Size b) { return eigenvalues_[a] > eigenvalues_[b]; });

        std::vector<Real> sorted_eigenvalues(n);
        std::vector<std::vector<Real>> sorted_eigenvectors(n, std::vector<Real>(n));
        for (Size i = 0; i < n; ++i) {
            sorted_eigenvalues[i] = eigenvalues_[order[i]];
            for (Size j = 0; j < n; ++j) {
                sorted_eigenvectors[i][j] = V[j][order[i]];
            }
        }
        eigenvalues_ = std::move(sorted_eigenvalues);
        eigenvectors_ = std::move(sorted_eigenvectors);
    }

    std::vector<std::vector<Real>> cov_;  // 协方差矩阵 (工作用)
};

}  // namespace v1
}  // namespace cpphub
