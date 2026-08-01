// RBF-FD (Radial Basis Function - Finite Difference) method for multi-asset PDEs
// Reference:
//   Wright & Fornberg (2016) "Scattered node compact finite difference-type
//     formulas generated from radial basis functions"
//   Fornberg & Flyer (2015) "A Primer on Radial Basis Functions with Applications
//     to the Geosciences"
//   Wright, Fornberg, Flyer (2018) "A stable algorithm for flat RBFs"
//
// 多资产 BSM PDE (log-S 变换, x_i = ln(S_i/K_i)):
//   ∂V/∂τ = 0.5 * Σ_ij a_ij ∂²V/∂x_i∂x_j + Σ_i b_i ∂V/∂x_i - r V
//   其中 a_ij = ρ_ij σ_i σ_j, b_i = r - q_i - 0.5 σ_i², τ = T - t
//
// RBF-FD 离散:
//   1. 在 d 维空间生成节点集 {x_1, ..., x_N}
//   2. 对每个节点 x_i, 找到 k 个最近邻 (stencil)
//   3. 用 RBF + 多项式增强计算局部微分权重
//   4. 组装稀疏全局微分矩阵 L
//   5. θ-scheme 时间步进: (I - θ dt L) V^{n+1} = (I + (1-θ) dt L) V^n
//
// Polyharmonic Spline (PHS) RBF φ(r) = r^(2m-1) (奇数次):
//   - 无形状参数, 避免 ε 选择病态
//   - 配合多项式增强 (degree p), 收敛阶 = O(h^min(p+1, 2m))
//   - 推荐配置: PHS m=3 (r^5) + poly_degree=2
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <map>
#include <limits>
#include <iostream>
#include <cstdlib>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"

namespace cpphub {
inline namespace v1 {

// ============ RBF 核函数类型 ============
enum class RBFType {
    Polyharmonic3,    // φ(r) = r³        (推荐 2D/3D)
    Polyharmonic5,    // φ(r) = r⁵        (推荐高精度)
    Polyharmonic7,    // φ(r) = r⁷
    Gaussian,         // φ(r) = exp(-(εr)²)
    Multiquadric,     // φ(r) = sqrt(1 + (εr)²)
    InverseMultiquadric  // φ(r) = 1/sqrt(1 + (εr)²)
};

// ============ RBF-FD 配置 ============
struct RBFFDConfig {
    // 空间离散
    Size n_per_dim = 25;          // 每维节点数 (张量积, 总节点数 = n_per_dim^d)
    Size stencil_size = 20;       // 每个 stencil 的邻居数 (>= poly_terms + dim + 1)
    RBFType rbf = RBFType::Polyharmonic5;
    Real epsilon = 1.0;           // 形状参数 (仅 Gaussian/MQ/IMQ)
    Size poly_degree = 2;         // 附加多项式阶数 (0/1/2, 推荐 2)

    // 时间离散
    Size n_time = 100;            // 时间步数
    Real theta = 0.5;             // θ-scheme: 0.5=CN, 1.0=Implicit Euler

    // 计算域 (每维 log-S 空间的半宽, 单位 σ_i√T)
    Real domain_width = 5.0;

    // 求解器参数
    Size max_iter = 1000;         // 迭代求解器最大迭代数
    Real solver_tol = 1e-10;      // 迭代求解器容差
};

// ============ 多资产 BSM 模型参数 ============
struct MultiAssetBSMParams {
    std::vector<Real> S0;         // 初始价格 [d]
    std::vector<Real> sigma;      // 波动率 [d]
    std::vector<Real> q;          // 股息率 [d]
    std::vector<std::vector<Real>> corr;  // d×d 相关矩阵
    Real r = 0.05;                // 无风险利率
    Real T = 1.0;                 // 到期时间
    std::vector<Real> K;          // 行权价 [d] (用于 log-S 变换中心, 默认 = S0)
};

// ============ RBF-FD 引擎 ============
class RBFFDEngine {
public:
    explicit RBFFDEngine(RBFFDConfig config = RBFFDConfig{})
        : config_(config) {
        validate_config();
    }

    // 多资产欧式期权定价
    // payoff: 接收 d 维现货价格向量 S, 返回 payoff(S)
    // 返回: t=0 时刻 S=S0 处的期权价格
    Real price_european(
        const std::function<Real(const std::vector<Real>&)>& payoff,
        const MultiAssetBSMParams& params) const {

        Size d = params.S0.size();
        if (d == 0 || d > 3)
            throw std::invalid_argument("RBFFD: dimension must be 1, 2, or 3");
        if (params.sigma.size() != d || params.q.size() != d)
            throw std::invalid_argument("RBFFD: sigma/q size mismatch");
        if (params.corr.size() != d)
            throw std::invalid_argument("RBFFD: corr matrix size mismatch");
        for (Size i = 0; i < d; ++i) {
            if (params.corr[i].size() != d)
                throw std::invalid_argument("RBFFD: corr row size mismatch");
        }

        // 行权价 K (默认 = S0, 即 ATM 中心)
        std::vector<Real> K_vec = params.K.empty() ? params.S0 : params.K;
        if (K_vec.size() != d)
            throw std::invalid_argument("RBFFD: K size mismatch");

        // PDE 系数: a_ij = ρ_ij σ_i σ_j, b_i = r - q_i - 0.5 σ_i²
        std::vector<std::vector<Real>> a_matrix(d, std::vector<Real>(d, 0.0));
        std::vector<Real> b_vector(d);
        for (Size i = 0; i < d; ++i) {
            b_vector[i] = params.r - params.q[i] - 0.5 * params.sigma[i] * params.sigma[i];
            for (Size j = 0; j < d; ++j) {
                a_matrix[i][j] = params.corr[i][j] * params.sigma[i] * params.sigma[j];
            }
        }

        // 1. 生成节点集 (log-S 空间)
        auto nodes = generate_tensor_nodes(params.sigma, params.T, config_.domain_width);
        Size N = nodes.size();

        // 2. 组装稀疏微分矩阵 L (先按行收集, 再转 CSR)
        // 注意: compute_stencil_weights 返回的权重 w_m 已编码完整算子 L
        //   (包括 -rate 项, 因为 rbf_operator 和 poly_operator 都包含 -rate*phi)
        //   满足 L f(x_i) ≈ Σ_m w_m f(x_stencil_m)
        // 因此 L[i, stencil_m] = w_m, L[i, i] = 0 (i 不在自己的 stencil 中)
        // 行和 = Σ_m w_m = -rate (由多项式约束 Σ w_m = L[1] = -rate 保证)
        // 故无需再对对角线额外加 -rate
        std::vector<std::vector<std::pair<Size, Real>>> L_rows(N);

        for (Size i = 0; i < N; ++i) {
            auto stencil = find_knn(nodes, i, config_.stencil_size);
            auto weights = compute_stencil_weights(nodes, i, stencil,
                                                    a_matrix, b_vector, params.r, d);

            // stencil 包含 i 自身, 权重 w_m 满足 L f(x_i) ≈ Σ_m w_m f(x_{stencil_m})
            // 对角线 L[i,i] = w_{对应 i 的位置} (通常非零, 提供隐式格式的对角占优)
            std::map<Size, Real> row_map;
            for (Size m = 0; m < stencil.size(); ++m) {
                row_map[stencil[m]] += weights[m];
            }

            for (const auto& kv : row_map) {
                L_rows[i].push_back({kv.first, kv.second});
            }
            // 已按 col 排序 (std::map 保证)
        }

        // 转换为 CSR
        std::vector<Size> L_row_ptr(N + 1, 0);
        std::vector<Size> L_col_idx;
        std::vector<Real> L_values;
        for (Size i = 0; i < N; ++i) {
            L_row_ptr[i + 1] = L_row_ptr[i] + L_rows[i].size();
            for (const auto& e : L_rows[i]) {
                L_col_idx.push_back(e.first);
                L_values.push_back(e.second);
            }
        }

        // 3. 初始条件: V(x, τ=0) = payoff(S(x))
        std::vector<Real> V(N), V_new(N);
        std::vector<Real> payoff_at_nodes(N);
        for (Size n = 0; n < N; ++n) {
            std::vector<Real> S(d);
            for (Size i = 0; i < d; ++i) {
                S[i] = K_vec[i] * std::exp(nodes[n][i]);
            }
            payoff_at_nodes[n] = payoff(S);
            V[n] = payoff_at_nodes[n];
        }

        // 3b. 识别边界节点 (张量积网格上至少一维处于边界的节点)
        // 对这些节点强制 Dirichlet BC: V(τ) = payoff * exp(-rate*τ)
        // (折现内在价值近似, 远场精确 for OTM, 近似 for ITM)
        std::vector<bool> is_boundary(N, false);
        std::vector<Size> boundary_indices;
        identify_boundary_nodes(nodes, config_.n_per_dim, d, is_boundary, boundary_indices);

        // 4. 时间步进: (I - θ dt L) V^{n+1} = (I + (1-θ) dt L) V^n
        Real dt = params.T / static_cast<Real>(config_.n_time);
        Real theta = config_.theta;

        // 构建 A = I - θ dt L 的 CSR
        std::vector<Size> A_row_ptr = L_row_ptr;
        std::vector<Size> A_col_idx = L_col_idx;
        std::vector<Real> A_values = L_values;
        for (Size i = 0; i < A_values.size(); ++i) {
            A_values[i] *= -theta * dt;
        }
        // 加上单位矩阵 (对角线 +1)
        for (Size i = 0; i < N; ++i) {
            for (Size idx2 = A_row_ptr[i]; idx2 < A_row_ptr[i + 1]; ++idx2) {
                if (A_col_idx[idx2] == i) {
                    A_values[idx2] += 1.0;
                    break;
                }
            }
        }
        // 边界行替换为单位向量: A[boundary, :] = e_boundary
        // 这样 V_new[boundary] = rhs[boundary] = payoff * exp(-r*tau) (Dirichlet BC)
        for (Size bi : boundary_indices) {
            for (Size idx = A_row_ptr[bi]; idx < A_row_ptr[bi + 1]; ++idx) {
                A_values[idx] = (A_col_idx[idx] == bi) ? 1.0 : 0.0;
            }
        }

        std::vector<Real> rhs(N), LV(N);
        bool dbg = std::getenv("RBF_DEBUG") != nullptr;
        if (dbg) {
            Size center = (N - 1) / 2;
            std::cout << "[RBF-DBG] N=" << N << " center=" << center
                      << " dt=" << dt << " theta=" << theta
                      << " n_boundary=" << boundary_indices.size() << "\n";
            std::cout << "[RBF-DBG] L row " << center << " ("
                      << (L_row_ptr[center+1]-L_row_ptr[center]) << " nnz): ";
            for (Size idx = L_row_ptr[center]; idx < L_row_ptr[center+1]; ++idx)
                std::cout << L_col_idx[idx] << ":" << L_values[idx] << " ";
            std::cout << "\n";
            Real rowsum = 0;
            for (Size idx = L_row_ptr[center]; idx < L_row_ptr[center+1]; ++idx)
                rowsum += L_values[idx];
            std::cout << "[RBF-DBG] L row sum at center = " << rowsum
                      << " (expect ~" << -params.r << ")\n";
            std::cout << "[RBF-DBG] V[center]=" << V[center]
                      << " V[center+1]=" << V[center+1] << "\n";
        }
        for (Size step = 0; step < config_.n_time; ++step) {
            // rhs = V + (1-θ) dt L V
            spmv_csr(L_row_ptr, L_col_idx, L_values, V, LV);
            for (Size i = 0; i < N; ++i) {
                rhs[i] = V[i] + (1.0 - theta) * dt * LV[i];
            }

            // 强制 Dirichlet BC: rhs[boundary] = payoff * exp(-rate*tau)
            Real tau = static_cast<Real>(step + 1) * dt;
            Real bc_discount = std::exp(-params.r * tau);
            for (Size bi : boundary_indices) {
                rhs[bi] = payoff_at_nodes[bi] * bc_discount;
            }

            // 解 A V_new = rhs
            bool ok = bicgstab(A_row_ptr, A_col_idx, A_values, rhs, V_new);
            if (dbg && (step < 5 || step % 10 == 0 || step == config_.n_time - 1 || !ok)) {
                Size center = (N - 1) / 2;
                // count negative interior nodes
                Size n_neg = 0;
                for (Size i = 0; i < N; ++i)
                    if (!is_boundary[i] && V_new[i] < 0.0) ++n_neg;
                std::cout << "[RBF-DBG] step " << step << ": ok=" << ok
                          << " rhs[center]=" << rhs[center]
                          << " V_new[center]=" << V_new[center]
                          << " n_neg_interior=" << n_neg << "/" << (N - boundary_indices.size()) << "\n";
            }
            if (!ok) {
                if (dbg) std::cout << "[RBF-DBG] BiCGSTAB failed at step " << step << "\n";
                break;
            }

            // 强制边界节点值 (防止 BiCGSTAB 微小残差累积)
            for (Size bi : boundary_indices) {
                V_new[bi] = payoff_at_nodes[bi] * bc_discount;
            }

            // 确保非负 (期权价格下界, 仅对内部节点)
            for (Size i = 0; i < N; ++i) {
                if (!is_boundary[i] && V_new[i] < 0.0) V_new[i] = 0.0;
            }
            V = V_new;
        }
        if (dbg) {
            Size center = (N - 1) / 2;
            std::cout << "[RBF-DBG] final V[center]=" << V[center] << "\n";
        }

        // 5. 插值得到 S0 处的价格
        return interpolate_at_s0(nodes, V, params.S0, K_vec, params.sigma, params.T);
    }

    const RBFFDConfig& config() const { return config_; }

private:
    RBFFDConfig config_;

    void validate_config() const {
        if (config_.n_per_dim < 5)
            throw std::invalid_argument("RBFFD: n_per_dim must be >= 5");
        if (config_.stencil_size < 5)
            throw std::invalid_argument("RBFFD: stencil_size must be >= 5");
        if (config_.n_time == 0)
            throw std::invalid_argument("RBFFD: n_time must be > 0");
        if (config_.theta < 0.0 || config_.theta > 1.0)
            throw std::invalid_argument("RBFFD: theta must be in [0, 1]");
        if (config_.poly_degree > 2)
            throw std::invalid_argument("RBFFD: poly_degree must be 0, 1, or 2");
    }

    // ============ RBF 核函数 ============
    Real rbf_eval(Real r) const {
        switch (config_.rbf) {
            case RBFType::Polyharmonic3:
                return r * r * r;
            case RBFType::Polyharmonic5:
                return r * r * r * r * r;
            case RBFType::Polyharmonic7:
                return r * r * r * r * r * r * r;
            case RBFType::Gaussian: {
                Real e = config_.epsilon * r;
                return std::exp(-e * e);
            }
            case RBFType::Multiquadric: {
                Real e = config_.epsilon * r;
                return std::sqrt(1.0 + e * e);
            }
            case RBFType::InverseMultiquadric: {
                Real e = config_.epsilon * r;
                return 1.0 / std::sqrt(1.0 + e * e);
            }
        }
        return 0.0;
    }

    // PHS 的导数 (用于计算微分算子右端项)
    // d/dr φ(r) for r > 0
    Real rbf_deriv_r(Real r) const {
        if (r < 1e-15) r = 1e-15;
        switch (config_.rbf) {
            case RBFType::Polyharmonic3:
                return 3.0 * r * r;
            case RBFType::Polyharmonic5:
                return 5.0 * r * r * r * r;
            case RBFType::Polyharmonic7:
                return 7.0 * r * r * r * r * r * r;
            case RBFType::Gaussian: {
                Real e = config_.epsilon * r;
                return -2.0 * config_.epsilon * config_.epsilon * r * std::exp(-e * e);
            }
            case RBFType::Multiquadric: {
                Real e = config_.epsilon * r;
                Real s = std::sqrt(1.0 + e * e);
                return config_.epsilon * config_.epsilon * r / s;
            }
            case RBFType::InverseMultiquadric: {
                Real e = config_.epsilon * r;
                Real s = 1.0 + e * e;
                return -config_.epsilon * config_.epsilon * r / (s * std::sqrt(s));
            }
        }
        return 0.0;
    }

    // PHS 的二阶导数 d²/dr² φ(r)
    Real rbf_deriv2_r(Real r) const {
        if (r < 1e-15) r = 1e-15;
        switch (config_.rbf) {
            case RBFType::Polyharmonic3:
                return 6.0 * r;
            case RBFType::Polyharmonic5:
                return 20.0 * r * r * r;
            case RBFType::Polyharmonic7:
                return 42.0 * r * r * r * r * r;
            case RBFType::Gaussian: {
                Real e2 = config_.epsilon * config_.epsilon * r * r;
                return 2.0 * config_.epsilon * config_.epsilon * (2.0 * e2 - 1.0) * std::exp(-e2);
            }
            // MQ/IMQ 二阶导较复杂, MVP 用数值差分
            default: {
                Real h = 1e-6;
                return (rbf_deriv_r(r + h) - rbf_deriv_r(r - h)) / (2.0 * h);
            }
        }
    }

    // ============ 节点生成 ============
    // d 维张量积网格, 每维 n_per_dim 个等距点
    // x_i ∈ [-width * sigma_i * sqrt(T), +width * sigma_i * sqrt(T)]
    std::vector<std::vector<Real>> generate_tensor_nodes(
        const std::vector<Real>& sigmas, Real T, Real width) const {

        Size d = sigmas.size();
        Size n_per = config_.n_per_dim;
        std::vector<std::vector<Real>> nodes;

        // 每维的坐标
        std::vector<std::vector<Real>> coords(d);
        for (Size i = 0; i < d; ++i) {
            Real half = width * sigmas[i] * std::sqrt(T);
            coords[i].resize(n_per);
            for (Size j = 0; j < n_per; ++j) {
                coords[i][j] = -half + 2.0 * half * static_cast<Real>(j) / static_cast<Real>(n_per - 1);
            }
        }

        // 张量积
        if (d == 1) {
            for (Size i0 = 0; i0 < n_per; ++i0) {
                nodes.push_back({coords[0][i0]});
            }
        } else if (d == 2) {
            for (Size i0 = 0; i0 < n_per; ++i0) {
                for (Size i1 = 0; i1 < n_per; ++i1) {
                    nodes.push_back({coords[0][i0], coords[1][i1]});
                }
            }
        } else if (d == 3) {
            for (Size i0 = 0; i0 < n_per; ++i0) {
                for (Size i1 = 0; i1 < n_per; ++i1) {
                    for (Size i2 = 0; i2 < n_per; ++i2) {
                        nodes.push_back({coords[0][i0], coords[1][i1], coords[2][i2]});
                    }
                }
            }
        } else {
            throw std::invalid_argument("RBFFD: dimension > 3 not supported");
        }
        return nodes;
    }

    // ============ 识别边界节点 ============
    // 张量积网格中, 至少一维处于该维边界 (index 0 或 n_per-1) 的节点为边界节点
    // 1D: 2 个边界节点 (首尾)
    // 2D: 4*(n_per-1) 个边界节点 (四条边)
    // 3D: 6*(n_per-1)^2 + 12*(n_per-1) + 8 个边界节点 (六个面)
    void identify_boundary_nodes(
        const std::vector<std::vector<Real>>& nodes,
        Size n_per, Size d,
        std::vector<bool>& is_boundary,
        std::vector<Size>& boundary_indices) const {

        Size N = nodes.size();
        is_boundary.assign(N, false);

        if (d == 1) {
            is_boundary[0] = true;
            is_boundary[n_per - 1] = true;
        } else if (d == 2) {
            for (Size i0 = 0; i0 < n_per; ++i0) {
                for (Size i1 = 0; i1 < n_per; ++i1) {
                    Size idx = i0 * n_per + i1;
                    if (i0 == 0 || i0 == n_per - 1 || i1 == 0 || i1 == n_per - 1) {
                        is_boundary[idx] = true;
                    }
                }
            }
        } else if (d == 3) {
            for (Size i0 = 0; i0 < n_per; ++i0) {
                for (Size i1 = 0; i1 < n_per; ++i1) {
                    for (Size i2 = 0; i2 < n_per; ++i2) {
                        Size idx = (i0 * n_per + i1) * n_per + i2;
                        if (i0 == 0 || i0 == n_per - 1 ||
                            i1 == 0 || i1 == n_per - 1 ||
                            i2 == 0 || i2 == n_per - 1) {
                            is_boundary[idx] = true;
                        }
                    }
                }
            }
        }

        for (Size i = 0; i < N; ++i) {
            if (is_boundary[i]) boundary_indices.push_back(i);
        }
    }

    // ============ k-NN (暴力搜索) ============
    // 返回节点 i 的 stencil_size 个最近邻索引 (含 i 自身, 标准做法)
    // 必须包含 i 自身: 否则 L 对角线 = 0, 隐式时间步进矩阵 A = I - θdt L 条件数差
    std::vector<Size> find_knn(
        const std::vector<std::vector<Real>>& nodes,
        Size i, Size k) const {

        Size N = nodes.size();
        Size d = nodes[i].size();

        std::vector<std::pair<Real, Size>> dists;
        dists.reserve(N);
        for (Size j = 0; j < N; ++j) {
            Real dist = 0.0;
            for (Size dim = 0; dim < d; ++dim) {
                Real diff = nodes[i][dim] - nodes[j][dim];
                dist += diff * diff;
            }
            dists.push_back({dist, j});
        }
        std::partial_sort(dists.begin(), dists.begin() + std::min(k, dists.size()),
                          dists.end());

        std::vector<Size> result;
        result.reserve(k);
        for (Size m = 0; m < std::min(k, dists.size()); ++m) {
            result.push_back(dists[m].second);
        }
        return result;
    }

    // ============ Gauss-Jordan 矩阵求逆 (小矩阵) ============
    bool invert_matrix(std::vector<std::vector<Real>>& A, Size n) const {
        std::vector<std::vector<Real>> I(n, std::vector<Real>(n, 0.0));
        for (Size i = 0; i < n; ++i) I[i][i] = 1.0;

        for (Size i = 0; i < n; ++i) {
            Size pivot = i;
            Real max_val = std::abs(A[i][i]);
            for (Size k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > max_val) {
                    max_val = std::abs(A[k][i]);
                    pivot = k;
                }
            }
            if (max_val < 1e-15) return false;
            if (pivot != i) {
                std::swap(A[i], A[pivot]);
                std::swap(I[i], I[pivot]);
            }
            Real a_ii = A[i][i];
            for (Size j = 0; j < n; ++j) {
                A[i][j] /= a_ii;
                I[i][j] /= a_ii;
            }
            for (Size k = 0; k < n; ++k) {
                if (k == i) continue;
                Real factor = A[k][i];
                for (Size j = 0; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                    I[k][j] -= factor * I[i][j];
                }
            }
        }
        A = I;
        return true;
    }

    // ============ 多项式项数 ============
    // degree 0: 1 项 (常数)
    // degree 1: 1 + d 项
    // degree 2: 1 + d + d(d+1)/2 项
    Size n_poly_terms(Size d, Size degree) const {
        if (degree == 0) return 1;
        if (degree == 1) return 1 + d;
        // degree == 2
        return 1 + d + d * (d + 1) / 2;
    }

    // 多项式基在点 x 处的取值
    std::vector<Real> poly_basis(const std::vector<Real>& x, Size d, Size degree) const {
        std::vector<Real> p;
        if (degree >= 0) {
            p.push_back(1.0);  // 常数项
        }
        if (degree >= 1) {
            for (Size i = 0; i < d; ++i) p.push_back(x[i]);
        }
        if (degree >= 2) {
            for (Size i = 0; i < d; ++i) {
                for (Size j = i; j < d; ++j) {
                    p.push_back(x[i] * x[j]);
                }
            }
        }
        return p;
    }

    // 线性算子 L 作用于多项式基在 x_center 处的值
    // L = 0.5 Σ a_ij ∂²/∂x_i∂x_j + Σ b_i ∂/∂x_i - rate
    std::vector<Real> poly_operator(
        const std::vector<Real>& x_center,
        const std::vector<std::vector<Real>>& a_matrix,
        const std::vector<Real>& b_vector,
        Real rate, Size d, Size degree) const {

        std::vector<Real> Lp;

        if (degree >= 0) {
            // L[1] = -rate * 1
            Lp.push_back(-rate);
        }
        if (degree >= 1) {
            // L[x_i] = b_i * 1 - rate * x_i
            for (Size i = 0; i < d; ++i) {
                Lp.push_back(b_vector[i] - rate * x_center[i]);
            }
        }
        if (degree >= 2) {
            // L = 0.5 Σ_kl a_kl ∂²/∂x_k∂x_l + Σ_k b_k ∂/∂x_k - rate
            // 对 x_i x_j (i ≤ j, 单次遍历):
            //   ∂²(x_i x_j)/∂x_k∂x_l = δ_ik δ_jl + δ_jk δ_il
            //   0.5 Σ_kl a_kl (δ_ik δ_jl + δ_jk δ_il) = 0.5*(a_ij + a_ji) = a_ij (对称)
            //   Σ_k b_k ∂(x_i x_j)/∂x_k = b_i x_j + b_j x_i
            for (Size i = 0; i < d; ++i) {
                for (Size j = i; j < d; ++j) {
                    // i==j: 系数 a_ii; i≠j: 系数 a_ij (= a_ji, 仅算一次)
                    Real second_order = a_matrix[i][j];
                    Real Lp_val = second_order
                                + b_vector[i] * x_center[j]
                                + b_vector[j] * x_center[i]
                                - rate * x_center[i] * x_center[j];
                    Lp.push_back(Lp_val);
                }
            }
        }
        return Lp;
    }

    // RBF 在 x 处的值, 中心在 x_c
    Real rbf_at(const std::vector<Real>& x, const std::vector<Real>& x_c) const {
        Real r2 = 0.0;
        for (Size i = 0; i < x.size(); ++i) {
            Real diff = x[i] - x_c[i];
            r2 += diff * diff;
        }
        return rbf_eval(std::sqrt(r2));
    }

    // L 作用于 RBF φ(||x - x_c||) 在 x = x_eval 处的值
    // L = 0.5 Σ a_ij ∂²/∂x_i∂x_j + Σ b_i ∂/∂x_i - rate
    Real rbf_operator(const std::vector<Real>& x_eval,
                       const std::vector<Real>& x_c,
                       const std::vector<std::vector<Real>>& a_matrix,
                       const std::vector<Real>& b_vector,
                       Real rate) const {
        Size d = x_eval.size();
        Real r2 = 0.0;
        for (Size i = 0; i < d; ++i) {
            Real diff = x_eval[i] - x_c[i];
            r2 += diff * diff;
        }
        Real dist = std::sqrt(r2);
        if (dist < 1e-15) {
            // 在中心点, 用数值差分
            Real h = 1e-5;
            Real phi_center = rbf_at(x_eval, x_c);
            std::vector<Real> val_plus(d), val_minus(d);
            for (Size i = 0; i < d; ++i) {
                std::vector<Real> xp = x_eval; xp[i] += h;
                std::vector<Real> xm = x_eval; xm[i] -= h;
                val_plus[i] = rbf_at(xp, x_c);
                val_minus[i] = rbf_at(xm, x_c);
            }
            Real result = -rate * phi_center;  // -rate V 项
            // 一阶导: b_i * ∂φ/∂x_i
            for (Size i = 0; i < d; ++i) {
                result += b_vector[i] * (val_plus[i] - val_minus[i]) / (2.0 * h);
            }
            // 二阶导: 0.5 Σ_ij a_ij ∂²φ/∂x_i∂x_j
            // 双重循环遍历 (i,j) 和 (j,i), 由 0.5 因子自动处理对称性
            for (Size i = 0; i < d; ++i) {
                for (Size j = 0; j < d; ++j) {
                    if (i == j) {
                        result += 0.5 * a_matrix[i][i] * (val_plus[i] - 2.0 * phi_center + val_minus[i]) / (h * h);
                    } else {
                        std::vector<Real> xpp = x_eval; xpp[i] += h; xpp[j] += h;
                        std::vector<Real> xpm = x_eval; xpm[i] += h; xpm[j] -= h;
                        std::vector<Real> xmp = x_eval; xmp[i] -= h; xmp[j] += h;
                        std::vector<Real> xmm = x_eval; xmm[i] -= h; xmm[j] -= h;
                        Real d2xy = (rbf_at(xpp, x_c) - rbf_at(xpm, x_c)
                                   - rbf_at(xmp, x_c) + rbf_at(xmm, x_c)) / (4.0 * h * h);
                        result += 0.5 * a_matrix[i][j] * d2xy;
                    }
                }
            }
            return result;
        }

        // dist > 0: 解析公式
        // ∂φ/∂x_i = φ'(r) * (x_i - x_c_i) / r
        // ∂²φ/∂x_i∂x_j = φ''(r) * (x_i-x_c_i)(x_j-x_c_j)/r²
        //               + φ'(r) * (δ_ij r² - (x_i-x_c_i)(x_j-x_c_j)) / r³
        Real phi_r = rbf_eval(dist);
        Real phi1 = rbf_deriv_r(dist);
        Real phi2 = rbf_deriv2_r(dist);

        Real result = -rate * phi_r;  // -rate V 项

        // 一阶导项
        for (Size i = 0; i < d; ++i) {
            Real dxi = x_eval[i] - x_c[i];
            Real dphi_dxi = phi1 * dxi / dist;
            result += b_vector[i] * dphi_dxi;
        }

        // 二阶导项
        // PDE 项: 0.5 * Σ_ij a_ij ∂²/∂x_i∂x_j
        // 双重循环遍历所有 (i,j) 对 (含 (i,j) 和 (j,i)), 由 0.5 因子自动处理对称性
        // 故 aij 直接取 a_matrix[i][j] (无需乘 2)
        for (Size i = 0; i < d; ++i) {
            for (Size j = 0; j < d; ++j) {
                Real dxi = x_eval[i] - x_c[i];
                Real dxj = x_eval[j] - x_c[j];
                Real d2phi = phi2 * dxi * dxj / (dist * dist)
                           + phi1 * ((i == j ? 1.0 : 0.0) * dist * dist - dxi * dxj) / (dist * dist * dist);
                Real aij = a_matrix[i][j];
                result += 0.5 * aij * d2phi;
            }
        }
        return result;
    }

    // ============ 计算单个节点的 stencil 权重 ============
    // 对节点 i 的 stencil, 计算 L 在 x_i 处的权重向量 w (长度 = stencil_size)
    // 使得 L f(x_i) ≈ Σ w_j f(x_{stencil_j})
    std::vector<Real> compute_stencil_weights(
        const std::vector<std::vector<Real>>& nodes,
        Size i, const std::vector<Size>& stencil,
        const std::vector<std::vector<Real>>& a_matrix,
        const std::vector<Real>& b_vector,
        Real rate, Size d) const {

        Size k = stencil.size();
        Size np = n_poly_terms(d, config_.poly_degree);
        Size n_total = k + np;

        // 增广矩阵 A: (k+np) × (k+np)
        std::vector<std::vector<Real>> A(n_total, std::vector<Real>(n_total, 0.0));

        for (Size m = 0; m < k; ++m) {
            for (Size n = 0; n < k; ++n) {
                A[m][n] = rbf_at(nodes[stencil[m]], nodes[stencil[n]]);
            }
            auto pm = poly_basis(nodes[stencil[m]], d, config_.poly_degree);
            for (Size n = 0; n < np; ++n) {
                A[m][k + n] = pm[n];
                A[k + n][m] = pm[n];
            }
        }

        // 右端项 b: (k+np)
        std::vector<Real> b(n_total, 0.0);
        for (Size m = 0; m < k; ++m) {
            b[m] = rbf_operator(nodes[i], nodes[stencil[m]], a_matrix, b_vector, rate);
        }
        auto Lp = poly_operator(nodes[i], a_matrix, b_vector, rate, d, config_.poly_degree);
        for (Size n = 0; n < np; ++n) {
            b[k + n] = Lp[n];
        }

        // 解 A w = b (Gauss-Jordan 消元)
        std::vector<std::vector<Real>> A_copy = A;
        for (Size iter = 0; iter < n_total; ++iter) {
            Size pivot = iter;
            Real max_val = std::abs(A_copy[iter][iter]);
            for (Size row = iter + 1; row < n_total; ++row) {
                if (std::abs(A_copy[row][iter]) > max_val) {
                    max_val = std::abs(A_copy[row][iter]);
                    pivot = row;
                }
            }
            if (max_val < 1e-15) {
                return std::vector<Real>(k, 0.0);
            }
            if (pivot != iter) {
                std::swap(A_copy[iter], A_copy[pivot]);
                std::swap(b[iter], b[pivot]);
            }
            Real a_ii = A_copy[iter][iter];
            for (Size j = 0; j < n_total; ++j) A_copy[iter][j] /= a_ii;
            b[iter] /= a_ii;
            for (Size row = 0; row < n_total; ++row) {
                if (row == iter) continue;
                Real factor = A_copy[row][iter];
                for (Size j = 0; j < n_total; ++j) A_copy[row][j] -= factor * A_copy[iter][j];
                b[row] -= factor * b[iter];
            }
        }

        std::vector<Real> weights(k);
        for (Size m = 0; m < k; ++m) weights[m] = b[m];
        return weights;
    }

    // ============ 稀疏矩阵向量乘 ============
    // CSR 格式: row_ptr, col_idx, values
    // y = A x
    static void spmv_csr(const std::vector<Size>& row_ptr,
                          const std::vector<Size>& col_idx,
                          const std::vector<Real>& values,
                          const std::vector<Real>& x,
                          std::vector<Real>& y) {
        Size n = row_ptr.size() - 1;
        for (Size i = 0; i < n; ++i) {
            Real sum = 0.0;
            for (Size idx = row_ptr[i]; idx < row_ptr[i + 1]; ++idx) {
                sum += values[idx] * x[col_idx[idx]];
            }
            y[i] = sum;
        }
    }

    // ============ BiCGSTAB 求解器 ============
    // 求解 A x = b, A 为 CSR 稀疏矩阵
    bool bicgstab(const std::vector<Size>& row_ptr,
                   const std::vector<Size>& col_idx,
                   const std::vector<Real>& values,
                   const std::vector<Real>& b,
                   std::vector<Real>& x) const {
        Size n = b.size();
        x.assign(n, 0.0);

        std::vector<Real> r(n), rhat(n), p(n), v(n), s(n), t(n), ph(n);
        spmv_csr(row_ptr, col_idx, values, x, r);
        for (Size i = 0; i < n; ++i) r[i] = b[i] - r[i];
        rhat = r;

        Real r0_norm = 0.0;
        for (Size i = 0; i < n; ++i) r0_norm += r[i] * r[i];
        r0_norm = std::sqrt(r0_norm);
        if (r0_norm < config_.solver_tol) return true;

        Real rho = 1.0, alpha = 1.0, omega = 1.0;
        std::fill(p.begin(), p.end(), 0.0);
        std::fill(v.begin(), v.end(), 0.0);

        for (Size iter = 0; iter < config_.max_iter; ++iter) {
            Real rho_new = 0.0;
            for (Size i = 0; i < n; ++i) rho_new += rhat[i] * r[i];
            if (std::abs(rho_new) < 1e-30) return false;

            Real beta = (rho_new / rho) * (alpha / omega);
            rho = rho_new;

            for (Size i = 0; i < n; ++i) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }

            spmv_csr(row_ptr, col_idx, values, p, v);
            Real pv = 0.0;
            for (Size i = 0; i < n; ++i) pv += rhat[i] * v[i];
            if (std::abs(pv) < 1e-30) return false;
            alpha = rho / pv;

            for (Size i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
            Real s_norm = 0.0;
            for (Size i = 0; i < n; ++i) s_norm += s[i] * s[i];
            s_norm = std::sqrt(s_norm);
            if (s_norm < config_.solver_tol) {
                for (Size i = 0; i < n; ++i) x[i] += alpha * p[i];
                return true;
            }

            spmv_csr(row_ptr, col_idx, values, s, t);
            Real ts = 0.0, tt = 0.0;
            for (Size i = 0; i < n; ++i) { ts += t[i] * s[i]; tt += t[i] * t[i]; }
            if (std::abs(tt) < 1e-30) return false;
            omega = ts / tt;

            for (Size i = 0; i < n; ++i) {
                x[i] += alpha * p[i] + omega * s[i];
                r[i] = s[i] - omega * t[i];
            }

            Real r_norm = 0.0;
            for (Size i = 0; i < n; ++i) r_norm += r[i] * r[i];
            r_norm = std::sqrt(r_norm);
            if (r_norm / r0_norm < config_.solver_tol) return true;
        }
        return false;
    }

    // ============ 线性插值: 在节点集中找到 S0 对应的 V 值 ============
    Real interpolate_at_s0(
        const std::vector<std::vector<Real>>& nodes,
        const std::vector<Real>& V,
        const std::vector<Real>& S0,
        const std::vector<Real>& K_vec,
        const std::vector<Real>& sigmas, Real T) const {

        Size d = S0.size();
        // 将 S0 转换到 log-S 空间: x_i = ln(S0_i / K_i)
        std::vector<Real> x0(d);
        for (Size i = 0; i < d; ++i) {
            x0[i] = std::log(S0[i] / K_vec[i]);
        }

        // 最近邻插值 (简化, MVP)
        Size nearest = 0;
        Real min_dist = std::numeric_limits<Real>::max();
        for (Size n = 0; n < nodes.size(); ++n) {
            Real dist = 0.0;
            for (Size i = 0; i < d; ++i) {
                Real diff = nodes[n][i] - x0[i];
                dist += diff * diff;
            }
            if (dist < min_dist) {
                min_dist = dist;
                nearest = n;
            }
        }

        // d=1 时用线性插值, d>=2 用逆距离加权
        if (d == 1) {
            // 找到 x0 两侧的节点
            Size n_below = 0, n_above = nodes.size() - 1;
            for (Size n = 0; n < nodes.size(); ++n) {
                if (nodes[n][0] <= x0[0] && n > n_below) n_below = n;
                if (nodes[n][0] >= x0[0] && n < n_above) n_above = n;
            }
            if (n_below == n_above) return V[n_below];
            Real w = (x0[0] - nodes[n_below][0]) / (nodes[n_above][0] - nodes[n_below][0]);
            return (1.0 - w) * V[n_below] + w * V[n_above];
        }

        // d>=2: 逆距离加权 (k 个最近邻)
        Size k = std::min(static_cast<Size>(8), static_cast<Size>(nodes.size()));
        std::vector<std::pair<Real, Size>> dists;
        dists.reserve(nodes.size());
        for (Size n = 0; n < nodes.size(); ++n) {
            Real dist = 0.0;
            for (Size i = 0; i < d; ++i) {
                Real diff = nodes[n][i] - x0[i];
                dist += diff * diff;
            }
            dists.push_back({std::sqrt(dist), n});
        }
        std::partial_sort(dists.begin(), dists.begin() + k, dists.end());

        Real weight_sum = 0.0, val_sum = 0.0;
        for (Size m = 0; m < k; ++m) {
            Real w = 1.0 / (dists[m].first + 1e-10);
            weight_sum += w;
            val_sum += w * V[dists[m].second];
        }
        return val_sum / weight_sum;
    }
};

}  // namespace v1
}  // namespace cpphub
