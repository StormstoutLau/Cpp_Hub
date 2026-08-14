#pragma once
// SOURCE: Longstaff & Schwartz (2001) "Valuing American Options by Simulation"
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.8
// 模块: Longstaff-Schwartz 蒙特卡洛 (LSMC) 美式/百慕大期权定价引擎
// 算法:
//   1. 生成 N 条前向路径，在行使点 t_1, ..., t_K 记录标的价格
//   2. 从到期日 T 反向迭代:
//      a. 计算内在价值 h_k = payoff(S_{t_k})
//      b. 对 ITM 路径 (h_k > 0), 回归条件期望 E[continuation | S_{t_k}]
//         - 回归基函数: Laguerre / Hermite / Monomial
//         - OLS: β = (X^T X + λI)^{-1} X^T Y
//      c. 若 h_k > continuation_estimate, 标记行使并停止该路径
//   3. 每条路径取最早行使点，贴现收益
//   4. 价格 = mean(discounted payoffs)
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/instruments/payoff/payoff.hpp"
#include "cpphub/pricing/monte_carlo/path_generator.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <numeric>

namespace cpphub {
inline namespace v1 {

enum class BasisType {
    Laguerre,   // Longstaff-Schwartz (2001) 经典选择
    Hermite,    // 概率学家 Hermite 多项式
    Monomial,   // 简单幂基 (1, x, x^2, ...)
    Chebyshev   // 第一类 Chebyshev 多项式
};

// RISK-007: K-fold 交叉验证配置
struct CVConfig {
    Size k_fold = 5;                                              // K-fold 数
    std::vector<Real> lambda_grid = {0.0, 0.01, 0.1, 1.0, 10.0};  // λ 候选
    uint64_t cv_seed = 12345;                                     // 分折随机种子
};

struct LSMCConfig {
    Real S0 = 100.0;
    Real K = 100.0;
    Real T = 1.0;
    Real r = 0.05;
    Real q = 0.0;
    Real sigma = 0.20;
    Size n_paths = 10000;       // 路径数 (建议 >= 10000)
    Size n_steps = 50;          // 行使点数 (不含 t=0)
    BasisType basis = BasisType::Laguerre;
    Size basis_order = 3;       // 基函数阶数 (建议 3-5)
    uint64_t seed = 42;
    bool antithetic = true;     // 反变量方差缩减
    Real ridge_lambda = 0.0;    // Ridge 正则化 (0=纯 OLS, >0 防 overfitting)
    PathScheme path_scheme = PathScheme::Exact;
    // RISK-007: 交叉验证配置
    CVConfig cv_config;                    // K-fold CV 配置
    bool use_cross_validation = false;     // false=用 ridge_lambda, true=用 CV 选 λ
};

struct LSMCResult {
    Real price = 0.0;           // 美式期权价格
    Real std_error = 0.0;       // 蒙特卡洛标准误差
    Real european_price = 0.0;  // 对应欧式期权价格 (参考)
    Real early_exercise_premium = 0.0;  // 美式 - 欧式
    Size n_exercised = 0;       // 实际行使的路径数
    Size n_paths = 0;           // 总路径数
    // RISK-007: 每个行使时点 CV 选择的 λ (若 use_cross_validation=true)
    std::vector<Real> selected_lambdas;
};

// Laguerre 多项式 L_n(x) (物理学家版本)
// L_0(x) = 1, L_1(x) = 1 - x, L_n(x) = ((2n-1-x)*L_{n-1}(x) - (n-1)*L_{n-2}(x)) / n
inline Real laguerre(Size n, Real x) {
    if (n == 0) return 1.0;
    if (n == 1) return 1.0 - x;
    Real L_prev2 = 1.0;
    Real L_prev1 = 1.0 - x;
    Real L_curr = 0.0;
    for (Size k = 2; k <= n; ++k) {
        L_curr = ((2.0 * static_cast<Real>(k) - 1.0 - x) * L_prev1
                  - (static_cast<Real>(k) - 1.0) * L_prev2) / static_cast<Real>(k);
        L_prev2 = L_prev1;
        L_prev1 = L_curr;
    }
    return L_curr;
}

// 概率学家 Hermite 多项式 He_n(x)
// He_0(x) = 1, He_1(x) = x, He_n(x) = x*He_{n-1}(x) - (n-1)*He_{n-2}(x)
inline Real hermite_prob(Size n, Real x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    Real H_prev2 = 1.0;
    Real H_prev1 = x;
    Real H_curr = 0.0;
    for (Size k = 2; k <= n; ++k) {
        H_curr = x * H_prev1 - (static_cast<Real>(k) - 1.0) * H_prev2;
        H_prev2 = H_prev1;
        H_prev1 = H_curr;
    }
    return H_curr;
}

// 第一类 Chebyshev 多项式 T_n(x) (定义域 [-1, 1])
// T_0(x) = 1, T_1(x) = x, T_n(x) = 2x*T_{n-1}(x) - T_{n-2}(x)
inline Real chebyshev_t(Size n, Real x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    Real T_prev2 = 1.0;
    Real T_prev1 = x;
    Real T_curr = 0.0;
    for (Size k = 2; k <= n; ++k) {
        T_curr = 2.0 * x * T_prev1 - T_prev2;
        T_prev2 = T_prev1;
        T_prev1 = T_curr;
    }
    return T_curr;
}

// 计算基函数向量 [phi_0(x), phi_1(x), ..., phi_{m-1}(x)]
inline std::vector<Real> basis_eval(BasisType type, Size order, Real x) {
    std::vector<Real> basis(order);
    for (Size n = 0; n < order; ++n) {
        switch (type) {
            case BasisType::Laguerre:
                basis[n] = laguerre(n, x);
                break;
            case BasisType::Hermite:
                basis[n] = hermite_prob(n, x);
                break;
            case BasisType::Monomial:
                basis[n] = (n == 0) ? 1.0 : std::pow(x, static_cast<Real>(n));
                break;
            case BasisType::Chebyshev:
                basis[n] = chebyshev_t(n, x);
                break;
        }
    }
    return basis;
}

// 解线性方程组 A x = b (A 为 n x n 矩阵, Gauss-Jordan 消元)
// 返回是否成功 (矩阵奇异时返回 false)
inline bool solve_linear_system(std::vector<std::vector<Real>>& A,
                                 std::vector<Real>& b, Size n) {
    for (Size i = 0; i < n; ++i) {
        // 寻找主元
        Size pivot = i;
        Real max_val = std::abs(A[i][i]);
        for (Size k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > max_val) {
                max_val = std::abs(A[k][i]);
                pivot = k;
            }
        }
        if (max_val < 1e-15) return false;  // 奇异矩阵
        if (pivot != i) {
            std::swap(A[i], A[pivot]);
            std::swap(b[i], b[pivot]);
        }
        // 消元
        for (Size k = i + 1; k < n; ++k) {
            Real factor = A[k][i] / A[i][i];
            for (Size j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    // 回代
    for (Size i = n; i-- > 0;) {
        for (Size j = i + 1; j < n; ++j) {
            b[i] -= A[i][j] * b[j];
        }
        b[i] /= A[i][i];
    }
    return true;
}

class LSMCEngine {
public:
    explicit LSMCEngine(LSMCConfig cfg) : cfg_(cfg) {
        validate_config();
    }

    // 美式期权定价 (所有行使点均可行使)
    LSMCResult price_american(const PayOff& payoff) const {
        return price_bermudan(payoff, all_exercise_times());
    }

    // 百慕大期权定价 (指定行使时间点)
    // exercise_times: 行使时间点向量 (0 < t_1 < t_2 < ... < T)
    LSMCResult price_bermudan(const PayOff& payoff,
                              const std::vector<Real>& exercise_times) const {
        if (exercise_times.empty()) {
            throw std::invalid_argument("LSMCEngine: exercise_times must not be empty");
        }
        cv_selected_lambdas_.clear();  // RISK-007: 重置 CV 记录
        // 生成路径
        GBMConfig gbm_cfg{cfg_.S0, cfg_.r, cfg_.q, cfg_.sigma, cfg_.T, cfg_.n_steps};
        GBMPathGenerator gen(gbm_cfg, cfg_.path_scheme);
        std::vector<Real> paths = gen.generate_paths(cfg_.n_paths, cfg_.seed, cfg_.antithetic);
        Size path_len = gen.path_length();  // n_steps + 1

        // 将行使时间转换为路径索引
        // 路径索引 k 对应时间 t_k = k * dt, dt = T / n_steps
        std::vector<Size> exercise_indices;
        exercise_indices.reserve(exercise_times.size());
        Real dt = cfg_.T / static_cast<Real>(cfg_.n_steps);
        for (Real t : exercise_times) {
            Size idx = static_cast<Size>(std::round(t / dt));
            if (idx == 0) idx = 1;  // t=0 不行使
            if (idx > cfg_.n_steps) idx = cfg_.n_steps;
            exercise_indices.push_back(idx);
        }
        // 去重并排序
        std::sort(exercise_indices.begin(), exercise_indices.end());
        exercise_indices.erase(std::unique(exercise_indices.begin(), exercise_indices.end()),
                               exercise_indices.end());

        Size n_paths = cfg_.n_paths;
        Size n_ex = exercise_indices.size();
        Real dt_step = cfg_.T / static_cast<Real>(cfg_.n_steps);

        // 每条路径的现金流 (贴现到 t=0)
        std::vector<Real> cashflows(n_paths, 0.0);
        // 每条路径的行使时间索引 (n_ex 表示未行使)
        std::vector<Size> exercise_step(n_paths, n_ex);
        // 到期日: 所有路径行使 payoff(S_T)
        Size T_idx = exercise_indices.back();
        for (Size p = 0; p < n_paths; ++p) {
            Real S_T = paths[p * path_len + T_idx];
            Real h = payoff(S_T);
            cashflows[p] = h * std::exp(-cfg_.r * cfg_.T);
            if (h > 0.0) exercise_step[p] = n_ex - 1;  // 标记在最后一个点行使
        }

        // 反向迭代: 从倒数第二个行使点到第一个
        for (Size ei = n_ex - 1; ei-- > 0;) {
            Size k = exercise_indices[ei];  // 当前行使点的路径索引
            Real t_k = static_cast<Real>(k) * dt_step;
            Real t_next = static_cast<Real>(exercise_indices[ei + 1]) * dt_step;
            Real dt_disc = t_next - t_k;  // 到下一个行使点的时间间隔
            Real disc = std::exp(-cfg_.r * dt_disc);

            // 收集 ITM 路径
            std::vector<Size> itm_indices;
            std::vector<Real> S_itm;
            std::vector<Real> Y_itm;  // 贴现的未来现金流
            for (Size p = 0; p < n_paths; ++p) {
                if (exercise_step[p] > ei) {  // 尚未行使
                    Real S = paths[p * path_len + k];
                    Real h = payoff(S);
                    if (h > 0.0) {  // ITM
                        itm_indices.push_back(p);
                        S_itm.push_back(S);
                        // Y = 当前 cashflows[p] 贴现到 t_k
                        // cashflows[p] 是贴现到 t=0 的值, 需要先 forward 到 t_next, 再 disc 到 t_k
                        // 但 cashflows[p] 已经是 t=0 贴现值, forward 到 t_k: * exp(r * t_k)
                        // 实际上: continuation value at t_k = E[discounted future CF | S_{t_k}]
                        // 未来 CF 在 t_next 的值 = cashflows[p] * exp(r * t_next)
                        // 贴现到 t_k = * exp(-r * dt_disc) = cashflows[p] * exp(r * t_next) * exp(-r * dt_disc)
                        //             = cashflows[p] * exp(r * t_k)
                        Y_itm.push_back(cashflows[p] * std::exp(cfg_.r * t_next) * disc);
                    }
                }
            }

            if (itm_indices.size() < cfg_.basis_order + 1) continue;  // 样本不足, 跳过回归

            // 构造回归矩阵 X (n_itm x m) 和 Y (n_itm)
            Size n_itm = itm_indices.size();
            Size m = cfg_.basis_order;
            // 归一化: x = S / K (Longstaff-Schwartz 2001)
            std::vector<std::vector<Real>> X(n_itm, std::vector<Real>(m));
            for (Size i = 0; i < n_itm; ++i) {
                Real x = S_itm[i] / cfg_.K;
                std::vector<Real> phi = basis_eval(cfg_.basis, m, x);
                for (Size j = 0; j < m; ++j) {
                    X[i][j] = phi[j];
                }
            }

            // 正规方程: (X^T X + λI) β = X^T Y
            std::vector<std::vector<Real>> XtX(m, std::vector<Real>(m, 0.0));
            std::vector<Real> XtY(m, 0.0);
            for (Size i = 0; i < n_itm; ++i) {
                for (Size j = 0; j < m; ++j) {
                    XtY[j] += X[i][j] * Y_itm[i];
                    for (Size k = 0; k < m; ++k) {
                        XtX[j][k] += X[i][j] * X[i][k];
                    }
                }
            }
            // Ridge 正则化: 选择 λ (固定或 CV)
            Real lambda_used = cfg_.ridge_lambda;
            if (cfg_.use_cross_validation) {
                lambda_used = select_lambda_cv(X, Y_itm, m);
                cv_selected_lambdas_.push_back(lambda_used);
            }
            if (lambda_used > 0.0) {
                for (Size j = 0; j < m; ++j) {
                    XtX[j][j] += lambda_used;
                }
            }

            std::vector<Real> beta = XtY;
            bool ok = solve_linear_system(XtX, beta, m);
            if (!ok) continue;  // 回归失败, 跳过此行使点

            // 对 ITM 路径估计 continuation value 并决策
            for (Size i = 0; i < n_itm; ++i) {
                Size p = itm_indices[i];
                Real x = S_itm[i] / cfg_.K;
                std::vector<Real> phi = basis_eval(cfg_.basis, m, x);
                Real continuation = 0.0;
                for (Size j = 0; j < m; ++j) {
                    continuation += beta[j] * phi[j];
                }
                Real h = payoff(S_itm[i]);
                if (h > continuation) {
                    // 行使: 更新该路径的现金流
                    cashflows[p] = h * std::exp(-cfg_.r * t_k);
                    exercise_step[p] = ei;
                }
                // 否则保持原有现金流 (continuation)
            }
        }

        // 计算最终价格和统计量
        Real sum = 0.0, sum_sq = 0.0;
        Size n_exercised = 0;
        for (Size p = 0; p < n_paths; ++p) {
            sum += cashflows[p];
            sum_sq += cashflows[p] * cashflows[p];
            if (exercise_step[p] < n_ex) ++n_exercised;
        }
        Real mean = sum / static_cast<Real>(n_paths);
        Real var = (sum_sq / static_cast<Real>(n_paths) - mean * mean);
        if (var < 0.0) var = 0.0;
        Real std_err = std::sqrt(var / static_cast<Real>(n_paths));

        // 欧式参考价格
        Real euro_price = 0.0;
        if (cfg_.q == 0.0 && payoff.name() == "Put") {
            euro_price = bsm_put_price(cfg_.S0, cfg_.K, cfg_.T, cfg_.r, cfg_.q, cfg_.sigma);
        } else if (cfg_.q == 0.0 && payoff.name() == "Call") {
            euro_price = bsm_call_price(cfg_.S0, cfg_.K, cfg_.T, cfg_.r, cfg_.q, cfg_.sigma);
        } else {
            // 用 MC 欧式价格 (取 S_T 的 payoff 贴现均值)
            Real euro_sum = 0.0;
            for (Size p = 0; p < n_paths; ++p) {
                Real S_T = paths[p * path_len + cfg_.n_steps];
                euro_sum += payoff(S_T) * std::exp(-cfg_.r * cfg_.T);
            }
            euro_price = euro_sum / static_cast<Real>(n_paths);
        }

        LSMCResult result;
        result.price = mean;
        result.std_error = std_err;
        result.european_price = euro_price;
        result.early_exercise_premium = mean - euro_price;
        result.n_exercised = n_exercised;
        result.n_paths = n_paths;
        result.selected_lambdas = std::move(cv_selected_lambdas_);  // RISK-007
        return result;
    }

    const LSMCConfig& config() const { return cfg_; }

private:
    LSMCConfig cfg_;
    mutable std::vector<Real> cv_selected_lambdas_;  // RISK-007: CV 选择的 λ 记录

    // RISK-007: K-fold 交叉验证选择最优 λ
    // 返回最优 λ; 若样本不足 (n_itm < k_fold * basis_order) 则 fallback 到 ridge_lambda
    Real select_lambda_cv(const std::vector<std::vector<Real>>& X,
                            const std::vector<Real>& Y, Size m) const {
        Size n_itm = X.size();
        Size k = cfg_.cv_config.k_fold;

        // Fallback: 样本不足以分折
        if (n_itm < k * m) {
            return cfg_.ridge_lambda;
        }

        // 生成分折索引 (随机排列后均分)
        std::vector<Size> perm(n_itm);
        std::iota(perm.begin(), perm.end(), 0);
        Philox4x64 cv_rng(cfg_.cv_config.cv_seed);
        // Fisher-Yates 洗牌 (用 Philox 代替 std::shuffle 保证跨平台一致)
        for (Size i = n_itm - 1; i > 0; --i) {
            Size j = static_cast<Size>(cv_rng() % (i + 1));
            std::swap(perm[i], perm[j]);
        }

        // 每折大小 (最后一个折吸收余数)
        Size fold_size = n_itm / k;
        Size fold_remainder = n_itm % k;

        Real best_mse = std::numeric_limits<Real>::max();
        Real best_lambda = cfg_.ridge_lambda;

        for (Real lambda : cfg_.cv_config.lambda_grid) {
            Real total_mse = 0.0;
            Size fold_start = 0;
            for (Size fold = 0; fold < k; ++fold) {
                Size this_fold_size = fold_size + (fold < fold_remainder ? 1 : 0);
                Size fold_end = fold_start + this_fold_size;

                // 训练集: perm[0..fold_start) + perm[fold_end..n_itm)
                // 验证集: perm[fold_start..fold_end)
                std::vector<std::vector<Real>> X_train;
                std::vector<Real> Y_train;
                for (Size idx = 0; idx < n_itm; ++idx) {
                    if (idx >= fold_start && idx < fold_end) continue;
                    X_train.push_back(X[perm[idx]]);
                    Y_train.push_back(Y[perm[idx]]);
                }

                // 训练: (X_train^T X_train + λI) β = X_train^T Y_train
                std::vector<std::vector<Real>> XtX_train(m, std::vector<Real>(m, 0.0));
                std::vector<Real> XtY_train(m, 0.0);
                for (Size i = 0; i < X_train.size(); ++i) {
                    for (Size j = 0; j < m; ++j) {
                        XtY_train[j] += X_train[i][j] * Y_train[i];
                        for (Size jj = 0; jj < m; ++jj) {
                            XtX_train[j][jj] += X_train[i][j] * X_train[i][jj];
                        }
                    }
                }
                if (lambda > 0.0) {
                    for (Size j = 0; j < m; ++j) XtX_train[j][j] += lambda;
                }
                std::vector<Real> beta = XtY_train;
                bool ok = solve_linear_system(XtX_train, beta, m);
                if (!ok) {
                    total_mse = std::numeric_limits<Real>::max();
                    break;
                }

                // 验证集 MSE
                for (Size idx = fold_start; idx < fold_end; ++idx) {
                    Real pred = 0.0;
                    for (Size j = 0; j < m; ++j) {
                        pred += beta[j] * X[perm[idx]][j];
                    }
                    Real err = Y[perm[idx]] - pred;
                    total_mse += err * err;
                }
                fold_start = fold_end;
            }

            if (total_mse < best_mse) {
                best_mse = total_mse;
                best_lambda = lambda;
            }
        }
        return best_lambda;
    }

    void validate_config() const {
        if (cfg_.S0 <= 0.0) throw std::invalid_argument("LSMCEngine: S0 must be positive");
        if (cfg_.K <= 0.0) throw std::invalid_argument("LSMCEngine: K must be positive");
        if (cfg_.T <= 0.0) throw std::invalid_argument("LSMCEngine: T must be positive");
        if (cfg_.sigma < 0.0) throw std::invalid_argument("LSMCEngine: sigma must be non-negative");
        if (cfg_.n_paths < 100) throw std::invalid_argument("LSMCEngine: n_paths must be >= 100");
        if (cfg_.n_steps == 0) throw std::invalid_argument("LSMCEngine: n_steps must be positive");
        if (cfg_.basis_order == 0) throw std::invalid_argument("LSMCEngine: basis_order must be positive");
    }

    std::vector<Real> all_exercise_times() const {
        std::vector<Real> times;
        Real dt = cfg_.T / static_cast<Real>(cfg_.n_steps);
        times.reserve(cfg_.n_steps);
        for (Size k = 1; k <= cfg_.n_steps; ++k) {
            times.push_back(static_cast<Real>(k) * dt);
        }
        return times;
    }
};

}  // namespace v1
}  // namespace cpphub
