// =============================================================================
// har_model.hpp
// Phase 5 v1.4.2 Wave C - Heterogeneous Autoregressive (HAR) Model
//
// R 对照: HARmodel(data, periods=c(1,5,22), periodsJ=c(1,5,22), periodsQ=c(1),
//                  leverage=NULL, RVest=c("rCov","rBPCov","rQuar"),
//                  type="HAR", inputType="RM", jumpTest="ABDJumptest",
//                  alpha=0.05, h=1, transform=NULL, externalRegressor=NULL,
//                  periodsExternal=c(1))
//
// 文献: Corsi (2009), JFE 4(2), 174-196, doi:10.1093/jjfinec/nbp001
//       Andersen, Bollerslev, Diebold (2007), Restat 89, 701-720
//       Bollerslev, Patton, Quaedvlieg (2016), J. Econometrics 192, 1-18 (BPQ RQ 变换)
//
// SOURCE: PHASE5_HFE_SPEC §5.2, §5.3 D7/D9, §5.4
//   R highfrequency 1.0.3 src/HARmodel.cpp har_agg (L7-21)
//   R highfrequency 1.0.3 R/HARmodel.R HARmodel (L189-514)
//   R highfrequency 1.0.3 R/HARmodel.R estimhar (L5-16)
//   R highfrequency 1.0.3 R/HARmodel.R harInsanityFilter (L49-52)
//
// 关键幻觉排除 (spec §5.3 D7/D9):
//   D7: RQ 变换用 BPQ 2016 (sqrt(RQ) - sqrt(mean(RM3))), 非 Corsi 2009 原始 RQ
//   D9: har_agg 索引边界 [j-p, j-1] 含 j-1 不含 j (R 1-based → C++ 0-based 需调整)
//
// 实现范围:
//   - type="HAR":   基础 HAR-RV (Corsi 2009)
//   - type="HARJ":  HAR-RV-J (含跳跃, Andersen et al. 2007)
//   - type="HARQ":  HAR-Q (含四阶矩, BPQ 2016)
//   - type="CHAR":  CHAR (连续成分, Andersen et al. 2007)
//   - type="HARCJ": HAR-CJ (条件跳跃)
//   - type="HARQJ": HAR-QJ (四阶矩+跳跃)
//   - type="CHARQ": CHARQ (连续成分+四阶矩)
//   - transform: "log", "sqrt", NULL
//   - h: 多步聚合
//   - leverage: 杠杆效应 (需要收益率输入)
//
// 算法 (R 源码实测):
//   1. har_agg(RM, periods): 滚动窗口平均
//      mHARData[row][i] = mean(RM[row-p+1 .. row]), row = p-1..iT-1
//   2. y = har_agg(RM1, h, 1)[(maxp+h):n]
//   3. x = RVmatrix1[(maxp:(n-h)), ]
//   4. 根据 type 添加 J/RQ/C 组件
//   5. OLS: beta = (X^T X)^{-1} X^T y
//   6. harInsanityFilter: fitted < lower | fitted > upper → replacement = mean(y)
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <functional>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// HAR 模型类型枚举
// =============================================================================
enum class HarType {
    HAR,     // 基础 HAR-RV
    HARJ,    // HAR-RV-J (含跳跃)
    HARCJ,   // HAR-CJ (条件跳跃)
    HARQ,    // HAR-Q (含四阶矩)
    HARQJ,   // HAR-QJ (四阶矩+跳跃)
    CHAR,    // CHAR (连续成分)
    CHARQ    // CHARQ (连续成分+四阶矩)
};

// =============================================================================
// HAR 模型变换类型
// =============================================================================
enum class HarTransform {
    None,    // 无变换
    Log,     // 对数变换
    Sqrt     // 平方根变换
};

// =============================================================================
// HAR 模型估计结果
// =============================================================================
struct HarModelResult {
    std::vector<Real> coefficients;   // OLS 系数 (含常数项)
    std::vector<std::string> coef_names;  // 系数名称
    std::vector<Real> fitted_values;  // 拟合值 (经 insanity filter)
    std::vector<Real> residuals;      // 残差
    Real r_squared;                   // R^2
    Real adj_r_squared;               // 调整 R^2
    Real llh;                         // 对数似然
    HarType type;                     // 模型类型
    HarTransform transform;           // 变换类型
    Size h;                           // 预测步长
    Size n_obs;                       // 观测数 (估计样本大小)
    Size maxp;                        // 最大聚合阶数
    std::vector<Real> y;              // 因变量 (变换后)
    std::vector<std::vector<Real>> x; // 设计矩阵 (不含常数列)
};

// =============================================================================
// har_agg: 滚动窗口聚合 (R har_agg 对标)
//
// R 源码 (HARmodel.cpp L7-21):
//   for (i in 0..nperiods-1) {
//     for (j in periods[i]..iT) {  // j 是 1-based
//       mHARData[j-1][i] = sum(RM[j-periods[i] .. j-1]) / periods[i]
//     }
//   }
//
// C++ 0-based 等价:
//   for (i in 0..nperiods-1) {
//     p = periods[i]
//     for (row in p-1..iT-1) {
//       mHARData[row][i] = sum(RM[row-p+1 .. row]) / p
//     }
//   }
//
// 输入: RM (长度 iT), periods (nperiods 个聚合周期)
// 输出: iT x nperiods 矩阵 (row-major), 前 p-1 行为 NaN
// =============================================================================
inline std::vector<std::vector<Real>> har_agg(
        const std::vector<Real>& RM,
        const std::vector<Size>& periods) {
    const Size iT = RM.size();
    const Size nperiods = periods.size();
    std::vector<std::vector<Real>> mHARData(iT,
        std::vector<Real>(nperiods, std::numeric_limits<Real>::quiet_NaN()));

    for (Size i = 0; i < nperiods; ++i) {
        const Size p = periods[i];
        if (p == 0 || p > iT) continue;
        // row 从 p-1 开始 (0-based), 对应 R 的 j=p (1-based)
        for (Size row = p - 1; row < iT; ++row) {
            Real s = 0.0;
            // RM[row-p+1 .. row], 共 p 个元素
            for (Size k = row - p + 1; k <= row; ++k) s += RM[k];
            mHARData[row][i] = s / static_cast<Real>(p);
        }
    }
    return mHARData;
}

// =============================================================================
// har_agg_single: 单周期聚合 (h 参数用)
// =============================================================================
inline std::vector<Real> har_agg_single(
        const std::vector<Real>& RM, Size period) {
    const Size iT = RM.size();
    std::vector<Real> result(iT, std::numeric_limits<Real>::quiet_NaN());
    if (period == 0 || period > iT) return result;
    for (Size row = period - 1; row < iT; ++row) {
        Real s = 0.0;
        for (Size k = row - period + 1; k <= row; ++k) s += RM[k];
        result[row] = s / static_cast<Real>(period);
    }
    return result;
}

// =============================================================================
// harInsanityFilter: BPQ insanity filter
//
// R 源码 (HARmodel.R L49-52):
//   fittedValues[(fittedValues < lower | fittedValues > upper)] = replacement
// =============================================================================
inline std::vector<Real> har_insanity_filter(
        const std::vector<Real>& fitted,
        Real lower, Real upper, Real replacement) {
    std::vector<Real> result = fitted;
    for (Size i = 0; i < result.size(); ++i) {
        if (result[i] < lower || result[i] > upper) {
            result[i] = replacement;
        }
    }
    return result;
}

// =============================================================================
// OLS 估计: y = X beta + epsilon (含常数项)
//
// 构造设计矩阵 [1, X], 求解 (X^T X) beta = X^T y
// 使用 Gauss-Jordan 消元 (来自 optimizer.hpp 的 solve_linear_system 逻辑)
// =============================================================================
inline std::vector<Real> ols_estimate(
        const std::vector<Real>& y,
        const std::vector<std::vector<Real>>& X,
        std::vector<Real>& fitted_values,
        std::vector<Real>& residuals,
        Real& r_squared,
        Real& adj_r_squared,
        Real& llh) {

    const Size n = y.size();
    if (n == 0) {
        throw std::invalid_argument("ols_estimate: empty y");
    }
    const Size k = X.empty() ? 0 : X[0].size();
    const Size p = k + 1;  // 含常数项

    if (n < p) {
        throw std::invalid_argument(
            "ols_estimate: insufficient observations");
    }

    // 构造设计矩阵 [1, X] (n x p)
    std::vector<std::vector<Real>> design(n, std::vector<Real>(p));
    for (Size i = 0; i < n; ++i) {
        design[i][0] = 1.0;  // 常数项
        for (Size j = 0; j < k; ++j) {
            design[i][j + 1] = X[i][j];
        }
    }

    // 正规方程: (X^T X) beta = X^T y
    std::vector<std::vector<Real>> XtX(p, std::vector<Real>(p, 0.0));
    std::vector<Real> Xty(p, 0.0);
    for (Size i = 0; i < p; ++i) {
        for (Size j = 0; j < p; ++j) {
            Real s = 0.0;
            for (Size t = 0; t < n; ++t) {
                s += design[t][i] * design[t][j];
            }
            XtX[i][j] = s;
        }
        Real s = 0.0;
        for (Size t = 0; t < n; ++t) {
            s += design[t][i] * y[t];
        }
        Xty[i] = s;
    }

    // Gauss-Jordan 消元求解 (partial pivoting)
    std::vector<std::vector<Real>> A = XtX;
    std::vector<Real> b = Xty;
    for (Size col = 0; col < p; ++col) {
        // Partial pivot
        Size piv = col;
        Real maxv = std::fabs(A[col][col]);
        for (Size r = col + 1; r < p; ++r) {
            Real v = std::fabs(A[r][col]);
            if (v > maxv) { maxv = v; piv = r; }
        }
        if (maxv < 1e-15) {
            throw std::runtime_error(
                "ols_estimate: singular design matrix");
        }
        if (piv != col) {
            std::swap(A[piv], A[col]);
            std::swap(b[piv], b[col]);
        }
        Real akk = A[col][col];
        for (Size r = col + 1; r < p; ++r) {
            Real f = A[r][col] / akk;
            if (f == 0.0) continue;
            for (Size c = col; c < p; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    // Back substitution
    std::vector<Real> beta(p);
    for (Size i = p; i-- > 0;) {
        Real s = b[i];
        for (Size j = i + 1; j < p; ++j) s -= A[i][j] * beta[j];
        beta[i] = s / A[i][i];
    }

    // 计算 fitted values 和 residuals
    fitted_values.resize(n);
    residuals.resize(n);
    for (Size t = 0; t < n; ++t) {
        Real pred = beta[0];
        for (Size j = 0; j < k; ++j) {
            pred += beta[j + 1] * X[t][j];
        }
        fitted_values[t] = pred;
        residuals[t] = y[t] - pred;
    }

    // R^2
    Real y_mean = 0.0;
    for (Size t = 0; t < n; ++t) y_mean += y[t];
    y_mean /= static_cast<Real>(n);
    Real ss_tot = 0.0, ss_res = 0.0;
    for (Size t = 0; t < n; ++t) {
        Real dy = y[t] - y_mean;
        ss_tot += dy * dy;
        ss_res += residuals[t] * residuals[t];
    }
    r_squared = (ss_tot > 1e-300) ? (1.0 - ss_res / ss_tot) : 0.0;
    adj_r_squared = (n > p && ss_tot > 1e-300)
        ? (1.0 - (1.0 - r_squared) * static_cast<Real>(n - 1) /
                   static_cast<Real>(n - p))
        : 0.0;

    // 对数似然 (Gaussian): llh = -n/2 * (log(2pi) + log(ss_res/n) + 1)
    if (n > 0 && ss_res > 0.0) {
        llh = -0.5 * static_cast<Real>(n) *
              (std::log(2.0 * PI) +
               std::log(ss_res / static_cast<Real>(n)) + 1.0);
    } else {
        llh = -0.5 * static_cast<Real>(n) * std::log(2.0 * PI);
    }

    return beta;
}

// =============================================================================
// HARmodel 估计器
// =============================================================================
class HarModel {
public:
    // =====================================================================
    // 主接口: HAR-RV 模型估计
    //
    // 输入:
    //   RM1       - 已实现测度 (realized variance), 长度 n
    //   RM2       - 跳跃稳健测度 (realized bipower variation), type=HARJ/HARCJ/CHAR 时需要
    //   RM3       - 已实现四阶矩 (realized quarticity), type=HARQ/HARQJ/CHARQ 时需要
    //   periods   - RV 聚合周期 (默认 {1, 5, 22})
    //   periodsJ  - 跳跃聚合周期 (默认 {1, 5, 22})
    //   periodsQ  - 四阶矩聚合周期 (默认 {1})
    //   type      - 模型类型
    //   h         - 预测步长 (默认 1)
    //   transform - 变换类型 (默认 None)
    //   alpha     - 跳跃检验显著性水平 (默认 0.05, type=HARCJ 时使用)
    //
    // 返回: HarModelResult
    // 异常: 数据不足 / 奇异矩阵 抛 invalid_argument/runtime_error
    // =====================================================================
    static HarModelResult estimate(
            const std::vector<Real>& RM1,
            const std::vector<Real>& RM2 = {},
            const std::vector<Real>& RM3 = {},
            const std::vector<Size>& periods = {1, 5, 22},
            const std::vector<Size>& periodsJ = {1, 5, 22},
            const std::vector<Size>& periodsQ = {1},
            HarType type = HarType::HAR,
            Size h = 1,
            HarTransform transform = HarTransform::None,
            Real alpha = 0.05) {

        const Size n = RM1.size();
        if (n < 2) {
            throw std::invalid_argument(
                "HarModel: need at least 2 observations");
        }

        // 验证 type 与输入数据的一致性
        const bool is_jump_model = (type == HarType::HARJ ||
                                     type == HarType::HARCJ ||
                                     type == HarType::HARQJ);
        const bool is_quarticity_model = (type == HarType::HARQ ||
                                           type == HarType::HARQJ ||
                                           type == HarType::CHARQ);
        const bool is_bpv_model = (type == HarType::CHAR ||
                                    type == HarType::CHARQ);

        if (is_jump_model && RM2.empty()) {
            throw std::invalid_argument(
                "HarModel: jump models require RM2 (bipower variation)");
        }
        if (is_quarticity_model && RM3.empty()) {
            throw std::invalid_argument(
                "HarModel: quarticity models require RM3 (realized quarticity)");
        }
        if (is_bpv_model && RM2.empty()) {
            throw std::invalid_argument(
                "HarModel: CHAR models require RM2 (bipower variation)");
        }

        // maxp: 最大聚合阶数
        Size maxp = 0;
        for (Size p : periods) maxp = std::max(maxp, p);
        if (is_jump_model) {
            for (Size p : periodsJ) maxp = std::max(maxp, p);
        }
        if (is_quarticity_model) {
            for (Size p : periodsQ) maxp = std::max(maxp, p);
        }
        maxp = std::max(maxp, h);

        // 检查数据长度
        if (n < maxp + h + 1) {
            throw std::invalid_argument(
                "HarModel: insufficient data for given periods and h");
        }

        // 步骤 1: 聚合 RM1
        std::vector<std::vector<Real>> RVmatrix1 = har_agg(RM1, periods);

        // 步骤 2: 构造 y = har_agg(RM1, h, 1)[(maxp+h):n]
        // R: y <- har_agg(RM1, c(h), 1L)[(maxp+h):(n)]
        // 注意 R 索引 1-based, (maxp+h):(n) 对应 C++ 0-based [maxp+h-1 .. n-1]
        std::vector<Real> y_h = har_agg_single(RM1, h);
        std::vector<Real> y;
        for (Size i = maxp + h - 1; i < n; ++i) {
            y.push_back(y_h[i]);
        }

        // 步骤 3: 构造 x1 = RVmatrix1[(maxp:(n-h)), ]
        // R: x1 <- RVmatrix1[(maxp:(n-h)), ]
        // R 1-based: maxp:(n-h) 对应 C++ 0-based [maxp-1 .. n-h-1]
        const Size n_est = n - h - maxp + 1;  // 估计样本大小
        if (n_est != y.size()) {
            throw std::runtime_error(
                "HarModel: dimension mismatch in y and x");
        }

        std::vector<std::vector<Real>> x1(n_est, std::vector<Real>(periods.size()));
        for (Size i = 0; i < n_est; ++i) {
            Size row = maxp - 1 + i;  // C++ 0-based
            for (Size j = 0; j < periods.size(); ++j) {
                x1[i][j] = RVmatrix1[row][j];
            }
        }

        // 步骤 4: 根据 type 构造完整设计矩阵
        std::vector<std::vector<Real>> x = x1;
        std::vector<std::string> x_names;
        for (Size p : periods) x_names.push_back("RV" + std::to_string(p));

        if (is_jump_model) {
            // J = pmax(RM1 - RM2, 0)
            std::vector<Real> J(RM1.size());
            for (Size i = 0; i < RM1.size(); ++i) {
                J[i] = std::max(RM1[i] - RM2[i], 0.0);
            }
            auto Jmatrix = har_agg(J, periodsJ);
            // J <- J[(maxp:(n-h)), ]
            for (Size i = 0; i < n_est; ++i) {
                Size row = maxp - 1 + i;
                std::vector<Real> j_row(periodsJ.size());
                for (Size j = 0; j < periodsJ.size(); ++j) {
                    j_row[j] = Jmatrix[row][j];
                }
                // 追加到 x[i]
                x[i].insert(x[i].end(), j_row.begin(), j_row.end());
            }
            for (Size p : periodsJ) x_names.push_back("J" + std::to_string(p));
        }

        if (is_quarticity_model) {
            // RQmatrix = har_agg(RM3, periodsQ)
            // BPQ 变换: sqrt(RQmatrix) - sqrt(mean(RM3))
            Real mean_RM3 = 0.0;
            for (Size i = 0; i < RM3.size(); ++i) mean_RM3 += RM3[i];
            mean_RM3 /= static_cast<Real>(RM3.size());
            const Real sqrt_mean_RM3 = std::sqrt(std::max(mean_RM3, 0.0));

            auto RQmatrix = har_agg(RM3, periodsQ);
            // RQmatrix <- sqrt(RQmatrix) - sqrt(mean(RM3))
            // 然后乘以 x1 的对应列 (RQ * RV)
            for (Size i = 0; i < n_est; ++i) {
                Size row = maxp - 1 + i;
                std::vector<Real> rq_row(periodsQ.size());
                for (Size j = 0; j < periodsQ.size(); ++j) {
                    Real rq = RQmatrix[row][j];
                    if (std::isnan(rq)) rq = 0.0;
                    rq = std::sqrt(std::max(rq, 0.0)) - sqrt_mean_RM3;
                    // RQ * RV (对应列相乘)
                    rq_row[j] = rq * x1[i][j];
                }
                x[i].insert(x[i].end(), rq_row.begin(), rq_row.end());
            }
            for (Size p : periodsQ) x_names.push_back("RQ" + std::to_string(p));
        }

        if (is_bpv_model) {
            // CHAR: 用 RM2 (BPV) 替代 RM1 作为回归量
            // x2 = har_agg(RM2, periods)
            auto RVmatrix2 = har_agg(RM2, periods);
            // x2 <- RVmatrix2[(maxp:(n-h)), ]
            // CHAR 模型用 x2 替代 x1 (不包含 RV, 只包含 BPV)
            x.clear();
            x_names.clear();
            for (Size i = 0; i < n_est; ++i) {
                Size row = maxp - 1 + i;
                std::vector<Real> row_data(periods.size());
                for (Size j = 0; j < periods.size(); ++j) {
                    row_data[j] = RVmatrix2[row][j];
                }
                x.push_back(row_data);
            }
            for (Size p : periods) x_names.push_back("RV" + std::to_string(p));

            // CHARQ: 额外添加 RQ * x2
            if (is_quarticity_model) {
                Real mean_RM3 = 0.0;
                for (Size i = 0; i < RM3.size(); ++i) mean_RM3 += RM3[i];
                mean_RM3 /= static_cast<Real>(RM3.size());
                const Real sqrt_mean_RM3 = std::sqrt(std::max(mean_RM3, 0.0));

                auto RQmatrix = har_agg(RM3, periodsQ);
                for (Size i = 0; i < n_est; ++i) {
                    Size row = maxp - 1 + i;
                    for (Size j = 0; j < periodsQ.size(); ++j) {
                        Real rq = RQmatrix[row][j];
                        if (std::isnan(rq)) rq = 0.0;
                        rq = std::sqrt(std::max(rq, 0.0)) - sqrt_mean_RM3;
                        x[i].push_back(rq * x[i][j]);
                    }
                }
                for (Size p : periodsQ) x_names.push_back("RQ" + std::to_string(p));
            }
        }

        // 步骤 5: 应用变换
        std::vector<Real> y_transformed = y;
        std::vector<std::vector<Real>> x_transformed = x;

        if (transform == HarTransform::Log) {
            for (Size i = 0; i < y_transformed.size(); ++i) {
                y_transformed[i] = std::log(std::max(y_transformed[i], 1e-300));
            }
            for (Size i = 0; i < x_transformed.size(); ++i) {
                for (Size j = 0; j < x_transformed[i].size(); ++j) {
                    x_transformed[i][j] = std::log(
                        std::max(x_transformed[i][j], 1e-300));
                }
            }
        } else if (transform == HarTransform::Sqrt) {
            for (Size i = 0; i < y_transformed.size(); ++i) {
                y_transformed[i] = std::sqrt(std::max(y_transformed[i], 0.0));
            }
            for (Size i = 0; i < x_transformed.size(); ++i) {
                for (Size j = 0; j < x_transformed[i].size(); ++j) {
                    x_transformed[i][j] = std::sqrt(
                        std::max(x_transformed[i][j], 0.0));
                }
            }
        }

        // 步骤 6: OLS 估计
        std::vector<Real> fitted, residuals;
        Real r_sq, adj_r_sq, llh;
        std::vector<Real> beta = ols_estimate(
            y_transformed, x_transformed,
            fitted, residuals, r_sq, adj_r_sq, llh);

        // 步骤 7: harInsanityFilter
        // R: lower = min(y, 0), upper = max(y), replacement = mean(y)
        Real y_min = y_transformed[0], y_max = y_transformed[0], y_mean = 0.0;
        for (Size i = 0; i < y_transformed.size(); ++i) {
            y_min = std::min(y_min, y_transformed[i]);
            y_max = std::max(y_max, y_transformed[i]);
            y_mean += y_transformed[i];
        }
        y_mean /= static_cast<Real>(y_transformed.size());
        Real lower = std::min(y_min, 0.0);
        Real upper = y_max;
        fitted = har_insanity_filter(fitted, lower, upper, y_mean);

        // 重新计算残差 (使用 filter 后的 fitted)
        for (Size i = 0; i < residuals.size(); ++i) {
            residuals[i] = y_transformed[i] - fitted[i];
        }

        // 步骤 8: 构造结果
        HarModelResult result;
        result.coefficients = beta;
        result.coef_names.push_back("beta0");
        for (const auto& name : x_names) {
            result.coef_names.push_back(name);
        }
        result.fitted_values = fitted;
        result.residuals = residuals;
        result.r_squared = r_sq;
        result.adj_r_squared = adj_r_sq;
        result.llh = llh;
        result.type = type;
        result.transform = transform;
        result.h = h;
        result.n_obs = n_est;
        result.maxp = maxp;
        result.y = y_transformed;
        result.x = x_transformed;

        return result;
    }

    // =====================================================================
    // 便捷接口: 基础 HAR-RV (type="HAR")
    // =====================================================================
    static HarModelResult estimate_har(
            const std::vector<Real>& RM1,
            const std::vector<Size>& periods = {1, 5, 22},
            Size h = 1,
            HarTransform transform = HarTransform::None) {
        return estimate(RM1, {}, {}, periods, {}, {}, HarType::HAR,
                        h, transform, 0.05);
    }

    // =====================================================================
    // 便捷接口: HAR-RV-J (type="HARJ")
    // =====================================================================
    static HarModelResult estimate_harj(
            const std::vector<Real>& RM1,
            const std::vector<Real>& RM2,
            const std::vector<Size>& periods = {1, 5, 22},
            const std::vector<Size>& periodsJ = {1, 5, 22},
            Size h = 1,
            HarTransform transform = HarTransform::None) {
        return estimate(RM1, RM2, {}, periods, periodsJ, {},
                        HarType::HARJ, h, transform, 0.05);
    }

    // =====================================================================
    // 便捷接口: HAR-Q (type="HARQ")
    // =====================================================================
    static HarModelResult estimate_harq(
            const std::vector<Real>& RM1,
            const std::vector<Real>& RM3,
            const std::vector<Size>& periods = {1, 5, 22},
            const std::vector<Size>& periodsQ = {1},
            Size h = 1,
            HarTransform transform = HarTransform::None) {
        return estimate(RM1, {}, RM3, periods, {}, periodsQ,
                        HarType::HARQ, h, transform, 0.05);
    }

    // =====================================================================
    // 便捷接口: CHAR (type="CHAR")
    // =====================================================================
    static HarModelResult estimate_char(
            const std::vector<Real>& RM1,
            const std::vector<Real>& RM2,
            const std::vector<Size>& periods = {1, 5, 22},
            Size h = 1,
            HarTransform transform = HarTransform::None) {
        return estimate(RM1, RM2, {}, periods, {}, {},
                        HarType::CHAR, h, transform, 0.05);
    }

    // =====================================================================
    // 预测: 使用最后一行数据预测未来 h 步
    //
    // R: predict.HARmodel (L576-925)
    //   newdata = NULL 时, 用最后观测值预测
    //   y_hat = beta0 + sum(beta_i * x_last_i)
    //   transform="log" 时, backtransform="simple" → exp(y_hat)
    //   transform="sqrt" 时, backtransform="simple" → y_hat^2
    // =====================================================================
    static Real predict_last(const HarModelResult& model) {
        if (model.x.empty() || model.x.back().empty()) {
            return model.coefficients[0];
        }
        Real pred = model.coefficients[0];
        for (Size j = 0; j < model.x.back().size(); ++j) {
            pred += model.coefficients[j + 1] * model.x.back()[j];
        }

        // Backtransform
        if (model.transform == HarTransform::Log) {
            return std::exp(pred);
        } else if (model.transform == HarTransform::Sqrt) {
            return pred * pred;
        }
        return pred;
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
