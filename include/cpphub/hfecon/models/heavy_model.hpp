// =============================================================================
// heavy_model.hpp
// Phase 5 v1.4.2 Wave C - High frEquency bAsed VolatilitY (HEAVY) Model
//
// R 对照: HEAVYmodel(data, startingValues=NULL)
//
// 文献: Shephard & Sheppard (2010), Restat 92(2), doi:10.1162/REST_a_00017
//       Noureldin, Shephard, Sheppard (2012), JAE 27(8), doi:10.1002/jae.1260
//
// SOURCE: PHASE5_HFE_SPEC §5.2, §5.3 D8, §5.4
//   R highfrequency 1.0.3 src/HEAVYmodel.cpp calcRecVarEq (L5-15)
//   R highfrequency 1.0.3 R/HEAVYmodel.R HEAVYmodel (L72-137)
//   R highfrequency 1.0.3 R/internalHEAVY.R heavyLLH (L32-39)
//
// 关键幻觉排除 (spec §5.3 D8):
//   强制去均值 (data - mean(data)), 无常数项
//   用数值优化 (R 用 solnp, C++ 用 NelderMead) 最大化 MLE
//
// 模型 (Shephard & Sheppard 2010):
//   r_t = h_t^{1/2} * Z_t,  Z_t ~ iid(0, 1)
//
//   方差方程 (Variance Equation):
//     h_t = omega + alpha * RM_{t-1} + beta * h_{t-1}
//     约束: omega > 0, alpha > 0, beta > 0, alpha + beta < 1
//
//   RM 方程 (Measurement Equation):
//     mu_t = omega_R + alpha_R * RM_{t-1} + beta_R * mu_{t-1}
//     约束: omega_R > 0, alpha_R > 0, beta_R > 0,
//           alpha_R < 1, beta_R < 1, alpha_R + beta_R < 1
//
// 算法 (R 源码实测):
//   1. ret <- data[,1]; rm <- data[,2]
//   2. ret <- ret - mean(ret)  # 强制去均值
//   3. 起始值:
//      Var: [mean(ret^2)*(1 - 0.3*mean(rm)/mean(ret^2) - 0.5), 0.3, 0.5]
//      RM:  [mean(rm)*(1 - 0.6 - 0.3), 0.6, 0.3]
//   4. MLE 估计 (分别估计两个方程):
//      Var: max sum(heavyLLH(par, ret=ret, rm=rm))
//      RM:  max sum(heavyLLH(par, rm=rm, RMEq=TRUE))
//   5. calcRecVarEq: 方差递归
//      g[0] = mean(rm)  # 注意: g[0] = mean(rm), 不是 mean(ret^2)
//      g[i] = par[0] + par[1] * rm[i-1] + par[2] * g[i-1]
//   6. heavyLLH:
//      RMEq=TRUE:  -1/2 * log(2pi) - 1/2 * (log(g) + rm/g)
//      RMEq=FALSE: -1/2 * log(2pi) - 1/2 * (log(g) + ret^2/g)
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/calibration/optimizer.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// HEAVY 模型估计结果
// =============================================================================
struct HeavyModelResult {
    // 方差方程系数: [omega, alpha, beta]
    std::vector<Real> var_coefficients;
    // RM 方程系数: [omega_R, alpha_R, beta_R]
    std::vector<Real> rm_coefficients;

    // 全部系数 (6 个): [omega, alpha, beta, omega_R, alpha_R, beta_R]
    std::vector<Real> coefficients;
    std::vector<std::string> coef_names;

    // 对数似然: [llhVar, llhRM]
    Real llh_var;
    Real llh_rm;

    // 条件方差序列
    std::vector<Real> var_cond_variances;  // h_t (方差方程)
    std::vector<Real> rm_cond_variances;   // mu_t (RM 方程)

    // 残差 (标准化)
    std::vector<Real> residuals;  // ret / sqrt(h_t)

    // 原始数据
    std::vector<Real> ret;  // 去均值后的收益率
    std::vector<Real> rm;   // 已实现测度

    Size n_obs;
    bool converged;
};

// =============================================================================
// calcRecVarEq: 方差递归 (R calcRecVarEq 对标)
//
// R 源码 (HEAVYmodel.cpp L5-15):
//   g[0] = mean(rm)
//   for (i = 1..n-1):
//     g[i] = par[0] + par[1] * rm[i-1] + par[2] * g[i-1]
//
// 输入: par (3 个参数), rm (n 个观测)
// 输出: g (n 个条件方差)
// =============================================================================
inline std::vector<Real> calc_rec_var_eq(
        const std::vector<Real>& par,
        const std::vector<Real>& rm) {
    const Size n = rm.size();
    if (n == 0) return {};
    if (par.size() < 3) {
        throw std::invalid_argument(
            "calc_rec_var_eq: need 3 parameters");
    }
    std::vector<Real> g(n);
    // g[0] = mean(rm)
    Real mean_rm = 0.0;
    for (Size i = 0; i < n; ++i) mean_rm += rm[i];
    mean_rm /= static_cast<Real>(n);
    g[0] = mean_rm;

    for (Size i = 1; i < n; ++i) {
        g[i] = par[0] + par[1] * rm[i - 1] + par[2] * g[i - 1];
        // 保证 g[i] > 0 (数值稳定性)
        if (g[i] <= 0.0) g[i] = 1e-10;
    }
    return g;
}

// =============================================================================
// heavy_llh: 对数似然 (R heavyLLH 对标)
//
// R 源码 (internalHEAVY.R L32-39):
//   RMEq=TRUE:  -1/2 * log(2pi) - 1/2 * (log(condVar) + rm / condVar)
//   RMEq=FALSE: -1/2 * log(2pi) - 1/2 * (log(condVar) + ret^2 / condVar)
//
// 输入: par (3 个参数), ret (收益率, RMEq=FALSE 时用), rm (RM),
//       rm_eq (TRUE=RM方程, FALSE=方差方程)
// 输出: 对数似然向量 (长度 n)
// =============================================================================
inline std::vector<Real> heavy_llh(
        const std::vector<Real>& par,
        const std::vector<Real>& ret,
        const std::vector<Real>& rm,
        bool rm_eq) {
    std::vector<Real> cond_var;
    if (rm_eq) {
        cond_var = calc_rec_var_eq(par, rm);
    } else {
        // 方差方程: condVar = calcRecVarEq(par, ret^2)
        std::vector<Real> ret_sq(ret.size());
        for (Size i = 0; i < ret.size(); ++i) ret_sq[i] = ret[i] * ret[i];
        cond_var = calc_rec_var_eq(par, ret_sq);
    }

    const Size n = rm_eq ? rm.size() : ret.size();
    std::vector<Real> llh(n);
    const Real log_2pi = std::log(2.0 * PI);
    for (Size i = 0; i < n; ++i) {
        const Real g = std::max(cond_var[i], 1e-10);
        const Real obs = rm_eq ? rm[i] : ret[i] * ret[i];
        llh[i] = -0.5 * log_2pi - 0.5 * (std::log(g) + obs / g);
    }
    return llh;
}

// =============================================================================
// HEAVY 模型估计器
// =============================================================================
class HeavyModel {
public:
    // =====================================================================
    // 主接口: HEAVY 模型估计
    //
    // 输入:
    //   ret            - 收益率序列 (长度 n)
    //   rm             - 已实现测度序列 (长度 n)
    //   starting_values - 起始值 (6 个: [omega, alpha, beta, omegaR, alphaR, betaR])
    //                     空向量表示自动计算
    //
    // 返回: HeavyModelResult
    // 异常: 数据不足 / 优化失败 抛 invalid_argument/runtime_error
    // =====================================================================
    static HeavyModelResult estimate(
            const std::vector<Real>& ret_in,
            const std::vector<Real>& rm_in,
            const std::vector<Real>& starting_values = {}) {

        const Size n = ret_in.size();
        if (n < 5) {
            throw std::invalid_argument(
                "HeavyModel: need at least 5 observations");
        }
        if (rm_in.size() != n) {
            throw std::invalid_argument(
                "HeavyModel: ret and rm must have equal length");
        }

        // 步骤 1: 强制去均值 (D8)
        std::vector<Real> ret = ret_in;
        Real mean_ret = 0.0;
        for (Size i = 0; i < n; ++i) mean_ret += ret[i];
        mean_ret /= static_cast<Real>(n);
        for (Size i = 0; i < n; ++i) ret[i] -= mean_ret;

        std::vector<Real> rm = rm_in;

        // 步骤 2: 起始值
        std::vector<Real> start_var(3), start_rm(3);
        if (starting_values.size() >= 6) {
            start_var[0] = starting_values[0];
            start_var[1] = starting_values[1];
            start_var[2] = starting_values[2];
            start_rm[0] = starting_values[3];
            start_rm[1] = starting_values[4];
            start_rm[2] = starting_values[5];
        } else {
            // R 默认起始值
            Real mean_ret_sq = 0.0;
            for (Size i = 0; i < n; ++i) mean_ret_sq += ret[i] * ret[i];
            mean_ret_sq /= static_cast<Real>(n);

            Real mean_rm = 0.0;
            for (Size i = 0; i < n; ++i) mean_rm += rm[i];
            mean_rm /= static_cast<Real>(n);

            // Var: [mean(ret^2)*(1 - 0.3*mean(rm)/mean(ret^2) - 0.5), 0.3, 0.5]
            Real ratio = (mean_ret_sq > 1e-300) ? (mean_rm / mean_ret_sq) : 0.0;
            start_var[0] = mean_ret_sq * (1.0 - 0.3 * ratio - 0.5);
            start_var[1] = 0.3;
            start_var[2] = 0.5;

            // RM: [mean(rm)*(1 - 0.6 - 0.3), 0.6, 0.3]
            start_rm[0] = mean_rm * (1.0 - 0.6 - 0.3);
            start_rm[1] = 0.6;
            start_rm[2] = 0.3;
        }

        // 步骤 3: MLE 估计 (NelderMead + penalty)
        // 方差方程: max sum(heavyLLH(par, ret=ret, rm=rm, RMEq=FALSE))
        auto neg_llh_var = [&](const std::vector<Real>& par) -> Real {
            // 约束 penalty: omega>0, alpha>0, beta>0, alpha+beta<1
            if (par[0] <= 0.0 || par[1] <= 0.0 || par[2] <= 0.0) {
                return 1e10;
            }
            if (par[1] + par[2] >= 1.0) return 1e10;
            auto llh = heavy_llh(par, ret, rm, false);
            Real s = 0.0;
            for (Real v : llh) s += v;
            return -s;  // 最小化负对数似然
        };

        // RM 方程: max sum(heavyLLH(par, rm=rm, RMEq=TRUE))
        auto neg_llh_rm = [&](const std::vector<Real>& par) -> Real {
            // 约束: omega_R>0, alpha_R>0, beta_R>0, alpha_R<1, beta_R<1, alpha_R+beta_R<1
            if (par[0] <= 0.0 || par[1] <= 0.0 || par[2] <= 0.0) {
                return 1e10;
            }
            if (par[1] >= 1.0 || par[2] >= 1.0) return 1e10;
            if (par[1] + par[2] >= 1.0) return 1e10;
            auto llh = heavy_llh(par, ret, rm, true);
            Real s = 0.0;
            for (Real v : llh) s += v;
            return -s;
        };

        NelderMead::Config cfg;
        cfg.max_iterations = 2000;
        cfg.ftol = 1e-10;
        cfg.xtol = 1e-10;

        OptimizationResult res_var = NelderMead::minimize(
            neg_llh_var, start_var, cfg);
        OptimizationResult res_rm = NelderMead::minimize(
            neg_llh_rm, start_rm, cfg);

        // 步骤 4: 计算条件方差
        std::vector<Real> var_cond = calc_rec_var_eq(res_var.x, rm);
        std::vector<Real> rm_cond = calc_rec_var_eq(res_rm.x, rm);

        // 步骤 5: 残差 (标准化)
        std::vector<Real> residuals(n);
        for (Size i = 0; i < n; ++i) {
            Real h = std::max(var_cond[i], 1e-10);
            residuals[i] = ret[i] / std::sqrt(h);
        }

        // 步骤 6: 构造结果
        HeavyModelResult result;
        result.var_coefficients = res_var.x;
        result.rm_coefficients = res_rm.x;
        result.coefficients.resize(6);
        for (Size i = 0; i < 3; ++i) {
            result.coefficients[i] = res_var.x[i];
            result.coefficients[i + 3] = res_rm.x[i];
        }
        result.coef_names = {"omega", "alpha", "beta",
                             "omegaR", "alphaR", "betaR"};
        result.llh_var = -res_var.fx;
        result.llh_rm = -res_rm.fx;
        result.var_cond_variances = var_cond;
        result.rm_cond_variances = rm_cond;
        result.residuals = residuals;
        result.ret = ret;
        result.rm = rm;
        result.n_obs = n;
        result.converged = res_var.converged && res_rm.converged;

        return result;
    }

    // =====================================================================
    // 预测: 迭代多步预测
    //
    // R: predict.HEAVYmodel (L209-221)
    //   stepsAhead 步迭代预测
    //   oneStep = [omega + alpha*last(rm) + beta*last(h),
    //              omegaR + alphaR*last(rm) + betaR*last(mu)]
    //   k 步: 使用 Bj 矩阵递推
    //
    // 简化实现: 1 步预测
    //   h_{T+1} = omega + alpha * rm_T + beta * h_T
    //   mu_{T+1} = omega_R + alpha_R * rm_T + beta_R * mu_T
    // =====================================================================
    static std::pair<Real, Real> predict_one_step(
            const HeavyModelResult& model) {
        if (model.var_cond_variances.empty() ||
            model.rm_cond_variances.empty() ||
            model.rm.empty()) {
            throw std::runtime_error(
                "HeavyModel::predict_one_step: empty model");
        }
        const Size n = model.rm.size();
        const Real last_rm = model.rm[n - 1];
        const Real last_h = model.var_cond_variances[n - 1];
        const Real last_mu = model.rm_cond_variances[n - 1];

        const Real omega = model.var_coefficients[0];
        const Real alpha = model.var_coefficients[1];
        const Real beta = model.var_coefficients[2];
        const Real omega_r = model.rm_coefficients[0];
        const Real alpha_r = model.rm_coefficients[1];
        const Real beta_r = model.rm_coefficients[2];

        Real h_next = omega + alpha * last_rm + beta * last_h;
        Real mu_next = omega_r + alpha_r * last_rm + beta_r * last_mu;

        return {h_next, mu_next};
    }

    // =====================================================================
    // 多步迭代预测 (简化版)
    //
    // 对于 k 步预测, 需要迭代:
    //   h_{T+k} = omega + alpha * mu_{T+k-1} + beta * h_{T+k-1}
    //   mu_{T+k} = omega_R + alpha_R * mu_{T+k-1} + beta_R * mu_{T+k-1}
    // =====================================================================
    static std::vector<std::pair<Real, Real>> predict_multi_step(
            const HeavyModelResult& model, Size steps_ahead) {
        if (model.var_cond_variances.empty()) {
            throw std::runtime_error(
                "HeavyModel::predict_multi_step: empty model");
        }

        const Real omega = model.var_coefficients[0];
        const Real alpha = model.var_coefficients[1];
        const Real beta = model.var_coefficients[2];
        const Real omega_r = model.rm_coefficients[0];
        const Real alpha_r = model.rm_coefficients[1];
        const Real beta_r = model.rm_coefficients[2];

        const Size n = model.rm.size();
        Real h = model.var_cond_variances[n - 1];
        Real mu = model.rm_cond_variances[n - 1];
        Real last_rm = model.rm[n - 1];

        std::vector<std::pair<Real, Real>> forecasts;
        for (Size k = 0; k < steps_ahead; ++k) {
            Real h_next = omega + alpha * last_rm + beta * h;
            Real mu_next = omega_r + alpha_r * last_rm + beta_r * mu;
            forecasts.push_back({h_next, mu_next});
            h = h_next;
            mu = mu_next;
            last_rm = mu;  // 迭代: 下一轮的 rm 用 mu 替代
        }
        return forecasts;
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
