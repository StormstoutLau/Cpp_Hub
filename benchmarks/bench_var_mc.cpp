// Phase 4 LITE - E2: 10000 场景 MC VaR Full Revaluation 基准
// 补齐 Phase 3 审计整改项 E2: "VaR 计算, 10000 场景 < 500ms"
//
// 基准设计:
//   - 构造一个 5 因子组合 (5 个 BSM 欧式期权的多空组合)
//   - 使用 MCVaR 引擎, Full Revaluation 模式, 10000 条路径
//   - 测量总耗时 (含路径生成 + 组合估值 + 分位数计算)
//   - 同时测量 Delta-Gamma 和 Delta 近似作为对照
//
// 验收门槛: Full Revaluation < 500ms (Phase 3 审计 E2)
//
// 输出格式:
//   [BENCH] E2 MC VaR (paths=10000, factors=5)
//   [BENCH] Full:        total=XXX ms, per_path=XX us, throughput=XXX paths/s
//   [BENCH] DeltaGamma:  total=XXX ms, ...
//   [BENCH] Delta:       total=XXX ms, ...
//   [PASS/FAIL] Full threshold 500ms: ...

#include "cpphub/risk/var/mc_var.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include <chrono>
#include <cstdio>
#include <vector>

using namespace cpphub::v1;

// 5 因子组合: 5 个 BSM 欧式期权,风险因子为 5 个标的资产价格
// 组合价值 = sum(权重 * 期权价格)
struct Portfolio {
    std::vector<Real> weights;     // 每个期权的数量 (正=多头, 负=空头)
    std::vector<Real> K;           // strike
    std::vector<Real> T;           // maturity
    std::vector<Real> r_vec;       // rate
    std::vector<Real> q_vec;       // dividend
    std::vector<Real> sigma_vec;   // vol
    std::vector<bool> is_call;
    Real base_value;

    Real value(const std::vector<Real>& spots) const {
        Real total = 0.0;
        for (Size i = 0; i < spots.size(); ++i) {
            Real price = is_call[i]
                ? bsm_call_price(spots[i], K[i], T[i], r_vec[i], q_vec[i], sigma_vec[i])
                : bsm_put_price(spots[i], K[i], T[i], r_vec[i], q_vec[i], sigma_vec[i]);
            total += weights[i] * price;
        }
        return total;
    }
};

int main() {
    constexpr Size N_PATHS = 10000;
    constexpr Size N_FACTORS = 5;
    constexpr Real THRESHOLD_FULL_MS = 500.0;

    // 构造组合: 5 个期权,当前 spot=100
    Portfolio pf;
    pf.weights  = {+100.0, -50.0, +75.0, -25.0, +50.0};
    pf.K        = {95.0, 105.0, 100.0, 110.0, 90.0};
    pf.T        = {0.25, 0.5, 1.0, 0.75, 0.3};
    pf.r_vec    = {0.03, 0.03, 0.03, 0.03, 0.03};
    pf.q_vec    = {0.01, 0.01, 0.01, 0.01, 0.01};
    pf.sigma_vec= {0.2, 0.25, 0.18, 0.3, 0.22};
    pf.is_call  = {true, false, true, false, true};

    std::vector<Real> current_spots(N_FACTORS, 100.0);
    pf.base_value = pf.value(current_spots);

    // 协方差矩阵: 5 个标的, 10% 波动率, 中等相关
    // 对角 0.1^2=0.01, 非对角 0.3*0.1*0.1=0.003
    std::vector<Real> cov(N_FACTORS * N_FACTORS, 0.0);
    for (Size i = 0; i < N_FACTORS; ++i) {
        cov[i * N_FACTORS + i] = 0.01;  // variance = 0.1^2
        for (Size j = i + 1; j < N_FACTORS; ++j) {
            cov[i * N_FACTORS + j] = 0.003;  // covariance
            cov[j * N_FACTORS + i] = 0.003;
        }
    }

    // MCVaR 配置
    MCVarConfig cfg;
    cfg.n_paths = N_PATHS;
    cfg.seed = 42;
    cfg.antithetic = true;
    cfg.use_control_variate = false;

    MCVaR mcvar(
        [&pf](const std::vector<Real>& spots) { return pf.value(spots); },
        current_spots, cov, N_FACTORS, cfg);

    // ===== Full Revaluation 基准 =====
    auto t0 = std::chrono::steady_clock::now();
    Real var_full = mcvar.var(0.99, VaRApproximation::Full);
    auto t1 = std::chrono::steady_clock::now();
    double full_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double full_per_us = full_ms * 1000.0 / static_cast<double>(N_PATHS);
    double full_throughput = static_cast<double>(N_PATHS) / (full_ms / 1000.0);

    // ===== Delta-Gamma 近似基准 (对照) =====
    // 计算 delta/gamma 向量 (每个期权对各自 spot 的二阶导)
    std::vector<Real> delta_vec(N_FACTORS), gamma_vec(N_FACTORS * N_FACTORS, 0.0);
    for (Size i = 0; i < N_FACTORS; ++i) {
        Real h = 1e-4 * current_spots[i];
        std::vector<Real> sp_plus = current_spots, sp_minus = current_spots;
        sp_plus[i] += h; sp_minus[i] -= h;
        Real v0 = pf.value(sp_minus);
        Real v1 = pf.value(current_spots);
        Real v2 = pf.value(sp_plus);
        delta_vec[i] = (v2 - v0) / (2.0 * h);
        gamma_vec[i * N_FACTORS + i] = (v2 - 2.0 * v1 + v0) / (h * h);
    }

    auto t2 = std::chrono::steady_clock::now();
    Real var_dg = mcvar.var(0.99, VaRApproximation::DeltaGamma, delta_vec, gamma_vec);
    auto t3 = std::chrono::steady_clock::now();
    double dg_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();
    double dg_per_us = dg_ms * 1000.0 / static_cast<double>(N_PATHS);
    double dg_throughput = static_cast<double>(N_PATHS) / (dg_ms / 1000.0);

    // ===== Delta 近似基准 (对照) =====
    auto t4 = std::chrono::steady_clock::now();
    Real var_d = mcvar.var(0.99, VaRApproximation::Delta, delta_vec, {});
    auto t5 = std::chrono::steady_clock::now();
    double d_ms = std::chrono::duration<double, std::milli>(t5 - t4).count();
    double d_per_us = d_ms * 1000.0 / static_cast<double>(N_PATHS);
    double d_throughput = static_cast<double>(N_PATHS) / (d_ms / 1000.0);

    // ===== 输出 =====
    std::printf("[BENCH] E2 MC VaR (paths=%zu, factors=%zu)\n",
                static_cast<size_t>(N_PATHS), static_cast<size_t>(N_FACTORS));
    std::printf("[BENCH] Full:       total=%7.2f ms, per_path=%7.2f us, throughput=%10.0f paths/s, VaR(99%%)=%.4f\n",
                full_ms, full_per_us, full_throughput, var_full);
    std::printf("[BENCH] DeltaGamma: total=%7.2f ms, per_path=%7.2f us, throughput=%10.0f paths/s, VaR(99%%)=%.4f\n",
                dg_ms, dg_per_us, dg_throughput, var_dg);
    std::printf("[BENCH] Delta:      total=%7.2f ms, per_path=%7.2f us, throughput=%10.0f paths/s, VaR(99%%)=%.4f\n",
                d_ms, d_per_us, d_throughput, var_d);

    bool full_pass = full_ms < THRESHOLD_FULL_MS;
    std::printf("[PASS] Full Revaluation threshold %.0fms: %s (%.2f ms)\n",
                THRESHOLD_FULL_MS, full_pass ? "PASS" : "FAIL", full_ms);

    // 正确性检查: Full 与 Delta-Gamma 的 VaR 应在同一数量级 (相对差 < 30%)
    Real rel_diff = std::abs(var_full - var_dg) / std::max(std::abs(var_full), Real(1e-10));
    bool sanity_pass = rel_diff < 0.3;
    std::printf("[PASS] Full vs DeltaGamma sanity (rel_diff<0.3): %s (%.3f)\n",
                sanity_pass ? "PASS" : "FAIL", rel_diff);

    int exit_code = 0;
    if (!full_pass) exit_code = 1;
    if (!sanity_pass) exit_code = 2;
    return exit_code;
}
