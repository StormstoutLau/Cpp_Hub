// Phase 4 LITE - E1: 万期权 AAD Greeks 批量基准
// 补齐 Phase 3 审计整改项 E1: "万期权全 Greeks, AAD 单次扫描 < 100ms"
//
// 基准设计:
//   - 生成 10000 个随机 BSM 欧式期权 (S, K, T, r, q, sigma 在合理范围内变化)
//   - 对每个期权调用 AADGreeksEngine::bsm_european (autodiff::var 反向模式)
//   - 测量总耗时、每期权耗时、吞吐量 (options/s)
//   - 同时测量 Analytic Greeks 作为对照 (展示 GreeksFactory Auto 策略的正确性:
//     vanilla 应使用 Analytic 而非 AAD)
//
// 验收门槛: AAD < 100ms (Phase 3 审计 E1), Analytic < 10ms (参考)
//
// 输出格式:
//   [BENCH] E1 Greeks Batch (N=10000)
//   [BENCH] AAD:   total=XXX ms, per_option=XX μs, throughput=XXX opt/s
//   [BENCH] Analytic: total=XXX ms, per_option=XX μs, throughput=XXX opt/s
//   [PASS/FAIL] AAD threshold 100ms: ...

#include "cpphub/risk/greeks/aad_greeks.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using namespace cpphub::v1;

struct OptionSpec {
    Real S, K, T, r, q, sigma;
    bool is_call;
};

int main() {
    constexpr Size N = 10000;
    constexpr Real THRESHOLD_AAD_MS = 100.0;
    constexpr Real THRESHOLD_ANALYTIC_MS = 10.0;

    // 生成确定性随机期权组合 (seed=42 保证复现)
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<Real> dist_S(80.0, 120.0);
    std::uniform_real_distribution<Real> dist_K(70.0, 130.0);
    std::uniform_real_distribution<Real> dist_T(0.25, 2.0);
    std::uniform_real_distribution<Real> dist_r(0.01, 0.06);
    std::uniform_real_distribution<Real> dist_q(0.0, 0.03);
    std::uniform_real_distribution<Real> dist_sigma(0.1, 0.4);
    std::bernoulli_distribution dist_call(0.5);

    std::vector<OptionSpec> options;
    options.reserve(N);
    for (Size i = 0; i < N; ++i) {
        options.push_back({
            dist_S(rng), dist_K(rng), dist_T(rng),
            dist_r(rng), dist_q(rng), dist_sigma(rng),
            dist_call(rng)
        });
    }

    // 预分配结果数组 (避免在循环内分配影响测量)
    std::vector<AADGreeks> aad_results(N);
    std::vector<AnalyticGreeks> analytic_results(N);

    // ===== AAD 基准 =====
    auto t0 = std::chrono::steady_clock::now();
    for (Size i = 0; i < N; ++i) {
        const auto& o = options[i];
        aad_results[i] = AADGreeksEngine::bsm_european(
            o.S, o.K, o.T, o.r, o.q, o.sigma, o.is_call);
    }
    auto t1 = std::chrono::steady_clock::now();
    double aad_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double aad_per_us = aad_ms * 1000.0 / static_cast<double>(N);
    double aad_throughput = static_cast<double>(N) / (aad_ms / 1000.0);

    // ===== Analytic 基准 (对照) =====
    auto t2 = std::chrono::steady_clock::now();
    for (Size i = 0; i < N; ++i) {
        const auto& o = options[i];
        analytic_results[i] = AnalyticGreeksEngine::bsm_european(
            o.S, o.K, o.T, o.r, o.q, o.sigma, o.is_call);
    }
    auto t3 = std::chrono::steady_clock::now();
    double analytic_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();
    double analytic_per_us = analytic_ms * 1000.0 / static_cast<double>(N);
    double analytic_throughput = static_cast<double>(N) / (analytic_ms / 1000.0);

    // ===== 正确性交叉验证 (AAD vs Analytic, 抽样 100 个) =====
    Real max_rel_err = 0.0;
    Size sample = std::min<Size>(100, N);
    for (Size i = 0; i < sample; ++i) {
        Real p_aad = aad_results[i].price;
        Real p_ana = analytic_results[i].price;
        Real denom = std::max(std::abs(p_ana), Real(1e-10));
        Real rel = std::abs(p_aad - p_ana) / denom;
        max_rel_err = std::max(max_rel_err, rel);
    }

    // ===== 输出 =====
    std::printf("[BENCH] E1 Greeks Batch (N=%zu)\n", static_cast<size_t>(N));
    std::printf("[BENCH] AAD:      total=%7.2f ms, per_option=%7.2f us, throughput=%10.0f opt/s\n",
                aad_ms, aad_per_us, aad_throughput);
    std::printf("[BENCH] Analytic: total=%7.2f ms, per_option=%7.2f us, throughput=%10.0f opt/s\n",
                analytic_ms, analytic_per_us, analytic_throughput);
    std::printf("[BENCH] AAD vs Analytic max_rel_err (price, %zu samples): %.3e\n",
                static_cast<size_t>(sample), max_rel_err);

    bool aad_pass = aad_ms < THRESHOLD_AAD_MS;
    bool analytic_pass = analytic_ms < THRESHOLD_ANALYTIC_MS;
    std::printf("[PASS] AAD threshold %.0fms:      %s (%.2f ms)\n",
                THRESHOLD_AAD_MS, aad_pass ? "PASS" : "FAIL", aad_ms);
    std::printf("[PASS] Analytic threshold %.0fms:  %s (%.2f ms)\n",
                THRESHOLD_ANALYTIC_MS, analytic_pass ? "PASS" : "FAIL", analytic_ms);

    // 正确性门禁: AAD vs Analytic 相对误差应 < 1e-10 (两者数学等价)
    bool correctness_pass = max_rel_err < 1e-10;
    std::printf("[PASS] AAD vs Analytic correctness (1e-10): %s (%.3e)\n",
                correctness_pass ? "PASS" : "FAIL", max_rel_err);

    // 返回 0 仅当性能与正确性均通过
    int exit_code = 0;
    if (!aad_pass) {
        std::printf("[WARN] AAD 未达门槛 — 注意 GreeksFactory.Auto 对 vanilla 已用 Analytic,AAD 仅用于复杂模型\n");
        exit_code = 1;
    }
    if (!correctness_pass) exit_code = 2;
    return exit_code;
}
