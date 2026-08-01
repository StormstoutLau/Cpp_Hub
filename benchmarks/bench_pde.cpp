// Phase 3 验收项: PDE 数值求解性能基准
// 补齐 PHASE3_SPEC §6 验收项: "PDE 400x2000 < 100ms"
//
// 基准设计:
//   - B1: 1D Crank-Nicolson 欧式期权, 400 空间 x 2000 时间步 (SPEC 验收项)
//   - B2: 1D Rannacher smoothing 数字期权, 400x2000 (非光滑 payoff)
//   - B3: 1D American PSOR, 400x2000 (美式期权)
//   - B4: 2D ADI Heston CraigSneyd, 100x80x200 (默认配置)
//   - B5: 2D ADI Heston CraigSneyd, 200x150x400 (精细配置)
//
// 验收门槛:
//   B1 < 100ms (PHASE3_SPEC §6 硬性要求)
//   B2-B5: 仅记录性能, 无硬性门槛 (供后续优化参考)
//
// 输出格式:
//   [BENCH] B1 1D CN European 400x2000: price=X.XXXX time=XX ms [PASS/FAIL <100ms]
//   ...

#include "cpphub/pricing/pde/pde_engine.hpp"
#include "cpphub/pricing/pde/pde_engine_2d.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include <chrono>
#include <cstdio>
#include <vector>

using namespace cpphub::v1;

namespace {

struct BenchResult {
    double ms;
    Real price;
    bool pass;
};

template <typename Fn>
BenchResult run_bench(Fn fn, int warmup = 1, int repeats = 3) {
    // warmup
    for (int i = 0; i < warmup; ++i) (void)fn();

    double best_ms = 1e18;
    Real last_price = 0.0;
    for (int i = 0; i < repeats; ++i) {
        auto t0 = std::chrono::high_resolution_clock::now();
        last_price = fn();
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best_ms) best_ms = ms;
    }
    return {best_ms, last_price, false};
}

}  // namespace

int main() {
    std::printf("[BENCH] PDE Performance Benchmark (Cpp_Hub v1.4)\n");
    std::printf("[BENCH] ============================================\n");

    int total_pass = 0, total_tests = 0;

    // B1: 1D Crank-Nicolson European Call, 400x2000 (PHASE3_SPEC 验收项)
    {
        PDEEngineConfig cfg;
        cfg.n_spatial = 400;
        cfg.n_time = 2000;
        cfg.scheme = FDMSchemeType::CrankNicolson;
        PDEEngine engine(cfg);

        CallPayOff payoff(100.0);
        Real S0 = 100, K = 100, T = 1, rate = 0.05, q = 0.02, sigma = 0.2;
        // 计算 BSM 解析解作为精度参考
        Real d1 = (std::log(S0/K) + (rate-q+0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
        Real d2 = d1 - sigma*std::sqrt(T);
        Real analytic = S0*std::exp(-q*T)*normal_cdf(d1) - K*std::exp(-rate*T)*normal_cdf(d2);

        auto fn = [&]() {
            return engine.price_european(payoff, S0, K, T, rate, q, sigma);
        };
        auto res = run_bench(fn);
        bool pass = res.ms < 100.0;
        std::printf("[BENCH] B1 1D CN European 400x2000: price=%.4f (BSM=%.4f) "
                    "time=%.2f ms %s\n",
                    res.price, analytic, res.ms,
                    pass ? "[PASS <100ms]" : "[FAIL >100ms]");
        ++total_tests;
        if (pass) ++total_pass;
    }

    // B2: 1D Rannacher smoothing Digital Call, 400x2000 (非光滑 payoff)
    {
        PDEEngineConfig cfg;
        cfg.n_spatial = 400;
        cfg.n_time = 2000;
        cfg.scheme = FDMSchemeType::RannacherSmoothing;
        PDEEngine engine(cfg);

        DigitalCallPayOff payoff(100.0, 1.0);
        Real S0 = 100, K = 100, T = 1, rate = 0.05, q = 0.02, sigma = 0.2;

        auto fn = [&]() {
            return engine.price_european(payoff, S0, K, T, rate, q, sigma);
        };
        auto res = run_bench(fn);
        std::printf("[BENCH] B2 1D Rannacher Digital 400x2000: price=%.4f "
                    "time=%.2f ms\n",
                    res.price, res.ms);
    }

    // B3: 1D American Put PSOR, 400x2000 (美式期权)
    {
        PDEEngineConfig cfg;
        cfg.n_spatial = 400;
        cfg.n_time = 2000;
        cfg.scheme = FDMSchemeType::CrankNicolson;
        PDEEngine engine(cfg);

        PutPayOff payoff(100.0);
        Real S0 = 100, K = 100, T = 1, rate = 0.05, q = 0.0, sigma = 0.2;

        auto fn = [&]() {
            return engine.price_american(payoff, S0, K, T, rate, q, sigma);
        };
        auto res = run_bench(fn);
        std::printf("[BENCH] B3 1D American PSOR 400x2000: price=%.4f "
                    "time=%.2f ms\n",
                    res.price, res.ms);
    }

    // B4: 2D ADI Heston CraigSneyd, 100x80x200 (默认配置)
    {
        PDEEngine2DConfig cfg;
        cfg.n_x = 100;
        cfg.n_v = 80;
        cfg.n_time = 200;
        cfg.scheme = ADISchemeType::CraigSneyd;
        PDEEngine2D engine(cfg);

        HestonPDEParams hp{1.5, 0.04, 0.3, -0.5, 0.04, 0.0, 1.0, 100.0, 100.0, 0.04};

        auto fn = [&]() { return engine.price(hp); };
        auto res = run_bench(fn);
        std::printf("[BENCH] B4 2D ADI CS 100x80x200 Heston: price=%.4f "
                    "time=%.2f ms\n",
                    res.price, res.ms);
    }

    // B5: 2D ADI Heston CraigSneyd, 200x150x400 (精细配置)
    {
        PDEEngine2DConfig cfg;
        cfg.n_x = 200;
        cfg.n_v = 150;
        cfg.n_time = 400;
        cfg.scheme = ADISchemeType::CraigSneyd;
        PDEEngine2D engine(cfg);

        HestonPDEParams hp{1.5, 0.04, 0.3, -0.5, 0.04, 0.0, 1.0, 100.0, 100.0, 0.04};

        auto fn = [&]() { return engine.price(hp); };
        auto res = run_bench(fn);
        std::printf("[BENCH] B5 2D ADI CS 200x150x400 Heston: price=%.4f "
                    "time=%.2f ms\n",
                    res.price, res.ms);
    }

    // B6: 2D ADI Heston HundsdorferVerwer, 100x80x200 (HV vs CS 对照)
    {
        PDEEngine2DConfig cfg;
        cfg.n_x = 100;
        cfg.n_v = 80;
        cfg.n_time = 200;
        cfg.scheme = ADISchemeType::HundsdorferVerwer;
        PDEEngine2D engine(cfg);

        HestonPDEParams hp{1.5, 0.04, 0.3, -0.5, 0.04, 0.0, 1.0, 100.0, 100.0, 0.04};

        auto fn = [&]() { return engine.price(hp); };
        auto res = run_bench(fn);
        std::printf("[BENCH] B6 2D ADI HV 100x80x200 Heston: price=%.4f "
                    "time=%.2f ms\n",
                    res.price, res.ms);
    }

    std::printf("[BENCH] ============================================\n");
    std::printf("[BENCH] Summary: %d/%d acceptance tests passed\n",
                total_pass, total_tests);
    return (total_pass == total_tests) ? 0 : 1;
}
