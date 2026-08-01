// Phase 4 LITE - M2: Cpp_Hub Python 绑定 (nanobind)
//
// 覆盖核心模块:
//   - core: normal_pdf/cdf/inv_cdf
//   - bsm: BSM 欧式期权定价
//   - heston: Heston 特征函数
//   - greeks: Analytic + AAD Greeks
//   - var: Historical/Parametric/MC VaR
//   - cos: COS 方法定价引擎 (GBM/Heston/Bates/VG)
//   - bates: Bates 模型 (CF + 过程模拟)
//   - vg: Variance Gamma 模型 (CF + 过程模拟)
//   - cev: CEV 模型 (解析定价 + 过程模拟)
//   - sabr: SABR Hagan 隐含波动率
//   - calibration: Heston/SABR 标定器 + BSM 隐含波动率反推
//
// 构建: pip install . (需 nanobind)
// 测试: pytest python/tests/

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/pricing/analytic/bates_cf.hpp"
#include "cpphub/pricing/analytic/vg_analytic.hpp"
#include "cpphub/pricing/analytic/cev_analytic.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/risk/greeks/aad_greeks.hpp"
#include "cpphub/risk/var/historical_var.hpp"
#include "cpphub/risk/var/parametric_var.hpp"
#include "cpphub/models/diffusion/process.hpp"
#include "cpphub/models/diffusion/bates.hpp"
#include "cpphub/models/diffusion/variance_gamma.hpp"
#include "cpphub/models/diffusion/cev.hpp"

// v1.3 扩展: PayOff / PDE / Tree / MC 引擎 + Levy CF 扩展
#include "cpphub/instruments/payoff/payoff.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
#include "cpphub/pricing/pde/pde_engine.hpp"
#include "cpphub/pricing/tree/binomial.hpp"
#include "cpphub/pricing/tree/trinomial.hpp"
#include "cpphub/pricing/monte_carlo/mc_engine.hpp"
#include "cpphub/pricing/monte_carlo/path_generator.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/pricing/monte_carlo/path_dependent_payoffs.hpp"
#include <nanobind/stl/tuple.h>

// v1.4 扩展: 模型过程 (Heston/HestonQE/RoughHeston/RoughBergomi/SABR)
//        + Risk 模块 (ES/MCVaR/Backtesting) + GreeksFactory
#include "cpphub/models/diffusion/heston.hpp"
#include "cpphub/models/diffusion/heston_qe.hpp"
#include "cpphub/models/diffusion/rough_heston.hpp"
#include "cpphub/models/diffusion/rough_bergomi.hpp"
#include "cpphub/models/diffusion/sabr.hpp"
#include "cpphub/risk/var/expected_shortfall.hpp"
#include "cpphub/risk/var/mc_var.hpp"
#include "cpphub/risk/var/backtesting.hpp"
#include "cpphub/risk/greeks/greeks_factory.hpp"
#include <nanobind/stl/optional.h>

// ---------------------------------------------------------------------------
// 重要: calibrator.hpp 在 cpphub::v1 中重新定义了 HestonParams, SABRParams,
// BatesParams, 与 heston.hpp / sabr_hagan.hpp / bates.hpp 中的同名 struct 冲突
// (字段数量和顺序不同)。解决方案: 通过宏重命名 calibrator.hpp 内部的三个 struct,
// 避免冲突。Python 侧仍然以原始名称暴露。
// ---------------------------------------------------------------------------
#define HestonParams CalibHestonParams
#define SABRParams CalibSABRParams
#define BatesParams CalibBatesParams
#include "cpphub/calibration/calibrator.hpp"
#undef HestonParams
#undef SABRParams
#undef BatesParams

namespace nb = nanobind;
using namespace cpphub::v1;

NB_MODULE(_core, m) {
    m.doc() = "Cpp_Hub: Quantitative finance library (BSM/Heston/MC/Greeks/VaR/COS/Bates/VG/CEV/SABR/Calibration)";

    // =========================================================================
    // 1. Core 数学函数
    // =========================================================================
    m.def("normal_pdf", &normal_pdf, nb::arg("x"),
          "Standard normal PDF: φ(x) = (1/√(2π)) exp(-x²/2)");
    m.def("normal_cdf", &normal_cdf, nb::arg("x"),
          "Standard normal CDF: Φ(x) = 0.5 * erfc(-x/√2)");
    m.def("inv_normal_cdf",
          [](double p) {
              if (p <= 0.0 || p >= 1.0)
                  throw std::invalid_argument("p must be in (0, 1)");
              return inv_normal_cdf(p);
          },
          nb::arg("p"),
          "Inverse standard normal CDF: Φ⁻¹(p)");

    // =========================================================================
    // 2. BSM 欧式期权定价
    // =========================================================================
    m.def("bsm_price",
          [](double S, double K, double T, double r, double q,
             double sigma, bool is_call) {
              auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call);
              return g.price;
          },
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("is_call") = true,
          "Black-Scholes-Merton European option price");

    m.def("bsm_delta",
          [](double S, double K, double T, double r, double q,
             double sigma, bool is_call) {
              return AnalyticGreeksEngine::bsm_delta(S, K, T, r, q, sigma, is_call);
          },
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("is_call") = true,
          "BSM Delta: Call=N(d1), Put=N(d1)-1");

    m.def("bsm_gamma",
          &AnalyticGreeksEngine::bsm_gamma,
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          "BSM Gamma (Call = Put)");

    m.def("bsm_vega",
          &AnalyticGreeksEngine::bsm_vega,
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          "BSM Vega (Call = Put)");

    m.def("bsm_theta",
          [](double S, double K, double T, double r, double q,
             double sigma, bool is_call) {
              return AnalyticGreeksEngine::bsm_theta(S, K, T, r, q, sigma, is_call);
          },
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("is_call") = true,
          "BSM Theta");

    m.def("bsm_rho",
          [](double S, double K, double T, double r, double q,
             double sigma, bool is_call) {
              return AnalyticGreeksEngine::bsm_rho(S, K, T, r, q, sigma, is_call);
          },
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("is_call") = true,
          "BSM Rho");

    // 完整 Greeks (返回 dict)
    m.def("bsm_greeks",
          [](double S, double K, double T, double r, double q,
             double sigma, bool is_call) {
              auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call);
              nb::dict d;
              d["price"] = g.price;
              d["delta"] = g.delta;
              d["gamma"] = g.gamma;
              d["vega"] = g.vega;
              d["theta"] = g.theta;
              d["rho"] = g.rho;
              d["vanna"] = g.vanna;
              d["vomma"] = g.vomma;
              return d;
          },
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("is_call") = true,
          "BSM full Greeks (returns dict)");

    // =========================================================================
    // 3. AAD Greeks (使用 autodiff 库)
    // =========================================================================
    m.def("aad_greeks_bsm",
          [](double S, double K, double T, double r, double q,
             double sigma, bool is_call) {
              auto g = AADGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call);
              nb::dict d;
              d["price"] = g.price;
              d["delta"] = g.delta;
              d["vega"] = g.vega;
              d["rho"] = g.rho;
              d["theta"] = g.theta;
              d["gamma"] = g.gamma;
              return d;
          },
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("is_call") = true,
          "AAD Greeks for BSM (reverse-mode autodiff)");

    // =========================================================================
    // 4. Heston 特征函数
    // =========================================================================
    m.def("heston_cf",
          [](std::complex<double> u, double tau, double S0,
             double v0, double kappa, double theta, double sigma,
             double rho, double r, double q) {
              HestonCFParams p{v0, kappa, theta, sigma, rho, r, q};
              return heston_characteristic_function(u, tau, S0, p);
          },
          nb::arg("u"), nb::arg("tau"), nb::arg("S0"),
          nb::arg("v0"), nb::arg("kappa"), nb::arg("theta"),
          nb::arg("sigma"), nb::arg("rho"),
          nb::arg("r") = 0.0, nb::arg("q") = 0.0,
          "Heston characteristic function φ(u; T, S)");

    // =========================================================================
    // 5. VaR 模块
    // =========================================================================

    // 历史模拟 VaR
    m.def("historical_var",
          [](const std::vector<double>& pnl, double confidence,
             const std::string& interpolation) {
              HistoricalVaR hv(pnl, confidence, 1);
              QuantileInterpolation interp;
              if (interpolation == "linear") interp = QuantileInterpolation::Linear;
              else if (interpolation == "conservative") interp = QuantileInterpolation::Conservative;
              else interp = QuantileInterpolation::Empirical;
              return hv.var(interp);
          },
          nb::arg("pnl"), nb::arg("confidence") = 0.99,
          nb::arg("interpolation") = "linear",
          "Historical VaR (linear/conservative/empirical)");

    // 参数化 VaR (Normal)
    m.def("parametric_var_normal",
          [](double mean, double std, double confidence) {
              PortfolioStats stats;
              stats.mean = mean;
              stats.variance = std * std;
              ParametricVaR pv(stats, confidence, 1);
              return pv.var(ParametricMethod::Normal);
          },
          nb::arg("mean"), nb::arg("std"),
          nb::arg("confidence") = 0.99,
          "Parametric VaR (Normal distribution assumption)");

    // 加权历史 VaR (BRM decay)
    m.def("weighted_var",
          [](const std::vector<double>& pnl, double decay,
             double confidence) {
              HistoricalVaR hv(pnl, confidence, 1);
              return hv.weighted_var(decay, QuantileInterpolation::Linear);
          },
          nb::arg("pnl"), nb::arg("decay") = 0.99,
          nb::arg("confidence") = 0.99,
          "Weighted historical VaR (BRM decay)");

    // =========================================================================
    // 6. Philox4x64 随机数生成器 (供过程模拟使用)
    // =========================================================================
    nb::class_<Philox4x64>(m, "Philox4x64",
        "Philox4x64 counter-based RNG (10 rounds, 64-bit output)")
        .def(nb::init<uint64_t, uint64_t>(),
             nb::arg("seed") = 0, nb::arg("stream") = 0,
             "Construct Philox4x64 with seed and stream id")
        .def("__call__", &Philox4x64::operator(),
             "Draw next 64-bit unsigned integer")
        .def("discard", &Philox4x64::discard, nb::arg("n"),
             "Advance the state by n steps")
        .def("next4", &Philox4x64::next4,
             "Draw four 64-bit unsigned integers (returns list)");

    // =========================================================================
    // 7. COS 定价引擎 (Fang-Oosterlee 2008)
    // =========================================================================

    // COSEngine::Config 嵌套结构
    nb::class_<COSEngine::Config>(m, "COSConfig",
        "COS method configuration: n_terms, truncation L, optional [a,b]")
        .def(nb::init<>())
        .def_rw("n_terms", &COSEngine::Config::n_terms,
                "Number of cosine terms (typical 128-512)")
        .def_rw("L", &COSEngine::Config::L,
                "Truncation multiple (L std devs, typical 10-15)")
        .def_rw("a", &COSEngine::Config::a,
                "Optional lower bound (NaN = auto-estimate)")
        .def_rw("b", &COSEngine::Config::b,
                "Optional upper bound (NaN = auto-estimate)");

    // COSEngine 类: 通过 Python callable 构造 (CharFn = std::function<Complex(Complex)>)
    nb::class_<COSEngine>(m, "COSEngine",
        "COS method pricing engine (Fang & Oosterlee 2008). "
        "Constructor accepts a Python callable phi(u) -> complex.")
        .def(nb::init<CharFn, double, double, double, double, COSEngine::Config>(),
             nb::arg("phi"), nb::arg("S0"), nb::arg("r"), nb::arg("q"),
             nb::arg("T"), nb::arg("cfg") = COSEngine::Config(),
             "Construct COS engine with characteristic function phi(u) -> complex")
        .def("price_call", &COSEngine::price_call, nb::arg("K"),
             "European call price at strike K")
        .def("price_put", &COSEngine::price_put, nb::arg("K"),
             "European put price at strike K")
        .def("price_calls", &COSEngine::price_calls, nb::arg("strikes"),
             "Batch call pricing (vector of strikes)")
        .def("price_puts", &COSEngine::price_puts, nb::arg("strikes"),
             "Batch put pricing (vector of strikes)")
        .def("a", &COSEngine::a, "Lower truncation bound")
        .def("b", &COSEngine::b, "Upper truncation bound")
        .def("S0", &COSEngine::S0, "Initial spot")
        .def("T", &COSEngine::T, "Time to maturity")
        .def("n_terms", &COSEngine::n_terms, "Number of cosine terms");

    // CF 工厂函数 (返回 Python callable, 可传给 COSEngine)
    m.def("make_gbm_cf",
          [](double S0, double r, double q, double sigma, double T) {
              return make_gbm_cf(S0, r, q, sigma, T);
          },
          nb::arg("S0"), nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("T"),
          "GBM (Black-Scholes) characteristic function factory");

    m.def("make_heston_cf",
          [](double S0, double r, double q,
             double v0, double kappa, double theta,
             double sigma, double rho, double T) {
              HestonCFParams p{v0, kappa, theta, sigma, rho, r, q};
              return make_heston_cf(S0, r, q, p, T);
          },
          nb::arg("S0"), nb::arg("r"), nb::arg("q"),
          nb::arg("v0"), nb::arg("kappa"), nb::arg("theta"),
          nb::arg("sigma"), nb::arg("rho"), nb::arg("T"),
          "Heston characteristic function factory");

    m.def("make_bates_cf",
          [](double S0, double r, double q,
             double v0, double kappa, double theta, double sigma, double rho,
             double lambda, double mu_J, double sigma_J, double T) {
              BatesCFParams p{v0, kappa, theta, sigma, rho,
                              lambda, mu_J, sigma_J, r, q};
              return make_bates_cf(S0, r, q, p, T);
          },
          nb::arg("S0"), nb::arg("r"), nb::arg("q"),
          nb::arg("v0"), nb::arg("kappa"), nb::arg("theta"),
          nb::arg("sigma"), nb::arg("rho"),
          nb::arg("lambda"), nb::arg("mu_J"), nb::arg("sigma_J"),
          nb::arg("T"),
          "Bates characteristic function factory");

    m.def("make_vg_cf",
          [](double S0, double r, double q,
             double sigma, double nu, double theta, double T) {
              return make_vg_cf(S0, r, q, sigma, nu, theta, T);
          },
          nb::arg("S0"), nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("nu"), nb::arg("theta"),
          nb::arg("T"),
          "Variance Gamma characteristic function factory (via make_vg_cf)");

    m.def("make_vg_cf_direct",
          [](double S0, double r, double q,
             double sigma, double nu, double theta, double T) {
              VGParams p{sigma, nu, theta};
              return make_vg_cf_direct(S0, r, q, p, T);
          },
          nb::arg("S0"), nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("nu"), nb::arg("theta"),
          nb::arg("T"),
          "Variance Gamma characteristic function factory (via make_vg_cf_direct)");

    // 便捷函数: COS 定价 (GBM)
    m.def("cos_call_gbm", &cos_call_gbm,
          nb::arg("S0"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("n_terms") = 256, nb::arg("L") = 10.0,
          "COS call price under GBM (validation against BSM)");

    m.def("cos_put_gbm", &cos_put_gbm,
          nb::arg("S0"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("n_terms") = 256, nb::arg("L") = 10.0,
          "COS put price under GBM (validation against BSM)");

    // 便捷函数: COS 定价 (Heston)
    m.def("cos_call_heston",
          [](double S0, double K, double T, double r, double q,
             double v0, double kappa, double theta, double sigma, double rho,
             cpphub::v1::Size n_terms, double L) {
              HestonCFParams hp{v0, kappa, theta, sigma, rho, r, q};
              return cos_call_heston(S0, K, T, r, q, hp, n_terms, L);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"),
          nb::arg("v0"), nb::arg("kappa"), nb::arg("theta"),
          nb::arg("sigma"), nb::arg("rho"),
          nb::arg("n_terms") = 512, nb::arg("L") = 12.0,
          "COS call price under Heston");

    m.def("cos_put_heston",
          [](double S0, double K, double T, double r, double q,
             double v0, double kappa, double theta, double sigma, double rho,
             cpphub::v1::Size n_terms, double L) {
              HestonCFParams hp{v0, kappa, theta, sigma, rho, r, q};
              return cos_put_heston(S0, K, T, r, q, hp, n_terms, L);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"),
          nb::arg("v0"), nb::arg("kappa"), nb::arg("theta"),
          nb::arg("sigma"), nb::arg("rho"),
          nb::arg("n_terms") = 512, nb::arg("L") = 12.0,
          "COS put price under Heston");

    // 便捷函数: COS 定价 (Bates) - 内部构造 CF 和 COSEngine
    m.def("cos_call_bates",
          [](double S0, double K, double T, double r, double q,
             double v0, double kappa, double theta, double sigma, double rho,
             double lambda, double mu_J, double sigma_J,
             cpphub::v1::Size n_terms, double L) {
              BatesCFParams p{v0, kappa, theta, sigma, rho,
                              lambda, mu_J, sigma_J, r, q};
              auto phi = make_bates_cf(S0, r, q, p, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_call(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"),
          nb::arg("v0"), nb::arg("kappa"), nb::arg("theta"),
          nb::arg("sigma"), nb::arg("rho"),
          nb::arg("lambda"), nb::arg("mu_J"), nb::arg("sigma_J"),
          nb::arg("n_terms") = 512, nb::arg("L") = 12.0,
          "COS call price under Bates (Heston + Merton jumps)");

    m.def("cos_put_bates",
          [](double S0, double K, double T, double r, double q,
             double v0, double kappa, double theta, double sigma, double rho,
             double lambda, double mu_J, double sigma_J,
             cpphub::v1::Size n_terms, double L) {
              BatesCFParams p{v0, kappa, theta, sigma, rho,
                              lambda, mu_J, sigma_J, r, q};
              auto phi = make_bates_cf(S0, r, q, p, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_put(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"),
          nb::arg("v0"), nb::arg("kappa"), nb::arg("theta"),
          nb::arg("sigma"), nb::arg("rho"),
          nb::arg("lambda"), nb::arg("mu_J"), nb::arg("sigma_J"),
          nb::arg("n_terms") = 512, nb::arg("L") = 12.0,
          "COS put price under Bates (Heston + Merton jumps)");

    // 便捷函数: COS 定价 (VG)
    m.def("cos_call_vg",
          [](double S0, double K, double T, double r, double q,
             double sigma, double nu, double theta,
             cpphub::v1::Size n_terms, double L) {
              VGParams p{sigma, nu, theta};
              auto phi = make_vg_cf_direct(S0, r, q, p, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_call(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("nu"), nb::arg("theta"),
          nb::arg("n_terms") = 512, nb::arg("L") = 12.0,
          "COS call price under Variance Gamma");

    m.def("cos_put_vg",
          [](double S0, double K, double T, double r, double q,
             double sigma, double nu, double theta,
             cpphub::v1::Size n_terms, double L) {
              VGParams p{sigma, nu, theta};
              auto phi = make_vg_cf_direct(S0, r, q, p, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_put(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("nu"), nb::arg("theta"),
          nb::arg("n_terms") = 512, nb::arg("L") = 12.0,
          "COS put price under Variance Gamma");

    // =========================================================================
    // 8. Bates 模型 (CF + 过程模拟)
    // =========================================================================

    nb::class_<BatesCFParams>(m, "BatesCFParams",
        "Bates model CF parameters (Heston + Merton jumps). "
        "Note: the C++ field 'lambda' is exposed as 'lambda_' in Python "
        "(lambda is a Python keyword).")
        .def(nb::init<>())
        .def_rw("v0", &BatesCFParams::v0, "Initial variance (>0)")
        .def_rw("kappa", &BatesCFParams::kappa, "Mean reversion speed (>0)")
        .def_rw("theta", &BatesCFParams::theta, "Long-term variance (>0)")
        .def_rw("sigma", &BatesCFParams::sigma, "Vol of vol (>0)")
        .def_rw("rho", &BatesCFParams::rho, "Correlation [-1, 1]")
        .def_rw("lambda_", &BatesCFParams::lambda, "Jump intensity (>=0) [C++ field 'lambda']")
        .def_rw("mu_J", &BatesCFParams::mu_J, "Mean of log(J)")
        .def_rw("sigma_J", &BatesCFParams::sigma_J, "Std of log(J) (>0)")
        .def_rw("r", &BatesCFParams::r, "Risk-free rate")
        .def_rw("q", &BatesCFParams::q, "Dividend yield");

    m.def("bates_characteristic_function",
          [](std::complex<double> u, double tau, double S0,
             double v0, double kappa, double theta, double sigma, double rho,
             double lambda, double mu_J, double sigma_J,
             double r, double q) {
              BatesCFParams p{v0, kappa, theta, sigma, rho,
                              lambda, mu_J, sigma_J, r, q};
              return bates_characteristic_function(u, tau, S0, p);
          },
          nb::arg("u"), nb::arg("tau"), nb::arg("S0"),
          nb::arg("v0"), nb::arg("kappa"), nb::arg("theta"),
          nb::arg("sigma"), nb::arg("rho"),
          nb::arg("lambda"), nb::arg("mu_J"), nb::arg("sigma_J"),
          nb::arg("r") = 0.0, nb::arg("q") = 0.0,
          "Bates characteristic function φ(u; τ, S)");

    m.def("bates_jump_compensation", &bates_jump_compensation,
          nb::arg("mu_J"), nb::arg("sigma_J"),
          "Jump compensation m = E[J-1] = exp(μ_J + σ_J²/2) - 1");

    m.def("merton_jump_cf",
          [](std::complex<double> u, double tau,
             double lambda, double mu_J, double sigma_J) {
              return merton_jump_cf(u, tau, lambda, mu_J, sigma_J);
          },
          nb::arg("u"), nb::arg("tau"),
          nb::arg("lambda"), nb::arg("mu_J"), nb::arg("sigma_J"),
          "Merton jump characteristic function φ_J(u, τ)");

    // HestonScheme 枚举 (供 BatesProcess 构造使用)
    nb::enum_<HestonScheme>(m, "HestonScheme",
        "Heston/Bates discretization scheme")
        .value("Euler", HestonScheme::Euler)
        .value("FullTruncation", HestonScheme::FullTruncation)
        .value("QE_M", HestonScheme::QE_M)
        .value("Exact", HestonScheme::Exact);

    // BatesParams 结构 (用于过程模拟, 与 BatesCFParams 字段不同: 含 S0)
    nb::class_<BatesParams>(m, "BatesParams",
        "Bates process parameters (includes S0, for path simulation). "
        "Note: the C++ field 'lambda' is exposed as 'lambda_' in Python.")
        .def(nb::init<>())
        .def_rw("S0", &BatesParams::S0, "Initial spot (>0)")
        .def_rw("v0", &BatesParams::v0, "Initial variance (>0)")
        .def_rw("kappa", &BatesParams::kappa, "Mean reversion speed (>0)")
        .def_rw("theta", &BatesParams::theta, "Long-term variance (>0)")
        .def_rw("sigma", &BatesParams::sigma, "Vol of vol (>0)")
        .def_rw("rho", &BatesParams::rho, "Correlation [-1, 1]")
        .def_rw("lambda_", &BatesParams::lambda, "Jump intensity (>=0) [C++ field 'lambda']")
        .def_rw("mu_J", &BatesParams::mu_J, "Mean of log(J)")
        .def_rw("sigma_J", &BatesParams::sigma_J, "Std of log(J) (>0)")
        .def_rw("r", &BatesParams::r, "Risk-free rate")
        .def_rw("q", &BatesParams::q, "Dividend yield");

    // BatesProcess 类
    nb::class_<BatesProcess>(m, "BatesProcess",
        "Bates stochastic process (Heston SV + Merton lognormal jumps)")
        .def(nb::init<BatesParams, HestonScheme>(),
             nb::arg("p"),
             nb::arg("scheme") = HestonScheme::FullTruncation,
             "Construct BatesProcess with params and discretization scheme")
        .def("dimension", &BatesProcess::dimension,
             "State dimension (always 2: spot + variance)")
        .def("spot", &BatesProcess::spot,
             "Initial spot price")
        .def("characteristic_function",
             [](const BatesProcess& self, std::complex<double> u, double tau) {
                 return self.characteristic_function(u, tau);
             },
             nb::arg("u"), nb::arg("tau"),
             "Characteristic function φ(u; τ) at the process's S0/r/q")
        .def("generate_path",
             [](const BatesProcess& self, double T, cpphub::v1::Size n_steps,
                std::vector<double>& path, Philox4x64& rng) {
                 if (path.size() < n_steps + 1) path.resize(n_steps + 1);
                 std::span<double> sp(path.data(), path.size());
                 self.generate_path(T, n_steps, sp, rng);
                 return path;
             },
             nb::arg("T"), nb::arg("n_steps"),
             nb::arg("path").none(false), nb::arg("rng"),
             "Generate price path; path buffer resized to n_steps+1, returns path")
        .def("evolve",
             [](const BatesProcess& self, double S, double v, double dt,
                double z1, double z2, double jump_log_sum) {
                 double v_mut = v;
                 double S_new = self.evolve(S, v_mut, dt, z1, z2, jump_log_sum);
                 return nb::make_tuple(S_new, v_mut);
             },
             nb::arg("S"), nb::arg("v"), nb::arg("dt"),
             nb::arg("z1"), nb::arg("z2"), nb::arg("jump_log_sum"),
             "Single-step evolve; returns (S_new, v_new)")
        .def("sample_lognormal_jump", &BatesProcess::sample_lognormal_jump,
             nb::arg("rng"),
             "Sample lognormal jump J = exp(μ_J + σ_J * Z)")
        .def_static("sample_poisson",
                    [](double lambda_dt, Philox4x64& rng) {
                        return BatesProcess::sample_poisson(lambda_dt, rng);
                    },
                    nb::arg("lambda_dt"), nb::arg("rng"),
                    "Sample Poisson(lambda*dt) via Knuth algorithm")
        .def("params",
             [](const BatesProcess& self) { return self.params(); },
             "Return a copy of the BatesParams")
        .def("scheme", &BatesProcess::scheme,
             "Return the discretization scheme");

    // =========================================================================
    // 9. Variance Gamma 模型
    // =========================================================================

    nb::class_<VGParams>(m, "VGParams",
        "Variance Gamma parameters: sigma, nu, theta")
        .def(nb::init<>())
        .def_rw("sigma", &VGParams::sigma, "Brownian volatility (>0)")
        .def_rw("nu", &VGParams::nu, "Gamma process variance rate (>0)")
        .def_rw("theta", &VGParams::theta, "Drift (skewness source)");

    m.def("vg_omega", &vg_omega,
          nb::arg("sigma"), nb::arg("nu"), nb::arg("theta"),
          "Martingale correction ω = (1/ν) ln(1 - θν - σ²ν/2)");

    m.def("vg_characteristic_function",
          [](std::complex<double> u, double tau, double S0,
             double r, double q,
             double sigma, double nu, double theta) {
              VGParams p{sigma, nu, theta};
              return vg_characteristic_function(u, tau, S0, r, q, p);
          },
          nb::arg("u"), nb::arg("tau"), nb::arg("S0"),
          nb::arg("r") = 0.0, nb::arg("q") = 0.0,
          nb::arg("sigma"), nb::arg("nu"), nb::arg("theta"),
          "Variance Gamma characteristic function φ(u; τ, S)");

    m.def("vg_cumulant_mean", &vg_cumulant_mean,
          nb::arg("tau"), nb::arg("params"),
          "VG cumulant mean E[X_T] = θT");

    m.def("vg_cumulant_variance", &vg_cumulant_variance,
          nb::arg("tau"), nb::arg("params"),
          "VG cumulant variance Var[X_T] = (σ² + θ²ν)T");

    m.def("vg_cumulant_skewness", &vg_cumulant_skewness,
          nb::arg("tau"), nb::arg("params"),
          "VG standardized skewness");

    m.def("vg_cumulant_kurtosis_excess", &vg_cumulant_kurtosis_excess,
          nb::arg("tau"), nb::arg("params"),
          "VG excess kurtosis");

    // VarianceGammaProcess 类
    nb::class_<VarianceGammaProcess>(m, "VarianceGammaProcess",
        "Variance Gamma pure-jump Lévy process (Madan-Carr-Chang 1998)")
        .def(nb::init<VGParams, double, double, double>(),
             nb::arg("p"), nb::arg("S0"),
             nb::arg("r") = 0.0, nb::arg("q") = 0.0,
             "Construct VG process with params, spot, rate, dividend")
        .def("dimension", &VarianceGammaProcess::dimension,
             "State dimension (always 1)")
        .def("spot", &VarianceGammaProcess::spot, "Initial spot")
        .def("characteristic_function",
             [](const VarianceGammaProcess& self, std::complex<double> u, double tau) {
                 return self.characteristic_function(u, tau);
             },
             nb::arg("u"), nb::arg("tau"),
             "Characteristic function φ(u; τ)")
        .def("generate_path",
             [](const VarianceGammaProcess& self, double T,
                cpphub::v1::Size n_steps,
                std::vector<double>& path, Philox4x64& rng) {
                 if (path.size() < n_steps + 1) path.resize(n_steps + 1);
                 std::span<double> sp(path.data(), path.size());
                 self.generate_path(T, n_steps, sp, rng);
                 return path;
             },
             nb::arg("T"), nb::arg("n_steps"),
             nb::arg("path").none(false), nb::arg("rng"),
             "Generate price path; returns path buffer (resized to n_steps+1)")
        .def("vg_increment",
             [](const VarianceGammaProcess& self, double dt, double dG, double Z) {
                 return self.vg_increment(dt, dG, Z);
             },
             nb::arg("dt"), nb::arg("dG"), nb::arg("Z"),
             "VG increment ΔX = θ·ΔG + σ·sqrt(ΔG)·Z")
        .def("sample_gamma_increment",
             [](const VarianceGammaProcess& self, double dt, Philox4x64& rng) {
                 return self.sample_gamma_increment(dt, rng);
             },
             nb::arg("dt"), nb::arg("rng"),
             "Sample Gamma(dt/ν, ν) increment")
        .def("params",
             [](const VarianceGammaProcess& self) { return self.params(); },
             "Return a copy of VGParams")
        .def("S0", &VarianceGammaProcess::S0, "Initial spot")
        .def("r", &VarianceGammaProcess::r, "Risk-free rate")
        .def("q", &VarianceGammaProcess::q, "Dividend yield")
        .def("omega", &VarianceGammaProcess::omega,
             "Martingale correction ω");

    // =========================================================================
    // 10. CEV 模型 (Constant Elasticity of Variance)
    // =========================================================================

    nb::class_<CEVParams>(m, "CEVParams",
        "CEV parameters: sigma (vol scale), beta (elasticity)")
        .def(nb::init<>())
        .def_rw("sigma", &CEVParams::sigma, "Volatility scale (>0)")
        .def_rw("beta", &CEVParams::beta, "Elasticity (<1 in this impl; 1=GBM)");

    m.def("cev_call_price",
          [](double S, double K, double T, double r, double q,
             double sigma, double beta) {
              CEVParams p{sigma, beta};
              return cev_call_price(S, K, T, r, q, p);
          },
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("beta"),
          "CEV European call price (Schroder 1989, noncentral chi2)");

    m.def("cev_put_price",
          [](double S, double K, double T, double r, double q,
             double sigma, double beta) {
              CEVParams p{sigma, beta};
              return cev_put_price(S, K, T, r, q, p);
          },
          nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("beta"),
          "CEV European put price (via call-put parity)");

    m.def("validate_cev_params",
          [](const CEVParams& p) { validate_cev_params(p); },
          nb::arg("params"),
          "Validate CEV params (throws on invalid)");

    // CEVScheme 枚举
    nb::enum_<CEVScheme>(m, "CEVScheme",
        "CEV discretization scheme")
        .value("EulerAbsorbing", CEVScheme::EulerAbsorbing,
               "Euler with absorbing barrier at 0 (default, general beta)")
        .value("LogEuler", CEVScheme::LogEuler,
               "log-Euler (exact for beta=1, approximate for beta<1)");

    // CEVProcess 类
    nb::class_<CEVProcess>(m, "CEVProcess",
        "CEV stochastic process: dS = (r-q)S dt + σ S^β dW")
        .def(nb::init<CEVParams, double, double, double, CEVScheme>(),
             nb::arg("p"), nb::arg("S0"),
             nb::arg("r") = 0.0, nb::arg("q") = 0.0,
             nb::arg("scheme") = CEVScheme::EulerAbsorbing,
             "Construct CEV process")
        .def("dimension", &CEVProcess::dimension, "State dimension (always 1)")
        .def("spot", &CEVProcess::spot, "Initial spot")
        .def("characteristic_function",
             [](const CEVProcess& self, std::complex<double> u, double tau) {
                 return self.characteristic_function(u, tau);
             },
             nb::arg("u"), nb::arg("tau"),
             "CEV has no closed-form CF (returns 0; use cev_call_price)")
        .def("generate_path",
             [](const CEVProcess& self, double T, cpphub::v1::Size n_steps,
                std::vector<double>& path, Philox4x64& rng) {
                 if (path.size() < n_steps + 1) path.resize(n_steps + 1);
                 std::span<double> sp(path.data(), path.size());
                 self.generate_path(T, n_steps, sp, rng);
                 return path;
             },
             nb::arg("T"), nb::arg("n_steps"),
             nb::arg("path").none(false), nb::arg("rng"),
             "Generate price path; returns path buffer (resized to n_steps+1)")
        .def("evolve",
             [](const CEVProcess& self, double S, double dt, double Z) {
                 return self.evolve(S, dt, Z);
             },
             nb::arg("S"), nb::arg("dt"), nb::arg("Z"),
             "Single-step evolve; returns next S")
        .def("params",
             [](const CEVProcess& self) { return self.params(); },
             "Return a copy of CEVParams")
        .def("S0", &CEVProcess::S0, "Initial spot")
        .def("r", &CEVProcess::r, "Risk-free rate")
        .def("q", &CEVProcess::q, "Dividend yield")
        .def("scheme", &CEVProcess::scheme, "Discretization scheme");

    // =========================================================================
    // 11. SABR 模型 (Hagan 2002 隐含波动率)
    //
    // 注: calibrator.hpp 在 cpphub::v1 中重新定义了 SABRParams (字段顺序:
    // alpha, beta, nu, rho), 与 sabr_hagan.hpp 不同。Python 绑定使用
    // calibrator.hpp 的版本 (此处宏替换后类型名为 CalibSABRParams)。
    // 字段通过名字访问, 顺序对 Python 用户透明。
    // =========================================================================

    nb::class_<CalibSABRParams>(m, "SABRParams",
        "SABR parameters: alpha, beta, rho, nu")
        .def(nb::init<>())
        .def_rw("alpha", &CalibSABRParams::alpha, "Initial vol σ(0) (>0)")
        .def_rw("beta", &CalibSABRParams::beta, "CEV elasticity ∈ [0, 1]")
        .def_rw("nu", &CalibSABRParams::nu, "Vol-of-vol (>0)")
        .def_rw("rho", &CalibSABRParams::rho, "Correlation ∈ (-1, 1)");

    // sabr_implied_vol_hagan: Python 接口签名 (K, F, T, params),
    // 内部调用 calibrator.hpp detail 版本 (F, K, T, params) (参数顺序相反)
    m.def("sabr_implied_vol_hagan",
          [](double K, double F, double T, const CalibSABRParams& p) {
              return cpphub::v1::detail::sabr_implied_vol_hagan(F, K, T, p);
          },
          nb::arg("K"), nb::arg("F"), nb::arg("T"), nb::arg("params"),
          "SABR Black implied vol (Hagan 2002 asymptotic formula)");

    // =========================================================================
    // 12. 标定器 (Heston / SABR) + BSM 隐含波动率
    //
    // 注: HestonParams 在 calibrator.hpp 内为 5 字段 (v0, kappa, theta,
    // sigma_v, rho), 与 heston.hpp 中的 8 字段版本不同。此处宏替换后
    // 类型名为 CalibHestonParams, Python 侧暴露为 HestonCalibParams。
    // =========================================================================

    nb::class_<CalibHestonParams>(m, "HestonCalibParams",
        "Heston calibration parameters (5-vector: v0, kappa, theta, sigma_v, rho)")
        .def(nb::init<>())
        .def_rw("v0", &CalibHestonParams::v0, "Initial variance")
        .def_rw("kappa", &CalibHestonParams::kappa, "Mean reversion speed")
        .def_rw("theta", &CalibHestonParams::theta, "Long-term variance")
        .def_rw("sigma_v", &CalibHestonParams::sigma_v, "Vol of vol")
        .def_rw("rho", &CalibHestonParams::rho, "Correlation");

    // MarketQuote 结构
    nb::class_<MarketQuote>(m, "MarketQuote",
        "Market quote: strike, maturity, price, IV, vega")
        .def(nb::init<>())
        .def_rw("strike", &MarketQuote::strike, "Strike price")
        .def_rw("maturity", &MarketQuote::maturity, "Time to maturity (years)")
        .def_rw("market_price", &MarketQuote::market_price, "Observed market price")
        .def_rw("implied_vol", &MarketQuote::implied_vol, "Market implied vol")
        .def_rw("vega", &MarketQuote::vega, "Optional vega (for vega weighting)");

    // CalibConfig 结构
    nb::class_<CalibConfig>(m, "CalibConfig",
        "Calibration configuration: DE + LM settings")
        .def(nb::init<>())
        .def_rw("de_pop_size", &CalibConfig::de_pop_size, "DE population size")
        .def_rw("de_generations", &CalibConfig::de_generations, "DE generations")
        .def_rw("lm_max_iter", &CalibConfig::lm_max_iter, "LM max iterations")
        .def_rw("ftol", &CalibConfig::ftol, "LM function tolerance")
        .def_rw("xtol", &CalibConfig::xtol, "LM parameter tolerance")
        .def_rw("seed", &CalibConfig::seed, "RNG seed (0 = use default 42)")
        .def_rw("use_de_init", &CalibConfig::use_de_init,
                "Use DE global search for initialization");

    // WeightingScheme 枚举
    nb::enum_<WeightingScheme>(m, "WeightingScheme",
        "Calibration objective weighting scheme")
        .value("PriceWeighted", WeightingScheme::PriceWeighted)
        .value("VegaWeighted", WeightingScheme::VegaWeighted)
        .value("RelativeError", WeightingScheme::RelativeError)
        .value("Mixed", WeightingScheme::Mixed);

    // CalibrationResult 结构
    nb::class_<CalibrationResult>(m, "CalibrationResult",
        "Calibration result: params, objective, residuals, diagnostics")
        .def_rw("params", &CalibrationResult::params, "Calibrated parameter vector")
        .def_rw("objective_value", &CalibrationResult::objective_value,
                "Final objective value")
        .def_rw("n_iterations", &CalibrationResult::n_iterations,
                "LM iteration count")
        .def_rw("converged", &CalibrationResult::converged,
                "Convergence flag")
        .def_rw("message", &CalibrationResult::message,
                "Status message")
        .def_rw("residuals", &CalibrationResult::residuals,
                "Per-quote residuals")
        .def("__repr__",
             [](const CalibrationResult& r) {
                 return std::string("<CalibrationResult converged=") +
                        (r.converged ? "True" : "False") +
                        " n_iter=" + std::to_string(r.n_iterations) +
                        " obj=" + std::to_string(r.objective_value) +
                        " msg='" + r.message + "'>";
             });

    // bsm_implied_vol: Newton-Raphson 反推 BSM 隐含波动率
    m.def("bsm_implied_vol",
          [](double C_market, double S, double K, double T,
             double r, double q, bool is_call,
             double tol, int max_iter) {
              return bsm_implied_vol(C_market, S, K, T, r, q, is_call, tol, max_iter);
          },
          nb::arg("C_market"), nb::arg("S"), nb::arg("K"), nb::arg("T"),
          nb::arg("r"), nb::arg("q"), nb::arg("is_call") = true,
          nb::arg("tol") = 1e-10, nb::arg("max_iter") = 50,
          "BSM implied volatility (Newton-Raphson with vega)");

    // HestonCalibrator 类
    nb::class_<HestonCalibrator>(m, "HestonCalibrator",
        "Heston model calibrator (DE global init + LM refine, IV objective)")
        .def(nb::init<>(), "Construct HestonCalibrator")
        .def("set_market", &HestonCalibrator::set_market,
             nb::arg("S"), nb::arg("r"), nb::arg("q"),
             "Set market context (spot, rate, dividend) before calibrate")
        .def("calibrate",
             [](HestonCalibrator& self,
                const std::vector<MarketQuote>& quotes,
                const CalibConfig& cfg) {
                 return self.calibrate(quotes, cfg);
             },
             nb::arg("quotes"), nb::arg("cfg") = CalibConfig{},
             "Calibrate Heston to market quotes; returns CalibrationResult")
        .def("extract_params",
             [](const HestonCalibrator& self, const std::vector<double>& x) {
                 return self.extract_params(x);
             },
             nb::arg("x"),
             "Convert flat parameter vector to HestonCalibParams")
        .def_static("default_bounds",
                    []() {
                        auto bounds = HestonCalibrator::default_bounds();
                        nb::list result;
                        for (const auto& b : bounds) {
                            result.append(nb::make_tuple(b.lower, b.upper));
                        }
                        return result;
                    },
                    "Default parameter bounds [v0, kappa, theta, sigma_v, rho] "
                    "as list of (lower, upper) tuples")
        .def_static("check_feller",
                    [](const CalibHestonParams& p) {
                        return HestonCalibrator::check_feller(p);
                    },
                    nb::arg("params"),
                    "Check Feller condition 2κθ > σ²")
        .def("name", &HestonCalibrator::name, "Calibrator name");

    // SABRCalibrator 类
    nb::class_<SABRCalibrator>(m, "SABRCalibrator",
        "SABR model calibrator (Hagan IV formula + DE + LM)")
        .def(nb::init<>(), "Construct SABRCalibrator")
        .def("set_market", &SABRCalibrator::set_market,
             nb::arg("F"), nb::arg("r"), nb::arg("q"),
             "Set market context (forward, rate, dividend) before calibrate")
        .def("calibrate",
             [](SABRCalibrator& self,
                const std::vector<MarketQuote>& quotes,
                const CalibConfig& cfg) {
                 return self.calibrate(quotes, cfg);
             },
             nb::arg("quotes"), nb::arg("cfg") = CalibConfig{},
             "Calibrate SABR to market quotes; returns CalibrationResult")
        .def("extract_params",
             [](const SABRCalibrator& self, const std::vector<double>& x) {
                 return self.extract_params(x);
             },
             nb::arg("x"),
             "Convert flat parameter vector to SABRParams")
        .def_static("default_bounds",
                    []() {
                        auto bounds = SABRCalibrator::default_bounds();
                        nb::list result;
                        for (const auto& b : bounds) {
                            result.append(nb::make_tuple(b.lower, b.upper));
                        }
                        return result;
                    },
                    "Default parameter bounds [alpha, beta, nu, rho] "
                    "as list of (lower, upper) tuples")
        .def("name", &SABRCalibrator::name, "Calibrator name");

    // =========================================================================
    // 13. v1.3 扩展: PayOff 类层次 (抽象基类 + 具体子类)
    // =========================================================================
    // PayOff 抽象基类 — 不可直接实例化, 通过子类创建
    nb::class_<PayOff>(m, "PayOff",
        "Abstract payoff base class. Use CallPayOff/PutPayOff/DigitalCallPayOff/etc.")
        .def("__call__", [](const PayOff& self, double spot) { return self(spot); },
             nb::arg("spot"), "Evaluate payoff at given spot")
        .def("name", &PayOff::name, "Payoff name string");

    nb::class_<CallPayOff, PayOff>(m, "CallPayOff",
        "Vanilla call payoff: max(S - K, 0)")
        .def(nb::init<double>(), nb::arg("strike"))
        .def("__call__", [](const CallPayOff& self, double spot) { return self(spot); },
             nb::arg("spot"))
        .def("name", &CallPayOff::name);

    nb::class_<PutPayOff, PayOff>(m, "PutPayOff",
        "Vanilla put payoff: max(K - S, 0)")
        .def(nb::init<double>(), nb::arg("strike"))
        .def("__call__", [](const PutPayOff& self, double spot) { return self(spot); },
             nb::arg("spot"))
        .def("name", &PutPayOff::name);

    nb::class_<DigitalCallPayOff, PayOff>(m, "DigitalCallPayOff",
        "Digital call: pays fixed amount if S > K")
        .def(nb::init<double, double>(), nb::arg("strike"), nb::arg("payment") = 1.0)
        .def("__call__", [](const DigitalCallPayOff& self, double spot) { return self(spot); },
             nb::arg("spot"))
        .def("name", &DigitalCallPayOff::name);

    nb::class_<DigitalPutPayOff, PayOff>(m, "DigitalPutPayOff",
        "Digital put: pays fixed amount if S < K")
        .def(nb::init<double, double>(), nb::arg("strike"), nb::arg("payment") = 1.0)
        .def("__call__", [](const DigitalPutPayOff& self, double spot) { return self(spot); },
             nb::arg("spot"))
        .def("name", &DigitalPutPayOff::name);

    nb::class_<DoubleDigitalPayOff, PayOff>(m, "DoubleDigitalPayOff",
        "Double digital: pays if lower <= S <= upper")
        .def(nb::init<double, double, double>(),
             nb::arg("lower"), nb::arg("upper"), nb::arg("payment") = 1.0)
        .def("__call__", [](const DoubleDigitalPayOff& self, double spot) { return self(spot); },
             nb::arg("spot"))
        .def("name", &DoubleDigitalPayOff::name);

    // =========================================================================
    // 14. v1.3 扩展: PDE 引擎 (有限差分法)
    // =========================================================================
    nb::enum_<FDMSchemeType>(m, "FDMSchemeType",
        "Finite difference scheme type")
        .value("ExplicitEuler", FDMSchemeType::ExplicitEuler)
        .value("ImplicitEuler", FDMSchemeType::ImplicitEuler)
        .value("CrankNicolson", FDMSchemeType::CrankNicolson);

    nb::class_<PDEEngineConfig>(m, "PDEEngineConfig",
        "Configuration for PDE finite difference engine")
        .def(nb::init<>())
        .def_rw("n_spatial", &PDEEngineConfig::n_spatial,
                "Number of spatial grid points (default 400)")
        .def_rw("n_time", &PDEEngineConfig::n_time,
                "Number of time steps (default 1000)")
        .def_rw("alpha", &PDEEngineConfig::alpha,
                "Sinh grid concentration parameter (default 0.2)")
        .def_rw("scheme", &PDEEngineConfig::scheme,
                "FDM scheme type (default CrankNicolson)")
        .def_rw("s_multiplier", &PDEEngineConfig::s_multiplier,
                "Spatial domain multiplier (default 5.0)");

    // PDE Greeks 结构体
    nb::class_<PDEEngine::Greeks>(m, "PDEGreeks",
        "PDE-computed Greeks (delta, gamma, theta)")
        .def(nb::init<>())
        .def_rw("delta", &PDEEngine::Greeks::delta)
        .def_rw("gamma", &PDEEngine::Greeks::gamma)
        .def_rw("theta", &PDEEngine::Greeks::theta)
        .def("__repr__", [](const PDEEngine::Greeks& g) {
            return "<PDEGreeks delta=" + std::to_string(g.delta) +
                   ", gamma=" + std::to_string(g.gamma) +
                   ", theta=" + std::to_string(g.theta) + ">";
        });

    nb::class_<PDEEngine>(m, "PDEEngine",
        "PDE finite difference pricing engine (1D Black-Scholes PDE)")
        .def(nb::init<PDEEngineConfig>(), nb::arg("config") = PDEEngineConfig{})
        .def("price_european",
             [](const PDEEngine& self, const PayOff& payoff,
                double S0, double K, double T, double r, double q, double sigma) {
                 return self.price_european(payoff, S0, K, T, r, q, sigma);
             },
             nb::arg("payoff"), nb::arg("S0"), nb::arg("K"), nb::arg("T"),
             nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
             "Price European option via PDE")
        .def("price_american",
             [](const PDEEngine& self, const PayOff& payoff,
                double S0, double K, double T, double r, double q, double sigma) {
                 return self.price_american(payoff, S0, K, T, r, q, sigma);
             },
             nb::arg("payoff"), nb::arg("S0"), nb::arg("K"), nb::arg("T"),
             nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
             "Price American option via PDE (PSOR)")
        .def("greeks",
             [](const PDEEngine& self, const PayOff& payoff,
                double S0, double K, double T, double r, double q, double sigma,
                bool american) {
                 return self.greeks(payoff, S0, K, T, r, q, sigma, american);
             },
             nb::arg("payoff"), nb::arg("S0"), nb::arg("K"), nb::arg("T"),
             nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
             nb::arg("american") = false,
             "Compute Greeks via numerical differentiation")
        .def("name", &PDEEngine::name)
        .def_prop_ro("config", &PDEEngine::config,
                     "Engine configuration");

    // =========================================================================
    // 15. v1.3 扩展: 二叉树/三叉树引擎
    // =========================================================================
    nb::enum_<BinomialType>(m, "BinomialType",
        "Binomial tree type")
        .value("CRR", BinomialType::CRR)
        .value("JarrowRudd", BinomialType::JarrowRudd)
        .value("Tian", BinomialType::Tian)
        .value("LeisenReimer", BinomialType::LeisenReimer);

    nb::class_<BinomialParams>(m, "BinomialParams",
        "Parameters for binomial tree engine")
        .def(nb::init<>())
        .def_rw("S0", &BinomialParams::S0)
        .def_rw("K", &BinomialParams::K)
        .def_rw("T", &BinomialParams::T)
        .def_rw("r", &BinomialParams::r)
        .def_rw("q", &BinomialParams::q)
        .def_rw("sigma", &BinomialParams::sigma)
        .def_rw("n_steps", &BinomialParams::n_steps);

    nb::class_<BinomialTreeEngine>(m, "BinomialTreeEngine",
        "Binomial tree pricing engine (CRR/JR/Tian/Leisen-Reimer)")
        .def(nb::init<BinomialParams, BinomialType>(),
             nb::arg("params"), nb::arg("type") = BinomialType::CRR)
        .def("price_european",
             [](const BinomialTreeEngine& self, const PayOff& payoff) {
                 return self.price_european(payoff);
             },
             nb::arg("payoff"), "Price European option")
        .def("price_american",
             [](const BinomialTreeEngine& self, const PayOff& payoff) {
                 return self.price_american(payoff);
             },
             nb::arg("payoff"), "Price American option (early exercise)")
        .def("price_bermudan",
             [](const BinomialTreeEngine& self, const PayOff& payoff,
                const std::vector<Size>& exercise_steps) {
                 return self.price_bermudan(payoff, exercise_steps);
             },
             nb::arg("payoff"), nb::arg("exercise_steps"),
             "Price Bermudan option with given exercise steps");

    nb::enum_<TrinomialType>(m, "TrinomialType",
        "Trinomial tree type")
        .value("Explicit", TrinomialType::Explicit)
        .value("Implicit", TrinomialType::Implicit)
        .value("Hybrid", TrinomialType::Hybrid);

    nb::class_<TrinomialParams>(m, "TrinomialParams",
        "Parameters for trinomial tree engine")
        .def(nb::init<>())
        .def_rw("S0", &TrinomialParams::S0)
        .def_rw("K", &TrinomialParams::K)
        .def_rw("T", &TrinomialParams::T)
        .def_rw("r", &TrinomialParams::r)
        .def_rw("q", &TrinomialParams::q)
        .def_rw("sigma", &TrinomialParams::sigma)
        .def_rw("n_steps", &TrinomialParams::n_steps);

    nb::class_<TrinomialTreeEngine>(m, "TrinomialTreeEngine",
        "Trinomial tree pricing engine (Explicit/Implicit/Hybrid)")
        .def(nb::init<TrinomialParams, TrinomialType>(),
             nb::arg("params"), nb::arg("type") = TrinomialType::Explicit)
        .def("price_european",
             [](const TrinomialTreeEngine& self, const PayOff& payoff) {
                 return self.price_european(payoff);
             },
             nb::arg("payoff"), "Price European option")
        .def("price_american",
             [](const TrinomialTreeEngine& self, const PayOff& payoff) {
                 return self.price_american(payoff);
             },
             nb::arg("payoff"), "Price American option");

    // =========================================================================
    // 16. v1.3 扩展: Levy 过程 CF 工厂 (CGMY / Kou / NIG)
    // =========================================================================
    // CGMY (Carr-Geman-Madan-Yor 2002) 纯跳跃 Levy 过程
    m.def("make_cgmy_cf",
          [](double S0, double r, double q,
             double C, double G, double M, double Y, double T) {
              return make_cgmy_cf(S0, r, q, C, G, M, Y, T);
          },
          nb::arg("S0"), nb::arg("r"), nb::arg("q"),
          nb::arg("C"), nb::arg("G"), nb::arg("M"), nb::arg("Y"), nb::arg("T"),
          "CGMY (Carr-Geman-Madan-Yor 2002) characteristic function factory.\n"
          "Pure-jump Levy process with parameters (C, G, M, Y), 0 < Y < 2.\n"
          "Returns a callable phi(u) -> complex.");

    // Kou (2002) 双指数跳跃扩散
    m.def("make_kou_cf",
          [](double S0, double r, double q, double sigma,
             double lambda, double p, double eta1, double eta2, double T) {
              return make_kou_cf(S0, r, q, sigma, lambda, p, eta1, eta2, T);
          },
          nb::arg("S0"), nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
          nb::arg("lambda"), nb::arg("p"), nb::arg("eta1"), nb::arg("eta2"),
          nb::arg("T"),
          "Kou (2002) double-exponential jump-diffusion CF factory.\n"
          "Parameters: sigma (diffusion), lambda (jump intensity),\n"
          "p (down-jump prob), eta1/eta2 (jump sizes). Requires eta2 < 1.\n"
          "Returns a callable phi(u) -> complex.");

    // NIG (Normal Inverse Gaussian)
    m.def("make_nig_cf",
          [](double S0, double r, double q,
             double alpha, double beta, double delta, double T) {
              return make_nig_cf(S0, r, q, alpha, beta, delta, T);
          },
          nb::arg("S0"), nb::arg("r"), nb::arg("q"),
          nb::arg("alpha"), nb::arg("beta"), nb::arg("delta"), nb::arg("T"),
          "Normal Inverse Gaussian (NIG) CF factory.\n"
          "Parameters: alpha > 0, |beta| < alpha, delta > 0.\n"
          "Returns a callable phi(u) -> complex.");

    // COS 便捷函数: CGMY / Kou / NIG
    m.def("cos_call_cgmy",
          [](double S0, double K, double T, double r, double q,
             double C, double G, double M, double Y,
             int n_terms, double L) {
              auto phi = make_cgmy_cf(S0, r, q, C, G, M, Y, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_call(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"), nb::arg("r"), nb::arg("q"),
          nb::arg("C"), nb::arg("G"), nb::arg("M"), nb::arg("Y"),
          nb::arg("n_terms") = 256, nb::arg("L") = 12.0,
          "Price European call via COS method under CGMY model");

    m.def("cos_put_cgmy",
          [](double S0, double K, double T, double r, double q,
             double C, double G, double M, double Y,
             int n_terms, double L) {
              auto phi = make_cgmy_cf(S0, r, q, C, G, M, Y, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_put(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"), nb::arg("r"), nb::arg("q"),
          nb::arg("C"), nb::arg("G"), nb::arg("M"), nb::arg("Y"),
          nb::arg("n_terms") = 256, nb::arg("L") = 12.0,
          "Price European put via COS method under CGMY model");

    m.def("cos_call_kou",
          [](double S0, double K, double T, double r, double q,
             double sigma, double lambda, double p, double eta1, double eta2,
             int n_terms, double L) {
              auto phi = make_kou_cf(S0, r, q, sigma, lambda, p, eta1, eta2, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_call(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"), nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("lambda"), nb::arg("p"),
          nb::arg("eta1"), nb::arg("eta2"),
          nb::arg("n_terms") = 256, nb::arg("L") = 12.0,
          "Price European call via COS method under Kou model");

    m.def("cos_put_kou",
          [](double S0, double K, double T, double r, double q,
             double sigma, double lambda, double p, double eta1, double eta2,
             int n_terms, double L) {
              auto phi = make_kou_cf(S0, r, q, sigma, lambda, p, eta1, eta2, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_put(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"), nb::arg("r"), nb::arg("q"),
          nb::arg("sigma"), nb::arg("lambda"), nb::arg("p"),
          nb::arg("eta1"), nb::arg("eta2"),
          nb::arg("n_terms") = 256, nb::arg("L") = 12.0,
          "Price European put via COS method under Kou model");

    m.def("cos_call_nig",
          [](double S0, double K, double T, double r, double q,
             double alpha, double beta, double delta,
             int n_terms, double L) {
              auto phi = make_nig_cf(S0, r, q, alpha, beta, delta, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_call(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"), nb::arg("r"), nb::arg("q"),
          nb::arg("alpha"), nb::arg("beta"), nb::arg("delta"),
          nb::arg("n_terms") = 256, nb::arg("L") = 12.0,
          "Price European call via COS method under NIG model");

    m.def("cos_put_nig",
          [](double S0, double K, double T, double r, double q,
             double alpha, double beta, double delta,
             int n_terms, double L) {
              auto phi = make_nig_cf(S0, r, q, alpha, beta, delta, T);
              COSEngine::Config cfg;
              cfg.n_terms = n_terms;
              cfg.L = L;
              COSEngine engine(phi, S0, r, q, T, cfg);
              return engine.price_put(K);
          },
          nb::arg("S0"), nb::arg("K"), nb::arg("T"), nb::arg("r"), nb::arg("q"),
          nb::arg("alpha"), nb::arg("beta"), nb::arg("delta"),
          nb::arg("n_terms") = 256, nb::arg("L") = 12.0,
          "Price European put via COS method under NIG model");

    // =========================================================================
    // 17. v1.3 扩展: MC 引擎 (单资产欧式 + 路径依赖)
    // =========================================================================
    nb::class_<MCConfig>(m, "MCConfig",
        "Monte Carlo engine configuration")
        .def(nb::init<>())
        .def_rw("n_paths", &MCConfig::n_paths,
                "Number of paths (antithetic: n_paths/2 pairs)")
        .def_rw("seed", &MCConfig::seed, "RNG seed")
        .def_rw("use_antithetic", &MCConfig::use_antithetic,
                "Enable antithetic variates")
        .def_rw("use_control_variate", &MCConfig::use_control_variate,
                "Enable control variates")
        .def_rw("df", &MCConfig::df,
                "Discount factor exp(-rT)");

    nb::class_<MCResult>(m, "MCResult",
        "Monte Carlo pricing result")
        .def(nb::init<>())
        .def_ro("price", &MCResult::price, "Discounted option price")
        .def_ro("std_error", &MCResult::std_error, "Standard error")
        .def_ro("ci_lower", &MCResult::ci_lower, "95% CI lower bound")
        .def_ro("ci_upper", &MCResult::ci_upper, "95% CI upper bound")
        .def_ro("n_paths", &MCResult::n_paths, "Actual path count")
        .def_ro("beta_cv", &MCResult::beta_cv,
                "Optimal control variate beta (if enabled)")
        .def_ro("variance_reduction", &MCResult::variance_reduction,
                "Variance reduction factor (vs no reduction)")
        .def("__repr__", [](const MCResult& r) {
            return "<MCResult price=" + std::to_string(r.price) +
                   ", se=" + std::to_string(r.std_error) +
                   ", n=" + std::to_string(r.n_paths) + ">";
        });

    // MC 便捷函数: 欧式期权 (单资产 GBM)
    m.def("mc_price_european",
          [](double S0, double sigma, double r, double q, double T,
             double K, bool is_call,
             uint64_t n_paths, uint64_t seed, bool use_antithetic) {
              // 构建单资产 GBM 路径生成器 (n_steps=1: 欧式只需终端)
              auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, 1);
              MultiAssetGBMPathGenerator gen(gen_cfg);

              // 构建终端 payoff
              PathPayoff payoff;
              if (is_call) {
                  payoff = [K](const std::vector<Real>& path) {
                      return std::max(path.back() - K, 0.0);
                  };
              } else {
                  payoff = [K](const std::vector<Real>& path) {
                      return std::max(K - path.back(), 0.0);
                  };
              }

              MCConfig cfg;
              cfg.n_paths = n_paths;
              cfg.seed = seed;
              cfg.use_antithetic = use_antithetic;
              cfg.df = std::exp(-r * T);

              return price_path_dependent(gen, payoff, cfg);
          },
          nb::arg("S0"), nb::arg("sigma"), nb::arg("r"), nb::arg("q"),
          nb::arg("T"), nb::arg("K"), nb::arg("is_call") = true,
          nb::arg("n_paths") = 100000, nb::arg("seed") = 42,
          nb::arg("use_antithetic") = false,
          "Price European option via Monte Carlo (single-asset GBM).\n"
          "Returns MCResult with price, std_error, CI, etc.");

    // MC 便捷函数: 算术亚式期权
    m.def("mc_price_asian_arithmetic",
          [](double S0, double sigma, double r, double q, double T,
             double K, bool is_call, Size n_steps,
             uint64_t n_paths, uint64_t seed, bool use_antithetic) {
              auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
              MultiAssetGBMPathGenerator gen(gen_cfg);

              auto payoff = make_asian_payoff(
                  K, is_call ? OptionType::Call : OptionType::Put,
                  AsianAverageType::Arithmetic);

              MCConfig cfg;
              cfg.n_paths = n_paths;
              cfg.seed = seed;
              cfg.use_antithetic = use_antithetic;
              cfg.df = std::exp(-r * T);

              return price_path_dependent(gen, payoff, cfg);
          },
          nb::arg("S0"), nb::arg("sigma"), nb::arg("r"), nb::arg("q"),
          nb::arg("T"), nb::arg("K"), nb::arg("is_call") = true,
          nb::arg("n_steps") = 252,
          nb::arg("n_paths") = 100000, nb::arg("seed") = 42,
          nb::arg("use_antithetic") = true,
          "Price arithmetic Asian option via MC (discrete monitoring).\n"
          "Returns MCResult.");

    // MC 便捷函数: 障碍期权 (Up-and-Out Call)
    m.def("mc_price_barrier_up_out_call",
          [](double S0, double sigma, double r, double q, double T,
             double K, double barrier,
             Size n_steps, uint64_t n_paths, uint64_t seed,
             bool use_antithetic) {
              auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
              MultiAssetGBMPathGenerator gen(gen_cfg);

              BarrierSpec spec;
              spec.barrier = barrier;
              spec.dir = BarrierDirection::Up;
              spec.knock = BarrierKnock::Out;
              spec.K = K;
              spec.inner_opt = OptionType::Call;
              auto payoff = make_barrier_payoff(spec);

              MCConfig cfg;
              cfg.n_paths = n_paths;
              cfg.seed = seed;
              cfg.use_antithetic = use_antithetic;
              cfg.df = std::exp(-r * T);

              return price_path_dependent(gen, payoff, cfg);
          },
          nb::arg("S0"), nb::arg("sigma"), nb::arg("r"), nb::arg("q"),
          nb::arg("T"), nb::arg("K"), nb::arg("barrier"),
          nb::arg("n_steps") = 252,
          nb::arg("n_paths") = 100000, nb::arg("seed") = 42,
          nb::arg("use_antithetic") = true,
          "Price Up-and-Out Call barrier option via MC.\n"
          "Returns MCResult.");

    // MC 便捷函数: 回望期权 (Floating Strike Call)
    m.def("mc_price_lookback_floating_call",
          [](double S0, double sigma, double r, double q, double T,
             Size n_steps, uint64_t n_paths, uint64_t seed,
             bool use_antithetic) {
              auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
              MultiAssetGBMPathGenerator gen(gen_cfg);

              // make_lookback_payoff(K, OptionType, LookbackType)
              // floating strike: K 不使用, 传 0.0
              auto payoff = make_lookback_payoff(
                  0.0, OptionType::Call, LookbackType::FloatingStrike);

              MCConfig cfg;
              cfg.n_paths = n_paths;
              cfg.seed = seed;
              cfg.use_antithetic = use_antithetic;
              cfg.df = std::exp(-r * T);

              return price_path_dependent(gen, payoff, cfg);
          },
          nb::arg("S0"), nb::arg("sigma"), nb::arg("r"), nb::arg("q"),
          nb::arg("T"), nb::arg("n_steps") = 252,
          nb::arg("n_paths") = 100000, nb::arg("seed") = 42,
          nb::arg("use_antithetic") = true,
          "Price floating-strike Lookback Call via MC.\n"
          "Returns MCResult.");

    // =========================================================================
    // 18. v1.4 扩展: 模型过程 (Heston / HestonQE / RoughHeston / RoughBergomi / SABR)
    // =========================================================================
    // Heston 过程参数 (heston.hpp 8 字段版, 与 calibrator 的 5 字段版不同)
    // 注: HestonScheme 枚举已在 §8 (Bates) 中暴露
    nb::class_<HestonParams>(m, "HestonProcessParams",
        "Heston process parameters (8 fields, includes S0/r/q for path simulation). "
        "Distinct from HestonCalibParams (5 fields, calibration-only).")
        .def(nb::init<>())
        .def_rw("S0", &HestonParams::S0, "Initial spot (>0)")
        .def_rw("v0", &HestonParams::v0, "Initial variance (>0)")
        .def_rw("kappa", &HestonParams::kappa, "Mean reversion speed (>0)")
        .def_rw("theta", &HestonParams::theta, "Long-term variance (>0)")
        .def_rw("sigma", &HestonParams::sigma, "Vol of vol (>0)")
        .def_rw("rho", &HestonParams::rho, "Correlation [-1, 1]")
        .def_rw("r", &HestonParams::r, "Risk-free rate")
        .def_rw("q", &HestonParams::q, "Dividend yield");

    // Heston 过程类 (StochasticProcess 子类)
    nb::class_<Heston>(m, "HestonProcess",
        "Heston stochastic volatility process: dV = κ(θ-V)dt + σ√V dW₂, "
        "dS/S = (r-q)dt + √V dW₁,  dW₁dW₂ = ρ dt. "
        "Supports Euler / FullTruncation / Exact schemes (QE_M via HestonQE subclass).")
        .def(nb::init<HestonParams, HestonScheme>(),
             nb::arg("p"),
             nb::arg("scheme") = HestonScheme::FullTruncation,
             "Construct Heston process with params and discretization scheme")
        .def("dimension", &Heston::dimension, "State dimension (always 2)")
        .def("spot", &Heston::spot, "Initial spot price")
        .def("characteristic_function",
             [](const Heston& self, std::complex<double> u, double tau) {
                 return self.characteristic_function(u, tau);
             },
             nb::arg("u"), nb::arg("tau"),
             "Heston CF φ(u; τ) (Little Trap correction, Albrecher 2007)")
        .def("generate_path",
             [](const Heston& self, double T, cpphub::v1::Size n_steps,
                std::vector<double>& path, Philox4x64& rng) {
                 if (path.size() < n_steps + 1) path.resize(n_steps + 1);
                 std::span<double> sp(path.data(), path.size());
                 self.generate_path(T, n_steps, sp, rng);
                 return path;
             },
             nb::arg("T"), nb::arg("n_steps"),
             nb::arg("path").none(false), nb::arg("rng"),
             "Generate spot path; path buffer resized to n_steps+1, returns path")
        .def("params",
             [](const Heston& self) { return self.params(); },
             "Return a copy of HestonProcessParams");

    // HestonQE (Quadratic-Exponential with martingale correction, Andersen 2008)
    nb::class_<HestonQE, Heston>(m, "HestonQEProcess",
        "Heston with Quadratic-Exponential discretization (Andersen 2008). "
        "Subclass of HestonProcess; overrides generate_path to use QE-M scheme "
        "(psi <= 1.5: noncentral chi², psi > 1.5: double exponential).")
        .def(nb::init<HestonParams>(), nb::arg("p"),
             "Construct HestonQE process (uses QE_M scheme internally)");

    // Rough Heston (El Euch-Rosenbaum 2018, fractional CIR)
    nb::class_<RoughHestonParams>(m, "RoughHestonParams",
        "Rough Heston parameters: H ∈ (0, 0.5] (Hurst), κ, θ, σ, ρ, v0, S0, r, q. "
        "α = H + 0.5. H=0.5 degenerates to standard Heston.")
        .def(nb::init<>())
        .def_rw("H", &RoughHestonParams::H, "Hurst exponent (0, 0.5]")
        .def_rw("kappa", &RoughHestonParams::kappa, "Mean reversion speed")
        .def_rw("theta", &RoughHestonParams::theta, "Long-term variance")
        .def_rw("sigma", &RoughHestonParams::sigma, "Vol of vol")
        .def_rw("rho", &RoughHestonParams::rho, "Correlation")
        .def_rw("v0", &RoughHestonParams::v0, "Initial variance")
        .def_rw("S0", &RoughHestonParams::S0, "Initial spot")
        .def_rw("r", &RoughHestonParams::r, "Risk-free rate")
        .def_rw("q", &RoughHestonParams::q, "Dividend yield");

    m.def("validate_rough_heston_params",
          [](const RoughHestonParams& p) { validate_rough_heston_params(p); },
          nb::arg("params"),
          "Validate Rough Heston params (throws on invalid)");

    m.def("rough_heston_kernel",
          [](double T, cpphub::v1::Size n_steps, double alpha) {
              return rough_heston_kernel(T, n_steps, alpha);
          },
          nb::arg("T"), nb::arg("n_steps"), nb::arg("alpha"),
          "Precompute Volterra kernel K[j][i] = [(t_{j+1}-t_i)^α - (t_{j+1}-t_{i+1})^α] / α. "
          "Returns N×N lower-triangular matrix.");

    nb::class_<RoughHestonProcess>(m, "RoughHestonProcess",
        "Rough Heston process (El Euch-Rosenbaum 2018). Fractional CIR variance + "
        "GBM price with leverage. O(N²) per step due to Volterra kernel convolution.")
        .def(nb::init<RoughHestonParams, double, cpphub::v1::Size>(),
             nb::arg("params"), nb::arg("T"), nb::arg("n_steps"),
             "Construct RoughHestonProcess; precomputes Volterra kernel for (T, n_steps)")
        .def("dimension", &RoughHestonProcess::dimension, "State dimension (always 2)")
        .def("spot", &RoughHestonProcess::spot, "Initial spot")
        .def("generate_path",
             [](const RoughHestonProcess& self, double T, cpphub::v1::Size n_steps,
                std::vector<double>& path, Philox4x64& rng) {
                 if (path.size() < n_steps + 1) path.resize(n_steps + 1);
                 std::span<double> sp(path.data(), path.size());
                 self.generate_path(T, n_steps, sp, rng);
                 return path;
             },
             nb::arg("T"), nb::arg("n_steps"),
             nb::arg("path").none(false), nb::arg("rng"),
             "Generate spot path; returns path buffer (resized to n_steps+1)")
        .def("generate_path_sv",
             [](const RoughHestonProcess& self, Philox4x64& rng, double dt_scale = 1.0) {
                 return self.generate_path(rng, dt_scale);
             },
             nb::arg("rng"), nb::arg("dt_scale") = 1.0,
             "Generate both spot and variance paths. "
             "Returns [S_path, v_path], each of length n_steps+1.")
        .def("params",
             [](const RoughHestonProcess& self) { return self.params(); },
             "Return a copy of RoughHestonParams")
        .def("T", &RoughHestonProcess::T, "Configured maturity")
        .def("n_steps", &RoughHestonProcess::n_steps, "Configured step count")
        .def("alpha", &RoughHestonProcess::alpha, "α = H + 0.5");

    // Rough Bergomi (Bayer-Friz-Gatheral 2016)
    nb::class_<RoughBergomiParams>(m, "RoughBergomiParams",
        "Rough Bergomi parameters: H ∈ (0, 0.5), η (vol-of-vol), ρ, ξ₀ (initial "
        "forward variance), S0, r, q. Uses RL-fBm (non-stationary).")
        .def(nb::init<>())
        .def_rw("H", &RoughBergomiParams::H, "Hurst exponent (0, 0.5)")
        .def_rw("eta", &RoughBergomiParams::eta, "Vol-of-vol (>0)")
        .def_rw("rho", &RoughBergomiParams::rho, "Correlation [-1, 1]")
        .def_rw("xi0", &RoughBergomiParams::xi0, "Initial forward variance (>0)")
        .def_rw("S0", &RoughBergomiParams::S0, "Initial spot (>0)")
        .def_rw("r", &RoughBergomiParams::r, "Risk-free rate")
        .def_rw("q", &RoughBergomiParams::q, "Dividend yield");

    m.def("validate_rough_bergomi_params",
          [](const RoughBergomiParams& p) { validate_rough_bergomi_params(p); },
          nb::arg("params"),
          "Validate Rough Bergomi params (throws on invalid)");

    // RLFbmSampler (Riemann-Liouville fractional Brownian motion sampler)
    nb::class_<RLFbmSampler>(m, "RLFbmSampler",
        "Riemann-Liouville fBm sampler: W̃^H_t = √(2H+1) ∫(t-s)^{H-1/2} dW_s. "
        "Uses exact Cholesky of covariance matrix (O(N³) setup, O(N²) per sample).")
        .def(nb::init<double, cpphub::v1::Size, double>(),
             nb::arg("T"), nb::arg("n_steps"), nb::arg("H"),
             "Construct sampler for grid [0, T] with n_steps intervals")
        .def("sample",
             [](const RLFbmSampler& self, const std::vector<double>& Z) {
                 return self.sample(Z);
             },
             nb::arg("Z"),
             "Sample W̃^H path from independent N(0,1) vector Z (length n_steps). "
             "Returns W̃^H_{t_{i+1}} for i=0..n_steps-1.")
        .def("H", &RLFbmSampler::H, "Configured Hurst exponent")
        .def("T", &RLFbmSampler::T, "Configured maturity")
        .def("n_steps", &RLFbmSampler::n_steps, "Configured step count")
        .def_static("log_v_mean", &RLFbmSampler::log_v_mean,
                    nb::arg("t"), nb::arg("xi0"), nb::arg("eta"), nb::arg("H"),
                    "Theoretical E[log v_t] = log(ξ₀) - 0.5·η²·t^{2H}")
        .def_static("log_v_var", &RLFbmSampler::log_v_var,
                    nb::arg("t"), nb::arg("eta"), nb::arg("H"),
                    "Theoretical Var(log v_t) = η²·t^{2H}");

    nb::class_<RoughBergomiProcess>(m, "RoughBergomiProcess",
        "Rough Bergomi process (Bayer-Friz-Gatheral 2016). "
        "v_t = ξ₀·exp(η·W̃^H_t - 0.5·η²·t^{2H}),  dS/S = √v·dZ. "
        "Non-Markovian; requires Cholesky of RL-fBm covariance.")
        .def(nb::init<RoughBergomiParams>(), nb::arg("params"),
             "Construct RoughBergomiProcess")
        .def("dimension", &RoughBergomiProcess::dimension, "State dimension (always 2)")
        .def("spot", &RoughBergomiProcess::spot, "Initial spot")
        .def("generate_path",
             [](const RoughBergomiProcess& self, double T, cpphub::v1::Size n_steps,
                std::vector<double>& path, Philox4x64& rng) {
                 if (path.size() < n_steps + 1) path.resize(n_steps + 1);
                 std::span<double> sp(path.data(), path.size());
                 self.generate_path(T, n_steps, sp, rng);
                 return path;
             },
             nb::arg("T"), nb::arg("n_steps"),
             nb::arg("path").none(false), nb::arg("rng"),
             "Generate spot path (recomputes Cholesky each call — use "
             "generate_path_with_sampler for repeated sampling)")
        .def("generate_path_with_sampler",
             [](const RoughBergomiProcess& self, double T, cpphub::v1::Size n_steps,
                std::vector<double>& path, Philox4x64& rng,
                const RLFbmSampler& sampler) {
                 if (path.size() < n_steps + 1) path.resize(n_steps + 1);
                 std::span<double> sp(path.data(), path.size());
                 self.generate_path_with_sampler(T, n_steps, sp, rng, sampler);
                 return path;
             },
             nb::arg("T"), nb::arg("n_steps"),
             nb::arg("path").none(false), nb::arg("rng"), nb::arg("sampler"),
             "Generate spot path using pre-built RLFbmSampler (avoids recomputing Cholesky)")
        .def("params",
             [](const RoughBergomiProcess& self) { return self.params(); },
             "Return a copy of RoughBergomiParams");

    // SABR 过程 (Hagan 2002; process simulation)
    // 注: calibrator 的 SABRParams (CalibSABRParams) 已在 §11 暴露为 "SABRParams"
    // 此处 SABRProcess 接受 CalibSABRParams 并内部转换为 sabr_hagan.hpp 的 SABRParams
    nb::enum_<SABRScheme>(m, "SABRScheme",
        "SABR process discretization scheme")
        .value("EulerAbsorbing", SABRScheme::EulerAbsorbing,
               "F Euler + absorbing wall, σ log-Euler (default, general beta)")
        .value("LogEuler", SABRScheme::LogEuler,
               "F and σ both log-Euler (exact for beta=1, approximate for beta<1)");

    nb::class_<SABRProcess>(m, "SABRProcess",
        "SABR stochastic process: dF = σ·F^β·dW₁, dσ = ν·σ·dW₂, dW₁dW₂ = ρ dt. "
        "Constructor accepts the same SABRParams used by the calibrator.")
        .def("__init__",
             [](nb::pointer_and_handle<SABRProcess> v,
                const CalibSABRParams& p, double F0,
                double r, double q, SABRScheme scheme) {
                 // CalibSABRParams {alpha, beta, nu, rho} → SABRParams {alpha, beta, rho, nu}
                 SABRParams hp{p.alpha, p.beta, p.rho, p.nu};
                 new (v.p) SABRProcess(hp, F0, r, q, scheme);
             },
             nb::arg("params"), nb::arg("F0"),
             nb::arg("r") = 0.0, nb::arg("q") = 0.0,
             nb::arg("scheme") = SABRScheme::EulerAbsorbing,
             "Construct SABRProcess from SABRParams + forward F0")
        .def("dimension", &SABRProcess::dimension, "State dimension (always 2)")
        .def("spot", &SABRProcess::spot, "Initial forward")
        .def("generate_path",
             [](const SABRProcess& self, double T, cpphub::v1::Size n_steps,
                std::vector<double>& path, Philox4x64& rng) {
                 if (path.size() < n_steps + 1) path.resize(n_steps + 1);
                 std::span<double> sp(path.data(), path.size());
                 self.generate_path(T, n_steps, sp, rng);
                 return path;
             },
             nb::arg("T"), nb::arg("n_steps"),
             nb::arg("path").none(false), nb::arg("rng"),
             "Generate forward path; returns path buffer (resized to n_steps+1)")
        .def("evolve",
             [](const SABRProcess& self, double F, double sigma,
                double dt, double Z1, double Z2) {
                 double sigma_mut = sigma;
                 double F_new = self.evolve(F, sigma_mut, dt, Z1, Z2);
                 return nb::make_tuple(F_new, sigma_mut);
             },
             nb::arg("F"), nb::arg("sigma"), nb::arg("dt"),
             nb::arg("Z1"), nb::arg("Z2"),
             "Single-step evolve; returns (F_new, sigma_new)")
        .def("forward0", &SABRProcess::forward0, "Initial forward")
        .def("r", &SABRProcess::r, "Risk-free rate")
        .def("q", &SABRProcess::q, "Dividend yield");

    // =========================================================================
    // 19. v1.4 扩展: Risk 模块 (Expected Shortfall + MC VaR + Backtesting)
    // =========================================================================
    nb::class_<ExpectedShortfall>(m, "ExpectedShortfall",
        "Expected Shortfall (ES / CVaR) calculator. "
        "Supports empirical, Normal, Student-t, Cornish-Fisher expansions.")
        .def(nb::init<>())
        .def("from_losses", &ExpectedShortfall::from_losses,
             nb::arg("losses"), nb::arg("confidence") = 0.99,
             "Empirical ES from loss vector (sorted tail mean)")
        .def("from_mc_paths", &ExpectedShortfall::from_mc_paths,
             nb::arg("pnl_paths"), nb::arg("confidence") = 0.99,
             "ES from MC P&L paths (alias of from_losses)")
        .def("normal_es", &ExpectedShortfall::normal_es,
             nb::arg("mean"), nb::arg("sigma"), nb::arg("confidence") = 0.99,
             "Closed-form Normal ES = -μ + σ·φ(z)/(1-c)")
        .def("student_t_es", &ExpectedShortfall::student_t_es,
             nb::arg("mean"), nb::arg("sigma"), nb::arg("dof"),
             nb::arg("confidence") = 0.99,
             "Student-t ES (Newton iteration for t-quantile)")
        .def("cornish_fisher_es", &ExpectedShortfall::cornish_fisher_es,
             nb::arg("mean"), nb::arg("sigma"), nb::arg("skew"), nb::arg("kurt"),
             nb::arg("confidence") = 0.99,
             "Cornish-Fisher ES (skew/kurtosis-adjusted quantile)")
        .def("tail_average", &ExpectedShortfall::tail_average,
             nb::arg("losses"), nb::arg("confidence"),
             nb::arg("n_tail_points") = 100,
             "Tail average over n_tail_points worst losses");

    // MCVarConfig
    nb::class_<MCVarConfig>(m, "MCVarConfig",
        "Monte Carlo VaR configuration")
        .def(nb::init<>())
        .def_rw("n_paths", &MCVarConfig::n_paths, "Number of paths (default 100000)")
        .def_rw("seed", &MCVarConfig::seed, "RNG seed (default 42)")
        .def_rw("antithetic", &MCVarConfig::antithetic,
                "Enable antithetic variates (default True)");

    // VaRApproximation enum
    nb::enum_<VaRApproximation>(m, "VaRApproximation",
        "MC VaR approximation level")
        .value("Full", VaRApproximation::Full,
               "Full revaluation (default, slowest but exact)")
        .value("DeltaGamma", VaRApproximation::DeltaGamma,
               "Delta-Gamma approximation (2nd-order Taylor)")
        .value("Delta", VaRApproximation::Delta,
               "Delta-only approximation (1st-order Taylor)");

    // MCVaR class
    // portfolio_value_fn 接受 risk factor vector, 返回组合价值
    nb::class_<MCVaR>(m, "MCVaR",
        "Monte Carlo VaR engine for multi-factor portfolios. "
        "Constructor takes a Python callable portfolio_value_fn(factors) -> value, "
        "current risk factor levels, and covariance matrix (flattened, row-major).")
        .def(nb::init<std::function<Real(const std::vector<Real>&)>,
                      std::vector<Real>, std::vector<Real>,
                      cpphub::v1::Size, MCVarConfig>(),
             nb::arg("portfolio_value_fn"),
             nb::arg("current_risk_factors"),
             nb::arg("covariance_matrix"),
             nb::arg("n_factors"),
             nb::arg("config") = MCVarConfig{},
             "Construct MCVaR: portfolio_value_fn(factors) -> value, "
             "covariance matrix is flattened row-major (n_factors² entries)")
        .def("simulate_pnl_full", &MCVaR::simulate_pnl_full,
             "Simulate P&L via full revaluation; returns vector of P&L")
        .def("simulate_pnl_delta_gamma",
             [](const MCVaR& self, const std::vector<Real>& delta,
                const std::vector<Real>& gamma) {
                 return self.simulate_pnl_delta_gamma(delta, gamma);
             },
             nb::arg("delta"), nb::arg("gamma"),
             "Simulate P&L via Delta-Gamma approximation; "
             "gamma is flattened n×n Hessian (row-major)")
        .def("simulate_pnl_delta", &MCVaR::simulate_pnl_delta,
             nb::arg("delta"),
             "Simulate P&L via Delta-only approximation")
        .def("var",
             [](MCVaR& self, double confidence, VaRApproximation approx,
                const std::vector<Real>& delta, const std::vector<Real>& gamma) {
                 return self.var(confidence, approx, delta, gamma);
             },
             nb::arg("confidence") = 0.99,
             nb::arg("approx") = VaRApproximation::Full,
             nb::arg("delta") = std::vector<Real>{},
             nb::arg("gamma") = std::vector<Real>{},
             "Compute VaR at given confidence level")
        .def("standard_error",
             [](const MCVaR& self, double confidence, cpphub::v1::Size n_bootstrap) {
                 return self.standard_error(confidence, n_bootstrap);
             },
             nb::arg("confidence") = 0.99,
             nb::arg("n_bootstrap") = 1000,
             "Bootstrap standard error of VaR estimate");

    // Backtesting: BacktestResult
    nb::class_<BacktestResult>(m, "BacktestResult",
        "VaR backtest result")
        .def_ro("n_violations", &BacktestResult::n_violations,
                "Number of VaR violations")
        .def_ro("n_observations", &BacktestResult::n_observations,
                "Number of observations")
        .def_ro("violation_rate", &BacktestResult::violation_rate,
                "Observed violation rate")
        .def_ro("expected_rate", &BacktestResult::expected_rate,
                "Expected violation rate (1 - confidence)")
        .def_ro("p_value", &BacktestResult::p_value,
                "Test p-value")
        .def_ro("reject_null", &BacktestResult::reject_null,
                "True if null rejected at 5% level")
        .def("__repr__", [](const BacktestResult& r) {
            return "<BacktestResult violations=" + std::to_string(r.n_violations) +
                   "/" + std::to_string(r.n_observations) +
                   " rate=" + std::to_string(r.violation_rate) +
                   " expected=" + std::to_string(r.expected_rate) +
                   " p=" + std::to_string(r.p_value) +
                   " reject=" + (r.reject_null ? "True" : "False") + ">";
        });

    // Kupiec POF (Proportion of Failures) test
    nb::class_<KupiecPOF>(m, "KupiecPOF",
        "Kupiec Proportion-of-Failures (POF) test. "
        "Likelihood ratio test of unconditional coverage (H₀: violation rate = 1-c).")
        .def_static("test",
                    [](Size n_violations, Size n_observations, double confidence) {
                        return KupiecPOF::test(n_violations, n_observations, confidence);
                    },
                    nb::arg("n_violations"), nb::arg("n_observations"),
                    nb::arg("confidence") = 0.99,
                    "Run Kupiec POF test from violation/observation counts")
        .def_static("test_series",
                    [](const std::vector<Real>& var_series,
                       const std::vector<Real>& realized_losses, double confidence) {
                        return KupiecPOF::test(var_series, realized_losses, confidence);
                    },
                    nb::arg("var_series"), nb::arg("realized_losses"),
                    nb::arg("confidence") = 0.99,
                    "Run Kupiec POF test from VaR series and realized losses");

    // Christoffersen IID test (independence of violations)
    nb::class_<ChristoffersenIID>(m, "ChristoffersenIID",
        "Christoffersen IID test for independence of VaR violations. "
        "Tests whether breaches cluster (H₀: breaches are IID).")
        .def_static("test",
                    [](const std::vector<Real>& var_series,
                       const std::vector<Real>& realized_losses) {
                        return ChristoffersenIID::test(var_series, realized_losses);
                    },
                    nb::arg("var_series"), nb::arg("realized_losses"),
                    "Run Christoffersen independence test")
        .def_static("joint_test",
                    [](const std::vector<Real>& var_series,
                       const std::vector<Real>& realized_losses, double confidence) {
                        return ChristoffersenIID::joint_test(var_series, realized_losses, confidence);
                    },
                    nb::arg("var_series"), nb::arg("realized_losses"),
                    nb::arg("confidence") = 0.99,
                    "Joint conditional coverage test (LR_cc = LR_pof + LR_ind)");

    // Basel traffic light
    nb::enum_<BaselZone>(m, "BaselZone",
        "Basel III traffic light zone")
        .value("Green", BaselZone::Green, "≤4 violations: model OK, multiplier 3.0")
        .value("Yellow", BaselZone::Yellow, "5-9 violations: model may be inaccurate")
        .value("Red", BaselZone::Red, "≥10 violations: model inaccurate, multiplier 4.0");

    nb::class_<BaselTrafficLightResult>(m, "BaselTrafficLightResult",
        "Basel traffic light assessment result")
        .def_ro("zone", &BaselTrafficLightResult::zone, "Traffic light zone")
        .def_ro("n_violations", &BaselTrafficLightResult::n_violations,
                "Number of violations")
        .def_ro("n_observations", &BaselTrafficLightResult::n_observations,
                "Number of observations (default 250)")
        .def_ro("capital_multiplier", &BaselTrafficLightResult::capital_multiplier,
                "Capital multiplier to apply")
        .def_ro("description", &BaselTrafficLightResult::description,
                "Human-readable description")
        .def("__repr__", [](const BaselTrafficLightResult& r) {
            return "<BaselTrafficLight zone=" + std::to_string(static_cast<int>(r.zone)) +
                   " violations=" + std::to_string(r.n_violations) +
                   " multiplier=" + std::to_string(r.capital_multiplier) + ">";
        });

    nb::class_<BaselTrafficLight>(m, "BaselTrafficLight",
        "Basel III traffic light test for 250-day VaR backtesting")
        .def_static("assess",
                    [](Size n_violations, Size n_observations) {
                        return BaselTrafficLight::assess(n_violations, n_observations);
                    },
                    nb::arg("n_violations"),
                    nb::arg("n_observations") = 250,
                    "Assess Basel zone from violation count");

    // =========================================================================
    // 20. v1.4 扩展: GreeksFactory (统一 Greeks 入口, 方法自动分派)
    // =========================================================================
    nb::enum_<GreeksMethod>(m, "GreeksMethod",
        "Greeks computation method (used by GreeksFactory)")
        .value("Auto", GreeksMethod::Auto,
               "Auto-select: vanilla→Analytic, digital→LR")
        .value("Analytic", GreeksMethod::Analytic,
               "BSM closed-form (fastest, exact; vanilla European only)")
        .value("Pathwise", GreeksMethod::Pathwise,
               "Pathwise method (smooth payoff, lowest variance)")
        .value("LR", GreeksMethod::LR,
               "Likelihood ratio (discontinuous payoff: digital/barrier)")
        .value("FD", GreeksMethod::FD,
               "Finite difference (universal fallback, centered diff)")
        .value("AAD", GreeksMethod::AAD,
               "Automatic differentiation (complex/basket; exact but slower)");

    nb::enum_<PayoffType>(m, "PayoffType",
        "Payoff type for GreeksFactory dispatch")
        .value("VanillaCall", PayoffType::VanillaCall, "max(S_T - K, 0)")
        .value("VanillaPut", PayoffType::VanillaPut, "max(K - S_T, 0)")
        .value("DigitalCall", PayoffType::DigitalCall, "1{S_T > K}")
        .value("DigitalPut", PayoffType::DigitalPut, "1{S_T < K}");

    nb::class_<UnifiedGreeks>(m, "UnifiedGreeks",
        "Unified Greeks result from GreeksFactory")
        .def_ro("price", &UnifiedGreeks::price, "Option price")
        .def_ro("delta", &UnifiedGreeks::delta, "dPrice/dS")
        .def_ro("gamma", &UnifiedGreeks::gamma, "d²Price/dS²")
        .def_ro("vega", &UnifiedGreeks::vega, "dPrice/dσ")
        .def_ro("theta", &UnifiedGreeks::theta, "-dPrice/dT")
        .def_ro("rho", &UnifiedGreeks::rho, "dPrice/dr")
        .def_ro("method_used", &UnifiedGreeks::method_used,
                "Method actually used (Auto resolves to one of the others)")
        .def_ro("note", &UnifiedGreeks::note,
                "Diagnostic note (e.g., fallback reason)")
        .def("__repr__", [](const UnifiedGreeks& g) {
            return "<UnifiedGreeks price=" + std::to_string(g.price) +
                   " delta=" + std::to_string(g.delta) +
                   " gamma=" + std::to_string(g.gamma) +
                   " vega=" + std::to_string(g.vega) +
                   " method=" + std::to_string(static_cast<int>(g.method_used)) + ">";
        });

    // GreeksFactory: static class, expose compute_bsm as static method
    nb::class_<GreeksFactory>(m, "GreeksFactory",
        "Unified Greeks entry point with automatic method dispatch. "
        "Auto: vanilla European → Analytic; discontinuous (digital) → LR.")
        .def_static("compute_bsm",
                    [](double S, double K, double T, double r, double q,
                       double sigma, PayoffType payoff,
                       GreeksMethod method, cpphub::v1::Size n_paths, uint64_t seed) {
                        return GreeksFactory::compute_bsm(
                            S, K, T, r, q, sigma, payoff, method, n_paths, seed);
                    },
                    nb::arg("S"), nb::arg("K"), nb::arg("T"),
                    nb::arg("r"), nb::arg("q"), nb::arg("sigma"),
                    nb::arg("payoff"),
                    nb::arg("method") = GreeksMethod::Auto,
                    nb::arg("n_paths") = 100000,
                    nb::arg("seed") = 42,
                    "Compute unified Greeks for BSM European option. "
                    "Auto dispatch: vanilla→Analytic, digital→LR.");

    // =========================================================================
    // 21. 版本信息
    // =========================================================================
    m.attr("__version__") = "1.4.0";
    m.attr("__author__") = "Scott (鹏)";
}
