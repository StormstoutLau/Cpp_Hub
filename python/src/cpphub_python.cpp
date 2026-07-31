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

// ---------------------------------------------------------------------------
// 重要: calibrator.hpp 在 cpphub::v1 中重新定义了 HestonParams 和 SABRParams,
// 与 heston.hpp (8 字段) 和 sabr_hagan.hpp (字段顺序不同) 冲突。
// 解决方案: 通过宏重命名 calibrator.hpp 内部的两个 struct, 避免冲突。
// Python 侧仍然以 HestonParams / SABRParams 名称暴露。
// ---------------------------------------------------------------------------
#define HestonParams CalibHestonParams
#define SABRParams CalibSABRParams
#include "cpphub/calibration/calibrator.hpp"
#undef HestonParams
#undef SABRParams

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
    // 13. 版本信息
    // =========================================================================
    m.attr("__version__") = "1.2.0";
    m.attr("__author__") = "Scott (鹏)";
}
