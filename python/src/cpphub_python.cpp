// Phase 4 LITE - M2: Cpp_Hub Python 绑定 (nanobind)
//
// 覆盖核心模块:
//   - core: normal_pdf/cdf/inv_cdf
//   - bsm: BSM 欧式期权定价
//   - heston: Heston 特征函数
//   - greeks: Analytic + AAD Greeks
//   - var: Historical/Parametric/MC VaR
//
// 构建: pip install . (需 nanobind)
// 测试: pytest python/tests/

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/complex.h>

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/risk/greeks/aad_greeks.hpp"
#include "cpphub/risk/var/historical_var.hpp"
#include "cpphub/risk/var/parametric_var.hpp"

namespace nb = nanobind;
using namespace cpphub::v1;

NB_MODULE(_core, m) {
    m.doc() = "Cpp_Hub: Quantitative finance library (BSM/Heston/MC/Greeks/VaR)";

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
    // 6. 版本信息
    // =========================================================================
    m.attr("__version__") = "1.1.0";
    m.attr("__author__") = "Scott (鹏)";
}
