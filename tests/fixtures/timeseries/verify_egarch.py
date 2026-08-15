# verify_egarch.py — EGARCH(1,1) 基准生成 (Phase 7B v1.6 M1, spec §3.1)
#
# 生成 tests/unit/timeseries/egarch_baseline.inc:
#   - 模拟数据 T=1000 (numpy default_rng(52), ω=0, α=-0.05, β=0.94, γ=0.15, μ=0.05)
#     含杠杆效应 (α<0: 负冲击放大波动, Nelson 1991 θ)
#   - arch mean='Zero' vol='EGARCH' p=1,o=1,q=1 Normal:
#     params/llf/aic/bic/backcast/h_t/robust cov/SE/1-step analytic forecast
#   - arch mean='Zero' dist='t': params(含 ν)/llf (t(6) 冲击数据, seed 53)
#   - simulation 多步预测的中位数路径 (median(h) ≈ exp(E[ln h]), Jensen 差
#     由 C++ 测试端用宽松容差处理)
#
# 参数映射 (G23, spec §2.0.3):
#   arch alpha[1] (对称 |z|-E|z| 项) = spec gamma
#   arch gamma[1] (非对称 z 项)      = spec alpha
#   spec 参数顺序: {omega, alpha(非对称), beta, gamma(对称), [nu]}
#
# 运行: python tests/fixtures/timeseries/verify_egarch.py
import numpy as np
from arch import arch_model

OUT = r"f:\Cpp_Hub\tests\unit\timeseries\egarch_baseline.inc"

E_ABS_Z = np.sqrt(2.0 / np.pi)

# --- 1. 模拟数据 (确定性, 含杠杆) ---
rng = np.random.default_rng(52)
T = 1000
# spec 约定: alpha = θ (非对称 z 项, 负值 → 杠杆), gamma = Nelson γ (对称项)
w, al, be, ga, mu = 0.0, -0.05, 0.94, 0.15, 0.05
ln_h = np.empty(T)
eps = np.empty(T)
ln_h_prev = w / (1.0 - be)  # = 0 → h₀ = 1 (中位数方差 1)
z_prev = 0.0
for t in range(T):
    if t > 0:
        ln_h_prev = (w + be * ln_h_prev + al * z_prev
                     + ga * (abs(z_prev) - E_ABS_Z))
    ln_h[t] = ln_h_prev
    sig = np.sqrt(np.exp(ln_h[t]))
    eps[t] = sig * rng.standard_normal()
    z_prev = eps[t] / sig
y = mu + eps  # 含非零均值, C++ 内部去均值 (G2)
y_dm = y - np.mean(y)

# --- 2. arch Normal EGARCH 基准 ---
am = arch_model(y_dm, mean='Zero', vol='EGARCH', p=1, o=1, q=1, dist='normal')
res = am.fit(disp='off', show_warning=False)
p = res.params
omega = float(p['omega'])
arch_alpha_sym = float(p['alpha[1]'])   # 对称项 = spec gamma
arch_gamma_asym = float(p['gamma[1]'])  # 非对称项 = spec alpha
beta = float(p['beta[1]'])
llf, aic, bic = float(res.loglikelihood), float(res.aic), float(res.bic)
bc = float(am.volatility.backcast(y_dm))
h_arch = np.asarray(res.conditional_volatility, dtype=float) ** 2
cov = np.asarray(res.param_cov, dtype=float)  # arch 顺序 [ω, α_sym, γ_asym, β]
se_arch = np.sqrt(np.diag(cov))

print("Normal: omega=%.12g alpha(asym)=%.12g beta=%.12g gamma(sym)=%.12g"
      % (omega, arch_gamma_asym, beta, arch_alpha_sym))
print("llf=%.12g aic=%.12g bic=%.12g backcast=%.12g" % (llf, aic, bic, bc))

# 参数重排 → spec 顺序 [ω, alpha(=γ_arch), beta, gamma(=α_arch)]
P = [0, 2, 3, 1]
cov_spec = cov[np.ix_(P, P)]
se_spec = np.sqrt(np.diag(cov_spec))

# 1-step analytic forecast + 末态量 (C++ forecast_egarch 输入)
fc1 = res.forecast(horizon=1, method='analytic', reindex=False)
fc1_var = float(np.asarray(fc1.variance.values[0], dtype=float)[0])
ln_h_T = float(np.log(h_arch[-1]))
z_T = float(y_dm[-1] / np.sqrt(h_arch[-1]))

# simulation 多步预测 (spec §2.0.5 注: arch EGARCH 多步仅 simulation)
# 对照量: fc.variance = 模拟路径 h 的均值 E[h] (Jensen: > exp(E ln h))。
# 注 (probe_egarch_sim.py 2026-08-15): fc.simulations.variances 内部存储
# 异常 (mean ≡ median 逐位相等, 非单调), 不可用作路径对照; fc.variance 经
# 矩母函数闭式验证: E[e^{tξ}] = e^{-tγc}[e^{t²(α-γ)²/2}Φ(-t(α-γ))
#                               + e^{t²(α+γ)²/2}Φ(t(α+γ))], ξ = αz+γ(|z|-c)
# 手算 E[h_{T+2}] = 0.6501 vs arch 0.6477 (MC 噪声 0.4%) ✓
fc_sim = res.forecast(horizon=10, method='simulation', reindex=False)
fc_sim_var = np.asarray(fc_sim.variance.values[0], dtype=float)

# --- 3. arch t 分布基准 (ν 联合估计, t(6) 冲击数据 seed 53) ---
rng2 = np.random.default_rng(53)
nu_true = 6.0
ln_h2 = np.empty(T)
eps2 = np.empty(T)
ln_h_prev = w / (1.0 - be)
z_prev = 0.0
for t in range(T):
    if t > 0:
        ln_h_prev = (w + be * ln_h_prev + al * z_prev
                     + ga * (abs(z_prev) - E_ABS_Z))
    ln_h2[t] = ln_h_prev
    sig = np.sqrt(np.exp(ln_h2[t]))
    z = rng2.standard_t(nu_true) / np.sqrt(nu_true / (nu_true - 2.0))
    eps2[t] = sig * z
    z_prev = z
y2 = mu + eps2
y2_dm = y2 - np.mean(y2)

am_t = arch_model(y2_dm, mean='Zero', vol='EGARCH', p=1, o=1, q=1, dist='t')
res_t = am_t.fit(disp='off', show_warning=False)
pt = res_t.params
print("t: omega=%.12g alpha(asym)=%.12g beta=%.12g gamma(sym)=%.12g nu=%.12g llf=%.12g"
      % (pt['omega'], pt['gamma[1]'], pt['beta[1]'], pt['alpha[1]'], pt['nu'],
         res_t.loglikelihood))


def fmt(v):
    return repr(float(v))


lines = []
lines.append("// 自动生成: tests/fixtures/timeseries/verify_egarch.py (arch 8.0.0)")
lines.append("// 勿手改 — 重新生成请运行该脚本")
lines.append("// 参数顺序 (spec §2.0.3 G23): {omega, alpha(非对称 z 项),")
lines.append("//   beta, gamma(对称 |z|-E|z| 项), [nu]} — 已从 arch 顺序映射")
lines.append("#pragma once")
lines.append("namespace cpphub { inline namespace v1 { namespace timeseries {")
lines.append("namespace garch { namespace egarch_baseline {")
lines.append("")
lines.append("constexpr Size T = %d;" % T)
lines.append("constexpr Real E_ABS_Z = %s;  // sqrt(2/pi)" % fmt(E_ABS_Z))
lines.append("")
lines.append("// 模拟数据 y_t = mu + eps_t (EGARCH: w=0,al=-0.05,be=0.94,ga=0.15,mu=0.05, seed 52)")
lines.append("constexpr Real DATA[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(y[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// arch mean='Zero' vol='EGARCH' p=1,o=1,q=1 dist='normal' 基准 (spec 顺序)")
lines.append("constexpr Real BACKCAST = %s;" % fmt(bc))
lines.append("constexpr Real OMEGA = %s;" % fmt(omega))
lines.append("constexpr Real ALPHA = %s;  // 非对称 (arch gamma[1])" % fmt(arch_gamma_asym))
lines.append("constexpr Real BETA = %s;" % fmt(beta))
lines.append("constexpr Real GAMMA = %s;  // 对称 (arch alpha[1])" % fmt(arch_alpha_sym))
lines.append("constexpr Real LLF = %s;" % fmt(llf))
lines.append("constexpr Real AIC = %s;" % fmt(aic))
lines.append("constexpr Real BIC = %s;" % fmt(bic))
lines.append("constexpr Real VCOV[4][4] = {")
for i in range(4):
    lines.append("    {%s}," % ", ".join(fmt(cov_spec[i][j]) for j in range(4)))
lines.append("};")
lines.append("constexpr Real SE[4] = {%s};"
             % ", ".join(fmt(se_spec[j]) for j in range(4)))
lines.append("")
lines.append("// arch conditional variance h_t (T=1000, 递归对照用)")
lines.append("constexpr Real H[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(h_arch[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// forecast 输入末态 + 基准 (analytic 1-step / simulation 10-step E[h])")
lines.append("constexpr Real LN_H_T = %s;" % fmt(ln_h_T))
lines.append("constexpr Real Z_T = %s;" % fmt(z_T))
lines.append("constexpr Real FC1 = %s;" % fmt(fc1_var))
lines.append("// simulation 路径 h 均值 E[h_{T+k}] (C++ 用矩母函数闭式对照, 3%%)")
lines.append("constexpr Real FC_SIM_VAR[10] = {%s};"
             % ", ".join(fmt(v) for v in fc_sim_var))
lines.append("")
lines.append("// t(6) 冲击模拟数据 (seed 53, EGARCH 同参数, mu=0.05)")
lines.append("constexpr Real TDATA[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(y2[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// arch mean='Zero' vol='EGARCH' dist='t' 基准 (ν 联合估计, spec 顺序)")
lines.append("constexpr Real T_OMEGA = %s;" % fmt(pt['omega']))
lines.append("constexpr Real T_ALPHA = %s;  // 非对称 (arch gamma[1])" % fmt(pt['gamma[1]']))
lines.append("constexpr Real T_BETA = %s;" % fmt(pt['beta[1]']))
lines.append("constexpr Real T_GAMMA = %s;  // 对称 (arch alpha[1])" % fmt(pt['alpha[1]']))
lines.append("constexpr Real T_NU = %s;" % fmt(pt['nu']))
lines.append("constexpr Real T_LLF = %s;" % fmt(res_t.loglikelihood))
lines.append("")
lines.append("}}}}}  // egarch_baseline::garch::timeseries::v1::cpphub")

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
print("written:", OUT)
