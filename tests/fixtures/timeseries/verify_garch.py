# verify_garch.py — GARCH(1,1) 基准生成 (Phase 7B v1.6 M1, spec §3.1)
#
# 生成 tests/unit/timeseries/garch_baseline.inc:
#   - 模拟数据 T=1000 (numpy default_rng(42), ω=0.05, α=0.10, β=0.85, μ=0.05)
#   - arch mean='Zero' Normal: params/llf/aic/bic/backcast/h_t/robust cov/SE/forecast(h=10)
#   - arch mean='Zero' dist='t': params(含 ν)/llf
#
# 对齐策略 (spec §2.0.2 G2): C++ 内部先去均值再估计, 故 Python 基准用
# 去均值数据 + mean='Zero' 拟合, 残差序列完全一致。
# 数据尺度 ~1 (无条件方差 = 0.05/0.05 = 1), 无 DataScaleWarning。
#
# 运行: python tests/fixtures/timeseries/verify_garch.py
import numpy as np
from arch import arch_model

OUT = r"f:\Cpp_Hub\tests\unit\timeseries\garch_baseline.inc"

# --- 1. 模拟数据 (确定性) ---
rng = np.random.default_rng(42)
T = 1000
w, a, b, mu = 0.05, 0.10, 0.85, 0.05
h = np.empty(T)
eps = np.empty(T)
h_prev = w / (1.0 - a - b)  # 无条件方差起步
for t in range(T):
    if t == 0:
        h[t] = h_prev
    else:
        h[t] = w + a * eps[t - 1] ** 2 + b * h_prev
    eps[t] = np.sqrt(h[t]) * rng.standard_normal()
    h_prev = h[t]
y = mu + eps  # 含非零均值, C++ 内部去均值 (G2)

y_dm = y - np.mean(y)  # 与 C++ G2 一致

# --- 2. arch Normal 基准 ---
am = arch_model(y_dm, mean='Zero', vol='GARCH', p=1, q=1, dist='normal')
res = am.fit(disp='off', show_warning=False)
p = res.params
omega, alpha, beta = float(p['omega']), float(p['alpha[1]']), float(p['beta[1]'])
llf, aic, bic = float(res.loglikelihood), float(res.aic), float(res.bic)
bc = float(am.volatility.backcast(y_dm))
h_arch = res.conditional_volatility ** 2
cov = np.asarray(res.param_cov, dtype=float)  # robust (默认 cov_type)
se = np.sqrt(np.diag(cov))

# forecast: analytic h_{T+1..T+10}
fc = res.forecast(horizon=10, method='analytic', reindex=False)
fc_var = np.asarray(fc.variance.values[0], dtype=float)

print("Normal: omega=%.12g alpha=%.12g beta=%.12g" % (omega, alpha, beta))
print("llf=%.12g aic=%.12g bic=%.12g backcast=%.12g" % (llf, aic, bic, bc))

# --- 3. arch t 分布基准 (ν 联合估计) ---
# Normal 数据下 ν 处于平坦方向 (ν→∞), 跨优化器不可对齐;
# 用 t(6) 冲击模拟数据 (seed 43) 专测 ν 联合估计
rng2 = np.random.default_rng(43)
nu_true = 6.0
h2 = np.empty(T)
eps2 = np.empty(T)
h_prev = w / (1.0 - a - b)
for t in range(T):
    if t == 0:
        h2[t] = h_prev
    else:
        h2[t] = w + a * eps2[t - 1] ** 2 + b * h_prev
    z = rng2.standard_t(nu_true) / np.sqrt(nu_true / (nu_true - 2.0))  # 标准化 t
    eps2[t] = np.sqrt(h2[t]) * z
    h_prev = h2[t]
y2 = mu + eps2
y2_dm = y2 - np.mean(y2)

am_t = arch_model(y2_dm, mean='Zero', vol='GARCH', p=1, q=1, dist='t')
res_t = am_t.fit(disp='off', show_warning=False)
pt = res_t.params
print("t: omega=%.12g alpha=%.12g beta=%.12g nu=%.12g llf=%.12g"
      % (pt['omega'], pt['alpha[1]'], pt['beta[1]'], pt['nu'], res_t.loglikelihood))


def fmt(v):
    return repr(float(v))


lines = []
lines.append("// 自动生成: tests/fixtures/timeseries/verify_garch.py (arch 8.0.0)")
lines.append("// 勿手改 — 重新生成请运行该脚本")
lines.append("#pragma once")
lines.append("namespace cpphub { inline namespace v1 { namespace timeseries {")
lines.append("namespace garch { namespace baseline {")
lines.append("")
lines.append("constexpr Size T = %d;" % T)
lines.append("")
lines.append("// 模拟数据 y_t = mu + eps_t (GARCH(1,1): w=0.05,a=0.10,b=0.85,mu=0.05, seed 42)")
lines.append("constexpr Real DATA[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(y[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// arch mean='Zero' vol='GARCH' dist='normal' 基准")
lines.append("constexpr Real MEAN_SHIFT = %s;  // mean(y), C++ G2 去均值用" % fmt(np.mean(y)))
lines.append("constexpr Real BACKCAST = %s;" % fmt(bc))
lines.append("constexpr Real OMEGA = %s;" % fmt(omega))
lines.append("constexpr Real ALPHA = %s;" % fmt(alpha))
lines.append("constexpr Real BETA = %s;" % fmt(beta))
lines.append("constexpr Real LLF = %s;" % fmt(llf))
lines.append("constexpr Real AIC = %s;" % fmt(aic))
lines.append("constexpr Real BIC = %s;" % fmt(bic))
lines.append("constexpr Real VCOV[3][3] = {")
for i in range(3):
    lines.append("    {%s}," % ", ".join(fmt(cov[i][j]) for j in range(3)))
lines.append("};")
lines.append("constexpr Real SE[3] = {%s};"
             % ", ".join(fmt(se[j]) for j in range(3)))
lines.append("")
lines.append("// arch conditional variance h_t (T=1000, 递归对照用)")
lines.append("constexpr Real H[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(h_arch[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// arch forecast analytic h_{T+1..T+10}")
lines.append("constexpr Real FC10[10] = {%s};"
             % ", ".join(fmt(v) for v in fc_var))
lines.append("")
lines.append("// t(6) 冲击模拟数据 (seed 43, GARCH 同参数, mu=0.05)")
lines.append("constexpr Real TDATA[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(y2[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// arch mean='Zero' vol='GARCH' dist='t' 基准 (ν 联合估计, t(6) 数据)")
lines.append("constexpr Real T_NU = %s;" % fmt(pt['nu']))
lines.append("constexpr Real T_OMEGA = %s;" % fmt(pt['omega']))
lines.append("constexpr Real T_ALPHA = %s;" % fmt(pt['alpha[1]']))
lines.append("constexpr Real T_BETA = %s;" % fmt(pt['beta[1]']))
lines.append("constexpr Real T_LLF = %s;" % fmt(res_t.loglikelihood))
lines.append("")
lines.append("}}}}}  // baseline::garch::timeseries::v1::cpphub")

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
print("written:", OUT)
