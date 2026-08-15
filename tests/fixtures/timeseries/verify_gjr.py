# verify_gjr.py — GJR-GARCH(1,1) 基准生成 (Phase 7B v1.6 M1, spec §3.1)
#
# 生成 tests/unit/timeseries/gjr_baseline.inc:
#   - 模拟数据 T=1000 (numpy default_rng(62), ω=0.05, α=0.05, γ=0.15, β=0.80,
#     μ=0.05; 含杠杆效应 γ>0)
#   - arch mean='Zero' vol='GARCH' p=1,o=1,q=1 Normal:
#     params/llf/aic/bic/backcast/h_t/robust cov/SE/forecast(analytic h=10)
#   - arch mean='Zero' dist='t': params(含 ν)/llf (t(6) 冲击数据, seed 63)
#
# 约定探测 (排幻觉点, 2026-08-15):
#   probe-1 h₁ 初始化: arch compute_variance t=0 对 o 项用 0.5*backcast,
#     即 h₁ = ω + (α + γ/2 + β)·σ²₀ (与 GARCH(1,1) 的 h₁=ω+(α+β)σ²₀ 同构)
#   probe-2 多步预测: E[I(z<0)·ε²] = h/2 (对称分布) → φ = α+γ/2+β;
#     对照 arch analytic forecast 确认 (排除零残差递归伪影 φ=α+β)
#
# 运行: python tests/fixtures/timeseries/verify_gjr.py
import numpy as np
from arch import arch_model

OUT = r"f:\Cpp_Hub\tests\unit\timeseries\gjr_baseline.inc"

# --- 1. 模拟数据 (确定性, 含杠杆 γ>0) ---
rng = np.random.default_rng(62)
T = 1000
w, al, ga, be, mu = 0.05, 0.05, 0.15, 0.80, 0.05
uc_var = w / (1.0 - al - ga / 2.0 - be)  # 无条件方差 ω/(1-α-γ/2-β)
h = np.empty(T)
eps = np.empty(T)
h_prev = uc_var
for t in range(T):
    if t == 0:
        h[t] = h_prev
    else:
        h[t] = w + al * eps[t - 1] ** 2 \
               + ga * (eps[t - 1] < 0) * eps[t - 1] ** 2 + be * h_prev
    eps[t] = np.sqrt(h[t]) * rng.standard_normal()
    h_prev = h[t]
y = mu + eps  # 含非零均值, C++ 内部去均值 (G2)
y_dm = y - np.mean(y)

# --- 2. arch Normal GJR 基准 ---
am = arch_model(y_dm, mean='Zero', vol='GARCH', p=1, o=1, q=1, dist='normal')
res = am.fit(disp='off', show_warning=False)
p = res.params
omega = float(p['omega'])
alpha = float(p['alpha[1]'])
gamma = float(p['gamma[1]'])
beta = float(p['beta[1]'])
llf, aic, bic = float(res.loglikelihood), float(res.aic), float(res.bic)
bc = float(am.volatility.backcast(y_dm))
h_arch = res.conditional_volatility ** 2
cov = np.asarray(res.param_cov, dtype=float)  # arch 顺序 [ω, α, γ, β] = spec 顺序
se = np.sqrt(np.diag(cov))

print("Normal: omega=%.12g alpha=%.12g gamma=%.12g beta=%.12g"
      % (omega, alpha, gamma, beta))
print("llf=%.12g aic=%.12g bic=%.12g backcast=%.12g" % (llf, aic, bic, bc))
print("persistence alpha+gamma/2+beta=%.12g" % (alpha + gamma / 2 + beta))

# --- probe-1: h₁ 初始化约定 (A: γ/2·bc vs B: γ·I(ε₀<0)·bc) ---
# 残差 = y_dm (mean='Zero' 模型残差, 非模拟真值 eps — probe_gjr_h1.py 教训)
r = y_dm
hA = np.empty(T); hB = np.empty(T)
hA_prev = hB_prev = bc
e2A_prev = e2B_prev = bc
indA_prev = 0.5           # A: t=0 用 0.5 (arch backcast 惯例)
indB_prev = 1.0 if r[0] < 0 else 0.0  # B: 用首残差符号
for t in range(T):
    hA[t] = omega + alpha * e2A_prev + gamma * indA_prev * e2A_prev + beta * hA_prev
    hB[t] = omega + alpha * e2B_prev + gamma * indB_prev * e2B_prev + beta * hB_prev
    e2A_prev = r[t] ** 2; hA_prev = hA[t]
    e2B_prev = r[t] ** 2; hB_prev = hB[t]
    indA_prev = indB_prev = 1.0 if r[t] < 0 else 0.0
dA = np.max(np.abs(hA - h_arch))
dB = np.max(np.abs(hB - h_arch))
print("probe-1 h1: A(gamma/2*bc) maxdiff=%.3e | B(I(eps0)) maxdiff=%.3e" % (dA, dB))
assert dA < 1e-10, "probe-1: convention A mismatch (dA=%.3e)" % dA

# --- forecast: analytic h_{T+1..T+10} + probe-2 φ 约定 ---
fc = res.forecast(horizon=10, method='analytic', reindex=False)
fc_var = np.asarray(fc.variance.values[0], dtype=float)

h_T = float(h_arch[-1]); e_T = float(y_dm[-1])
phi_half = alpha + gamma / 2.0 + beta   # probe-2 理论: E[I*eps2]=h/2
phi_zero = alpha + beta                 # 零残差递归伪影
f1 = omega + alpha * e_T ** 2 + gamma * (e_T < 0) * e_T ** 2 + beta * h_T
fc_half = [f1]; fc_zero = [f1]
for k in range(1, 10):
    fc_half.append(omega + phi_half * fc_half[-1])
    fc_zero.append(omega + phi_zero * fc_zero[-1])
d_half = np.max(np.abs(np.array(fc_half) - fc_var))
d_zero = np.max(np.abs(np.array(fc_zero) - fc_var))
print("probe-2 forecast: phi=alpha+gamma/2+beta maxdiff=%.3e | phi=alpha+beta maxdiff=%.3e"
      % (d_half, d_zero))
assert d_half < 1e-10, "probe-2: phi_half mismatch (d=%.3e)" % d_half

# --- 3. arch t 分布基准 (ν 联合估计, t(6) 冲击数据 seed 63) ---
rng2 = np.random.default_rng(63)
nu_true = 6.0
h2 = np.empty(T)
eps2 = np.empty(T)
h_prev = uc_var
for t in range(T):
    if t == 0:
        h2[t] = h_prev
    else:
        h2[t] = w + al * eps2[t - 1] ** 2 \
               + ga * (eps2[t - 1] < 0) * eps2[t - 1] ** 2 + be * h_prev
    z = rng2.standard_t(nu_true) / np.sqrt(nu_true / (nu_true - 2.0))  # 标准化 t
    eps2[t] = np.sqrt(h2[t]) * z
    h_prev = h2[t]
y2 = mu + eps2
y2_dm = y2 - np.mean(y2)

am_t = arch_model(y2_dm, mean='Zero', vol='GARCH', p=1, o=1, q=1, dist='t')
res_t = am_t.fit(disp='off', show_warning=False)
pt = res_t.params
print("t: omega=%.12g alpha=%.12g gamma=%.12g beta=%.12g nu=%.12g llf=%.12g"
      % (pt['omega'], pt['alpha[1]'], pt['gamma[1]'], pt['beta[1]'], pt['nu'],
         res_t.loglikelihood))


def fmt(v):
    return repr(float(v))


lines = []
lines.append("// 自动生成: tests/fixtures/timeseries/verify_gjr.py (arch 8.0.0)")
lines.append("// 勿手改 — 重新生成请运行该脚本")
lines.append("// 参数顺序 (spec §2.0.4): {omega, alpha, gamma, beta, [nu]} = arch 顺序")
lines.append("// probe-1: h1 = omega+(alpha+gamma/2+beta)*bc (arch t=0 用 0.5*bc)")
lines.append("// probe-2: phi = alpha+gamma/2+beta (E[I(z<0)*eps2]=h/2, 已验证)")
lines.append("#pragma once")
lines.append("namespace cpphub { inline namespace v1 { namespace timeseries {")
lines.append("namespace garch { namespace gjr_baseline {")
lines.append("")
lines.append("constexpr Size T = %d;" % T)
lines.append("")
lines.append("// 模拟数据 y_t = mu + eps_t (GJR: w=0.05,al=0.05,ga=0.15,be=0.80,mu=0.05, seed 62)")
lines.append("constexpr Real DATA[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(y[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// arch mean='Zero' vol='GARCH' p=1,o=1,q=1 dist='normal' 基准")
lines.append("constexpr Real BACKCAST = %s;" % fmt(bc))
lines.append("constexpr Real OMEGA = %s;" % fmt(omega))
lines.append("constexpr Real ALPHA = %s;" % fmt(alpha))
lines.append("constexpr Real GAMMA = %s;" % fmt(gamma))
lines.append("constexpr Real BETA = %s;" % fmt(beta))
lines.append("constexpr Real PERSIST = %s;  // alpha+gamma/2+beta (G10)" % fmt(alpha + gamma / 2 + beta))
lines.append("constexpr Real LLF = %s;" % fmt(llf))
lines.append("constexpr Real AIC = %s;" % fmt(aic))
lines.append("constexpr Real BIC = %s;" % fmt(bic))
lines.append("constexpr Real VCOV[4][4] = {")
for i in range(4):
    lines.append("    {%s}," % ", ".join(fmt(cov[i][j]) for j in range(4)))
lines.append("};")
lines.append("constexpr Real SE[4] = {%s};"
             % ", ".join(fmt(se[j]) for j in range(4)))
lines.append("")
lines.append("// arch conditional variance h_t (T=1000, 递归对照用)")
lines.append("constexpr Real H[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(h_arch[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// forecast 输入末态 + arch analytic h_{T+1..T+10}")
lines.append("constexpr Real H_T = %s;" % fmt(h_T))
lines.append("constexpr Real E_T = %s;" % fmt(e_T))
lines.append("constexpr Real UC_VAR = %s;  // omega/(1-alpha-gamma/2-beta)" % fmt(omega / (1.0 - alpha - gamma / 2.0 - beta)))
lines.append("constexpr Real FC10[10] = {%s};"
             % ", ".join(fmt(v) for v in fc_var))
lines.append("")
lines.append("// t(6) 冲击模拟数据 (seed 63, GJR 同参数, mu=0.05)")
lines.append("constexpr Real TDATA[%d] = {" % T)
for i in range(0, T, 4):
    row = ", ".join(fmt(y2[j]) for j in range(i, min(i + 4, T)))
    lines.append("    " + row + ",")
lines.append("};")
lines.append("")
lines.append("// arch mean='Zero' vol='GARCH' p=1,o=1,q=1 dist='t' 基准 (ν 联合估计)")
lines.append("constexpr Real T_OMEGA = %s;" % fmt(pt['omega']))
lines.append("constexpr Real T_ALPHA = %s;" % fmt(pt['alpha[1]']))
lines.append("constexpr Real T_GAMMA = %s;" % fmt(pt['gamma[1]']))
lines.append("constexpr Real T_BETA = %s;" % fmt(pt['beta[1]']))
lines.append("constexpr Real T_NU = %s;" % fmt(pt['nu']))
lines.append("constexpr Real T_LLF = %s;" % fmt(res_t.loglikelihood))
lines.append("")
lines.append("}}}}}  // gjr_baseline::garch::timeseries::v1::cpphub")

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
print("written:", OUT)
