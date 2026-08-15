# 探查 arch GARCH(1,1) 的 h1 约定: h1 = backcast 还是 h1 = w + (a+b)*backcast
# 同时验证: mean='Zero' 下 resids 是否等于输入原样
import numpy as np
from arch import arch_model

rng = np.random.default_rng(42)
T = 500
w, a, b = 2e-6, 0.10, 0.85
h = np.empty(T)
eps = np.empty(T)
h_prev = w / (1.0 - a - b)  # 从无条件方差起步
for t in range(T):
    h[t] = w + a * eps[t - 1] ** 2 + b * h_prev if t > 0 else h_prev
    eps[t] = np.sqrt(h[t]) * rng.standard_normal()
    h_prev = h[t]

# arch: demean 不涉及 (mean='Zero')
am = arch_model(eps, mean='Zero', vol='GARCH', p=1, q=1, dist='normal')
res = am.fit(disp='off')
p = res.params
print("params:", dict(p))

# backcast 值 (arch 内部)
bc = am.volatility.backcast(eps)
print("backcast:", bc)

# 方式 A: h1 = backcast
om, al, be = p['omega'], p['alpha[1]'], p['beta[1]']
hA = np.empty(T)
hA[0] = bc
for t in range(1, T):
    hA[t] = om + al * eps[t - 1] ** 2 + be * hA[t - 1]

# 方式 B: h1 = w + (a+b)*backcast
hB = np.empty(T)
hB[0] = om + (al + be) * bc
for t in range(1, T):
    hB[t] = om + al * eps[t - 1] ** 2 + be * hB[t - 1]

h_arch = res.conditional_volatility ** 2
print("max|hA - h_arch| (h1=backcast):", np.max(np.abs(hA - h_arch)))
print("max|hB - h_arch| (h1=w+(a+b)bc):", np.max(np.abs(hB - h_arch)))

# llf 复算 (确认常数项): -0.5*sum(log(2pi)+log(h)+e^2/h)
llf_manual = -0.5 * np.sum(np.log(2 * np.pi) + np.log(h_arch) + eps ** 2 / h_arch)
print("llf arch:", res.loglikelihood, "llf manual:", llf_manual)

# robust cov 形状
print("cov type:", res.cov_type)
print("robust cov:\n", res.param_cov)
