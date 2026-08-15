# probe_gjr_h1.py — 排查 arch GJR h 序列差异 (2026-08-15)
import numpy as np
from arch import arch_model

rng = np.random.default_rng(62)
T = 1000
w, al, ga, be, mu = 0.05, 0.05, 0.15, 0.80, 0.05
uc_var = w / (1.0 - al - ga / 2.0 - be)
h = np.empty(T); eps = np.empty(T)
h_prev = uc_var
for t in range(T):
    if t == 0:
        h[t] = h_prev
    else:
        h[t] = w + al * eps[t-1]**2 + ga * (eps[t-1] < 0) * eps[t-1]**2 + be * h_prev
    eps[t] = np.sqrt(h[t]) * rng.standard_normal()
    h_prev = h[t]
y = mu + eps
y_dm = y - np.mean(y)

am = arch_model(y_dm, mean='Zero', vol='GARCH', p=1, o=1, q=1, dist='normal')
res = am.fit(disp='off', show_warning=False)
p = res.params
omega, alpha, gamma, beta = (float(p['omega']), float(p['alpha[1]']),
                             float(p['gamma[1]']), float(p['beta[1]']))
bc = float(am.volatility.backcast(y_dm))
h_arch = res.conditional_volatility ** 2

print("params:", omega, alpha, gamma, beta, "bc:", bc)
print("eps[0] =", eps[0], "-> I =", 1.0 if eps[0] < 0 else 0.0)

# 关键: 模型残差 = y_dm 本身 (mean='Zero'), 非模拟真值 eps
r = y_dm
cands = {
    'A: om+(al+ga/2+be)*bc': omega + (alpha + gamma/2 + beta) * bc,
    'B: om+(al+be)*bc (ga 跳过)': omega + (alpha + beta) * bc,
    'C: om+al*bc+ga*I0*bc+be*bc': omega + alpha*bc + gamma*(r[0]<0)*bc + beta*bc,
    'D: bc 直接': bc,
}
print("h_arch[0] =", h_arch[0])
for k, v in cands.items():
    print(f"  {k}: {v:.12g}  diff={abs(v-h_arch[0]):.3e}")

# 用候选 A/B/C 分别递归全程 (残差用 r=y_dm), 打印首个显著偏离点
for name, h0 in [('A', cands['A: om+(al+ga/2+be)*bc']),
                 ('B', cands['B: om+(al+be)*bc (ga 跳过)']),
                 ('C', cands['C: om+al*bc+ga*I0*bc+be*bc'])]:
    hm = np.empty(T); hm[0] = h0
    for t in range(1, T):
        hm[t] = omega + alpha*r[t-1]**2 + gamma*(r[t-1]<0)*r[t-1]**2 + beta*hm[t-1]
    d = np.abs(hm - h_arch)
    i = int(np.argmax(d))
    print(f"conv {name}: maxdiff={d.max():.3e} at t={i}; first>1e-10 at "
          f"{np.argmax(d>1e-10)}; h_m[{i}]={hm[i]:.6g} h_arch[{i}]={h_arch[i]:.6g}")

# 对照 arch 源码: 打印 compute_variance 片段
import inspect
src = inspect.getsource(type(am.volatility).compute_variance)
print("--- arch compute_variance source ---")
print(src)
