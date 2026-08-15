# probe_egarch_sim.py — 排查 arch EGARCH simulation 路径 vs 确定性递归差异
import numpy as np
from arch import arch_model

# 复现 verify_egarch.py 的数据 (seed 52)
rng = np.random.default_rng(52)
T = 1000
E_ABS_Z = np.sqrt(2.0 / np.pi)
w, al, be, ga, mu = 0.0, -0.05, 0.94, 0.15, 0.05
ln_h = np.empty(T); eps = np.empty(T)
ln_h_prev = 0.0; z_prev = 0.0
for t in range(T):
    if t > 0:
        ln_h_prev = w + be * ln_h_prev + al * z_prev + ga * (abs(z_prev) - E_ABS_Z)
    ln_h[t] = ln_h_prev
    sig = np.sqrt(np.exp(ln_h[t]))
    eps[t] = sig * rng.standard_normal()
    z_prev = eps[t] / sig
y = mu + eps
y_dm = y - np.mean(y)

am = arch_model(y_dm, mean='Zero', vol='EGARCH', p=1, o=1, q=1, dist='normal')
res = am.fit(disp='off', show_warning=False)
p = res.params
om, a_sym, g_asym, be = float(p['omega']), float(p['alpha[1]']), float(p['gamma[1]']), float(p['beta[1]'])
print("params: w=%.6f a_sym=%.6f g_asym=%.6f b=%.6f" % (om, a_sym, g_asym, be))

h_arch = np.asarray(res.conditional_volatility) ** 2
ln_h_T = np.log(h_arch[-1])
z_T = y_dm[-1] / np.sqrt(h_arch[-1])

# --- 确定性递归 (mine) ---
ln_det = [om + be * ln_h_T + g_asym * z_T + a_sym * (abs(z_T) - E_ABS_Z)]
for k in range(1, 10):
    ln_det.append(om + be * ln_det[-1])
ln_det = np.array(ln_det)

# --- 手动 MC (arch 同款递归, N(0,1) draws) ---
S = 500000
rng2 = np.random.default_rng(7)
shocks = rng2.standard_normal((S, 10))
ln_paths = np.empty((S, 10))
ln_prev = np.full(S, ln_det[0])  # 第一期确定
for k in range(1, 10):
    z = shocks[:, k - 1]
    ln_prev = om + be * ln_prev + g_asym * z + a_sym * (np.abs(z) - E_ABS_Z)
    ln_paths[:, k] = ln_prev
ln_paths[:, 0] = ln_det[0]

print("\ndeterministic ln h:", np.round(ln_det, 4))
print("manual MC mean ln :", np.round(ln_paths.mean(0), 4))
print("manual MC med  ln :", np.round(np.median(ln_paths, 0), 4))

# --- arch simulation ---
fc = res.forecast(horizon=10, method='simulation', reindex=False)
sv = np.asarray(fc.simulations.variances, dtype=float)
print("\nsim.variances shape:", sv.shape)
ln_arch = np.log(sv[:, 0, :])
print("arch mean ln       :", np.round(ln_arch.mean(0), 4))
print("arch median ln     :", np.round(np.median(ln_arch, 0), 4))
print("arch mean h        :", np.round(sv[:, 0, :].mean(0), 4))
print("reported fc var    :", np.round(np.asarray(fc.variance.values[0], dtype=float), 4))
