# 探查 arch EGARCH t=0 约定 (log 空间 backcast 进入方式)
# 同时验证参数命名: alpha[1]=对称 |z|-E|z|, gamma[1]=非对称 z
import numpy as np
from arch import arch_model

rng = np.random.default_rng(7)
T = 600
om, be, al_sym, ga_asym = 0.10, 0.90, 0.15, -0.10  # arch: alpha=0.15, gamma=-0.10
c = np.sqrt(2 / np.pi)
lh = np.empty(T)
eps = np.empty(T)
lh_prev = om / (1 - be)  # 无条件 log 方差 = 1.0 → var = e
for t in range(T):
    if t == 0:
        lh[t] = lh_prev
    else:
        z_prev = eps[t - 1] / np.sqrt(np.exp(lh_prev))
        lh[t] = om + be * lh_prev + al_sym * (abs(z_prev) - c) + ga_asym * z_prev
    eps[t] = np.exp(lh[t] / 2) * rng.standard_normal()
    lh_prev = lh[t]
mu = 0.05
y = mu + eps
y_dm = y - np.mean(y)

am = arch_model(y_dm, mean='Zero', vol='EGARCH', p=1, o=1, q=1)
res = am.fit(disp='off', show_warning=False)
p = res.params
print("params:", dict(p))
print("llf:", res.loglikelihood)

bc = float(am.volatility.backcast(y_dm))
print("backcast:", bc)

om_a, al_a, ga_a, be_a = p['omega'], p['alpha[1]'], p['gamma[1]'], p['beta[1]']

# 方式 A: ln h1 = ln bc
def rec_a():
    lh = np.empty(T)
    lh[0] = np.log(bc)
    for t in range(1, T):
        z = eps_dm[t - 1] / np.sqrt(np.exp(lh[t - 1]))
        lh[t] = om_a + be_a * lh[t - 1] + al_a * (abs(z) - c) + ga_a * z
    return np.exp(lh)

# 方式 B: ln h1 = om + be*ln bc
def rec_b():
    lh = np.empty(T)
    lh[0] = om_a + be_a * np.log(bc)
    for t in range(1, T):
        z = eps_dm[t - 1] / np.sqrt(np.exp(lh[t - 1]))
        lh[t] = om_a + be_a * lh[t - 1] + al_a * (abs(z) - c) + ga_a * z
    return np.exp(lh)

# 方式 C: ln h1 = om + al*(sqrt(bc)... z0=√bc? no: |z0|-c with z0 from bc? 试 D
# 方式 D: ln h1 = om + al*(|√bc|-c)?? 不合理; 试 z0 = 0: al*(0-c) + ga*0 + be*ln bc + om
def rec_d():
    lh = np.empty(T)
    lh[0] = om_a + be_a * np.log(bc) + al_a * (0.0 - c)
    for t in range(1, T):
        z = eps_dm[t - 1] / np.sqrt(np.exp(lh[t - 1]))
        lh[t] = om_a + be_a * lh[t - 1] + al_a * (abs(z) - c) + ga_a * z
    return np.exp(lh)

h_arch = res.conditional_volatility ** 2
eps_dm = y_dm

# 方式 E: ln h1 = om + be*bc (arch 源码: beta * backcast, 非 log(bc)!)
def rec_e():
    lh = np.empty(T)
    lh[0] = om_a + be_a * bc
    for t in range(1, T):
        z = eps_dm[t - 1] / np.sqrt(np.exp(lh[t - 1]))
        lh[t] = om_a + be_a * lh[t - 1] + al_a * (abs(z) - c) + ga_a * z
    return np.exp(lh)

for name, f in [("A lnbc", rec_a), ("B om+be*lnbc", rec_b), ("E om+be*bc", rec_e)]:
    h = f()
    print(name, "maxdiff:", np.max(np.abs(h - h_arch)))
