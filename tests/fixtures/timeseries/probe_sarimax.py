# -*- coding: utf-8 -*-
# probe_sarimax.py - statsmodels SARIMAX 语义探针 (Phase 8 M2 前置)
#
# SPEC §3.2 Step 1 四项任务:
#   (a) 嵌入差分 T/R/Z 构造 (sarimax.py initial_transition/initial_design/
#       initial_selection + update L1682-1694 逐元素)
#   (b) loglikelihood_burn 默认公式 (initialize_default L655-659 实录:
#       burn = k_states − _k_order = _k_states_diff = d, enforce_
#       stationarity=True 前提; = k_states 当 False)
#   (c) transform_params 细节 (L1286-1388: AR=constrain_stationary_
#       univariate (Monahan); MA=**负号**同函数; σ²=平方)
#   (d) llf 口径 (burn 应用 + 常数 + 手工 DK 更新形 KF 重算对齐)
#
# 已源码实锤 (本轮 grep, 冻结待数值自证):
#   - 布局: [差分块 (d) | ARMA companion 块 (_k_order) | exog]
#     _k_order = max(p, q+1) (k_ar=0 且 k_ma=0 时归 0, L425-430)
#   - companion_matrix: AR 系数在第一列 (上 Hessenberg) = Harvey 形
#     (hamilton_representation=False 默认; True 才转置)
#   - Z = [1]*d + [1] + [0]*(_k_order−1); R = [0]*d + [1, θ...]';
#     T: 差分块上三角全 1 (triu_indices), T[:d, d] = 1 (链接列),
#        ARMA 块 companion + AR 填 T[d:d+p, d] (transition_ar_params_idx)
#   - 初始化 (Initialization): approximate_diffuse_variance 默认 1e6
#     (initialization.py L158/182 — SPEC κ=1e6 口径正确; state_regression
#     才用 1e10); 差分块 diffuse, ARMA 块 stationary (enforce=True),
#     均值 = solve(I−T, c) (L720)
#
# 数值实证冻结 (2026-08-22 全绿运行, 本脚本可复现):
#   [A] ARIMA(1,1,1) 黄金锚逐元素对齐 (dT=0, dR=0):
#       Z=[1,1,0]  T=[[1,1,0],[0,φ,1],[0,0,0]]  R=[0,1,θ]'
#       drift 经 c_t[_k_states_diff] (ARMA 块首行) 注入
#   [B] burn = d 八组全验 (101/200/001/111/212/011/021/321);
#       enforce_stationarity=False → burn = k_states (全 diffuse)
#   [C] MA = −constrain(u) 负号实锤 (小值区 θ≈+u); 往返 5.55e-17
#   [D] DK 更新形 KF (C++ ss_model.hpp 参照) 数值等价实证:
#       ARMA(1,1) llf diff = 2.81e-9 (T=120); ARIMA(1,1,1) = 2.81e-8
#       (diffuse κ=1e6 放大尾差; 相对 ~1e-10); v_t/F_t diff ~3e-10
#       burn 应用 = 贡献求和起点 (滤波全程运行, 非丢观测; diff 3.7e-14)
#   [E] 初始化结构: Pinf 全零 (approximate diffuse 无 exact 标记);
#       P1 = Pstar = diag(κ=1e6, 差分块) ⊕ Lyapunov(ARMA 块)
#       (手工 Lyapunov 2.22e-16 对齐); a1 = 0
#   [D2] 集中化陷阱: concentrate_scale=True 时 k_params=p+q (params 无 σ²,
#       params[-1] 是 θ!); False 时 [φ,θ,σ²] — 对拍必须同口径 (K4)
#
# 用法: & C:\Users\Peng\.conda\envs\open-webui\python.exe probe_sarimax.py
import numpy as np

np.set_printoptions(precision=8, suppress=False, linewidth=120)


def sq(x):
    """statespace 矩阵去 time-varying 尾维 (design/transition/... 均 (…, nobs=1))"""
    return np.asarray(x).squeeze()


def gen_arma11(T=120, seed=88):
    """ARMA(1,1): y = 0.6y_{t-1} + e + 0.4e_{t-1}"""
    rng = np.random.RandomState(seed)
    burn = 200
    n = T + burn
    e = rng.standard_normal(n)
    y = np.zeros(n)
    for t in range(1, n):
        y[t] = 0.6 * y[t - 1] + e[t] + 0.4 * e[t - 1]
    return y[burn:]


def probe_matrices():
    """[A] 系统矩阵逐元素: ARMA(1,1) + ARIMA(1,1,1) 嵌入差分"""
    from statsmodels.tsa.statespace.sarimax import SARIMAX

    y = gen_arma11()
    print("=" * 74)
    print("[A1] ARMA(1,1) (d=0): m = max(p, q+1) = 2")
    mod = SARIMAX(y, order=(1, 0, 1),
                  enforce_stationarity=True, enforce_invertibility=True,
                  concentrate_scale=False)
    m = mod.ssm
    Z = np.asarray(m.design)[0]
    Tt = np.asarray(m.transition)
    R = np.asarray(m.selection)
    Q = np.asarray(m.state_cov)
    H = np.asarray(m.obs_cov)
    print(f"  k_states={mod.k_states} k_posdef={mod.k_posdef} "
          f"_k_order={mod._k_order} _k_states_diff={mod._k_states_diff}")
    print(f"  Z = {Z}")
    print(f"  T =\n{Tt}")
    print(f"  R = {R.flatten()}")
    print(f"  Q = {Q.flatten()}   H = {H.flatten()}")
    print(f"  (骨架; AR/MA 系数在 update() 填入: T[0:1,0]=φ, R[1:2,0]=θ)")

    print()
    print("[A2] ARIMA(1,1,1) 嵌入差分 (simple_differencing=False): m = 1+2 = 3")
    mod11 = SARIMAX(y, order=(1, 1, 1),
                    enforce_stationarity=True, enforce_invertibility=True,
                    concentrate_scale=False)
    m11 = mod11.ssm
    Z = sq(m11.design)
    Tt = sq(m11.transition)
    R = sq(m11.selection)
    print(f"  k_states={mod11.k_states} _k_order={mod11._k_order} "
          f"_k_states_diff={mod11._k_states_diff}")
    print(f"  Z = {Z}   (期望 [1,1,0])")
    print(f"  T =\n{Tt}")
    print(f"  R = {R.flatten()}   (期望 [0,1,0] 骨架)")
    print("  黄金锚 (SPEC §3.2, 参数 φ=0.5, θ=0.3, σ²=1):")
    phi, theta = 0.5, 0.3
    T_gold = np.array([[1.0, 1.0, 0.0],
                       [0.0, phi, 1.0],
                       [0.0, 0.0, 0.0]])
    R_gold = np.array([0.0, 1.0, theta])
    print(f"  T_gold =\n{T_gold}")
    print(f"  R_gold = {R_gold}")
    # update 填入后对照
    mod11.update(np.array([phi, theta, 1.0]))
    T_u = sq(mod11.ssm.transition)
    R_u = sq(mod11.ssm.selection)
    print(f"  update([φ,θ,σ²]) 后 diff: "
          f"dT={np.max(np.abs(T_u - T_gold)):.2e} "
          f"dR={np.max(np.abs(R_u.flatten() - R_gold)):.2e}")
    assert np.max(np.abs(T_u - T_gold)) < 1e-14
    assert np.max(np.abs(R_u.flatten() - R_gold)) < 1e-14
    print("  ✓ 黄金锚逐元素对齐")

    print()
    print("[A3] ARIMA(0,1,1) + trend='c' (d=1, drift): state_intercept 注入")
    modd = SARIMAX(y, order=(0, 1, 1), trend="c",
                   enforce_stationarity=True, enforce_invertibility=True,
                   concentrate_scale=False)
    print(f"  k_states={modd.k_states} _k_order={modd._k_order} "
          f"(_k_diff={modd._k_diff})")
    modd.update(np.array([0.2, 0.5, 1.0]))  # [c, θ, σ²]
    Td = sq(modd.ssm.transition)
    Zd = sq(modd.ssm.design)
    Rd = sq(modd.ssm.selection).flatten()
    cd = np.asarray(modd.ssm.state_intercept)
    print(f"  Z = {Zd}   T =\n{Td}   R = {Rd}")
    print(f"  state_intercept shape={cd.shape} 非零首行 = "
          f"{cd[:, 0] if cd.ndim == 2 else cd}")
    print("  → drift c 经 c_t[_k_states_diff] (ARMA 块首行) 注入")


def probe_burn():
    """[B] loglikelihood_burn 默认公式"""
    from statsmodels.tsa.statespace.sarimax import SARIMAX
    y = gen_arma11(T=60)
    print("=" * 74)
    print("[B] loglikelihood_burn 默认公式 (enforce_stationarity=True)")
    print("    源码: burn = k_states − _k_order = _k_states_diff = d")
    for order in [(1, 0, 1), (2, 0, 0), (0, 0, 1), (1, 1, 1),
                  (2, 1, 2), (0, 1, 1), (0, 2, 1), (3, 2, 1)]:
        mod = SARIMAX(y, order=order, enforce_stationarity=True,
                      enforce_invertibility=True, concentrate_scale=False)
        p, d, q = order
        expect = d
        ok = "✓" if mod.loglikelihood_burn == expect else "✗✗✗"
        print(f"  ARIMA{order}: k_states={mod.k_states} "
              f"_k_order={mod._k_order} burn={mod.loglikelihood_burn} "
              f"expect={expect} {ok}")
        assert mod.loglikelihood_burn == expect
    mod_ns = SARIMAX(y, order=(1, 1, 1), enforce_stationarity=False,
                     enforce_invertibility=False, concentrate_scale=False)
    print(f"  enforce_stationarity=False: burn={mod_ns.loglikelihood_burn} "
          f"(= k_states={mod_ns.k_states}, 全 diffuse)")
    assert mod_ns.loglikelihood_burn == mod_ns.k_states


def probe_transform():
    """[C] transform/untransform 细节"""
    from statsmodels.tsa.statespace.sarimax import SARIMAX
    y = gen_arma11(T=60)
    print("=" * 74)
    print("[C] transform_params: AR=Monahan, MA=**负号**Monahan, σ²=平方")
    mod = SARIMAX(y, order=(2, 0, 2), trend="c",
                  enforce_stationarity=True, enforce_invertibility=True,
                  concentrate_scale=False)
    # params = [trend(1), ar(2), ma(2), σ²]
    unc = np.array([0.3, -0.2, 0.5, -0.1, 0.4, 0.8])
    con = mod.transform_params(unc)
    print(f"  unconstrained = {unc}")
    print(f"  constrained   = {con}")
    # 手工 Monahan (constrain_stationary_univariate)
    from statsmodels.tsa.statespace.tools import (
        constrain_stationary_univariate, unconstrain_stationary_univariate)
    ar_man = constrain_stationary_univariate(unc[1:3])
    ma_man = -constrain_stationary_univariate(unc[3:5])
    print(f"  手工: trend 直通 = {con[0] == unc[0]}; "
          f"AR diff = {np.max(np.abs(con[1:3] - ar_man)):.2e}; "
          f"MA diff = {np.max(np.abs(con[3:5] - ma_man)):.2e}; "
          f"σ² = unc[-1]² = {con[-1] == unc[-1]**2}")
    assert con[0] == unc[0]
    assert np.max(np.abs(con[1:3] - ar_man)) < 1e-15
    assert np.max(np.abs(con[3:5] - ma_man)) < 1e-15
    assert con[-1] == unc[-1] ** 2
    # 往返
    back = mod.untransform_params(con)
    print(f"  往返 untransform(transform(x)) diff = "
          f"{np.max(np.abs(back - unc)):.2e}")
    assert np.max(np.abs(back - unc)) < 1e-12
    # MA 负号约定实证: θ = −constrain(u) (源码 L1346 负号; Monahan 小值
    # 线性区 constrain(u) ≈ −u → θ ≈ +u)
    mod1 = SARIMAX(y, order=(0, 0, 1), enforce_stationarity=True,
                   enforce_invertibility=True, concentrate_scale=False)
    c1 = mod1.transform_params(np.array([0.1, 0.5]))
    expect1 = -constrain_stationary_univariate(np.array([0.1]))
    print(f"  MA(1) 单参数: u=0.1 → θ={c1[0]:.6f} "
          f"(= −constrain(u) = {expect1[0]:.6f}; 小值区 ≈ +u)")
    assert abs(c1[0] - expect1[0]) < 1e-15
    assert c1[0] > 0  # 小值区 θ ≈ +u (负号×负号)


def manual_kf(y, Z, T, R, Q, H, a1, P1, c=None, burn=0):
    """手工 Kalman 滤波 (DK Ch4 更新形) — C++ ss_model.hpp 的参照实现

    P_{t+1} = T P_t (T − K_t Z)' + R Q R'  (更新形; 与标准形代数等价)
    loglik 贡献 (t ≥ burn): −½(log F_t + v_t²/F_t + log 2π)
    """
    T_n = len(y)
    m = T.shape[0]
    if c is None:
        c = np.zeros(m)
    a = a1.astype(float).copy()
    P = P1.astype(float).copy()
    Zf = np.asarray(Z).flatten()
    Rm = np.asarray(R).reshape(m, -1)
    Qm = np.asarray(Q).reshape(Rm.shape[1], Rm.shape[1])
    h = float(np.asarray(H).flatten()[0])
    contribs = []
    ll = 0.0
    innov = np.zeros(T_n)
    Fv = np.zeros(T_n)
    for t in range(T_n):
        v = y[t] - Zf @ a
        F = float(Zf @ P @ Zf + h)
        innov[t] = v
        Fv[t] = F
        if t >= burn:
            ll_t = -0.5 * (np.log(F) + v * v / F + np.log(2 * np.pi))
            contribs.append(ll_t)
            ll += ll_t
        K = T @ P @ Zf / F           # (m,)
        a = T @ a + K * v + c
        P = T @ P @ (T - np.outer(K, Zf)).T + Rm @ Qm @ Rm.T
    return ll, np.array(contribs), innov, Fv


def probe_llf():
    """[D]+[E]+[F] llf 口径 + 初始化 + 手工 KF 对齐"""
    from statsmodels.tsa.statespace.sarimax import SARIMAX
    print("=" * 74)
    print("[D] llf 口径: burn 应用 + 常数 + 手工 DK 更新形 KF 重算")

    y = gen_arma11(T=120)

    # --- ARMA(1,1), burn=0, stationary 初始化 ---
    mod = SARIMAX(y, order=(1, 0, 1), enforce_stationarity=True,
                  enforce_invertibility=True, concentrate_scale=False)
    params = np.array([0.55, 0.35, 0.9])  # [φ, θ, σ²]
    mod.update(params)
    m = mod.ssm
    # 初始化 [E]: Initialization.__call__ → (a1, P_inf, P_star)
    init = m.initialization
    a1, Pinf, Pstar = init(model=m)  # stationary 块需 model 求解
    print("  [E] ARMA(1,1) stationary 初始化:")
    print(f"      a1 = {a1}")
    print(f"      Pinf =\n{Pinf}")
    print(f"      Pstar =\n{Pstar}")
    print(f"      (Pinf 全零 → 纯 stationary; Pstar = Lyapunov 解)")
    # 手工 Lyapunov 验证: T=[[φ,1],[0,0]], R=[1,θ]', Q=[σ²]
    phi, th, s2 = params
    Tm = np.array([[phi, 1.0], [0.0, 0.0]])
    Rm = np.array([[1.0], [th]])
    RQR = Rm @ np.array([[s2]]) @ Rm.T
    # 迭代 Lyapunov
    P = RQR.copy()
    for _ in range(200):
        P = Tm @ P @ Tm.T + RQR
    print(f"      手工 Lyapunov diff = {np.max(np.abs(P - Pstar)):.2e}")
    assert np.max(np.abs(P - Pstar)) < 1e-10

    ll_sm = mod.loglike(params)
    Z = sq(m.design)
    Tt = sq(m.transition)
    R = sq(m.selection)
    Q = sq(m.state_cov)
    H = np.asarray(m.obs_cov).reshape(1, 1)
    ll_man, contribs, innov, Fv = manual_kf(y, Z, Tt, R, Q, H, a1, Pstar,
                                            burn=mod.loglikelihood_burn)
    print(f"  ARMA(1,1): loglike(sm) = {ll_sm:.10f}  手工 = {ll_man:.10f}  "
          f"diff = {abs(ll_sm - ll_man):.2e}")
    # SPEC §3.2 测试要点: KF loglik vs SARIMAX llf 1e-8 (浮点尾差级:
    # statsmodels 内部 BLAS 顺序/P 更新形差异)
    assert abs(ll_sm - ll_man) < 1e-8
    # 逐贡献对照 (results.llf 与 loglike 同; contribs 经 filter 输出)
    res = mod.filter(params)
    v_sm = np.asarray(res.forecasts_error)[0]     # v_t
    F_sm = np.asarray(res.forecasts_error_cov)[0, 0]  # F_t
    print(f"      v_t diff = {np.max(np.abs(v_sm - innov)):.2e}  "
          f"F_t diff = {np.max(np.abs(F_sm - Fv)):.2e}")
    assert np.max(np.abs(v_sm - innov)) < 1e-6  # 误差随 T 累积, llf 对齐为主锚
    assert np.max(np.abs(F_sm - Fv)) < 1e-6

    # --- ARIMA(1,1,1), burn=1, 混合初始化 (diffuse d 块 + stationary ARMA) ---
    print()
    mod11 = SARIMAX(y, order=(1, 1, 1), enforce_stationarity=True,
                    enforce_invertibility=True, concentrate_scale=False)
    params11 = np.array([0.55, 0.35, 0.9])
    mod11.update(params11)
    m11 = mod11.ssm
    init11 = m11.initialization
    a1_11, Pinf_11, Pstar_11 = init11(model=m11)
    print("  [E] ARIMA(1,1,1) 混合初始化 (差分块 diffuse κ=1e6 + ARMA 块 stationary):")
    print(f"      a1 = {a1_11}")
    print(f"      Pinf =\n{Pinf_11}")
    print(f"      Pstar =\n{Pstar_11}")
    # Pinf 对角 = 差分块 1 (exact diffuse 标记)? approximate diffuse 时
    # Pinf 全零而 Pstar 对角含 κ → 数值实证
    print(f"      burn = {mod11.loglikelihood_burn} (= d = 1)")
    ll_sm11 = mod11.loglike(params11)
    # approximate diffuse: P1 = Pstar (含 κ=1e6 对角于差分块)
    Z11 = sq(m11.design)
    T11 = sq(m11.transition)
    R11 = sq(m11.selection)
    Q11 = np.asarray(m11.state_cov).reshape(1, 1)
    H11 = np.asarray(m11.obs_cov).reshape(1, 1)
    P1_11 = Pinf_11 + Pstar_11  # approximate diffuse: 全部进 Pstar
    ll_man11, _, innov11, Fv11 = manual_kf(y, Z11, T11, R11, Q11, H11,
                                           a1_11, P1_11,
                                           burn=mod11.loglikelihood_burn)
    print(f"  ARIMA(1,1,1): loglike(sm) = {ll_sm11:.10f}  手工 = "
          f"{ll_man11:.10f}  diff = {abs(ll_sm11 - ll_man11):.2e}")
    # diffuse κ=1e6 放大浮点尾差 (~1e-8); 相对 diff ~1e-10 → 等价性成立
    if abs(ll_sm11 - ll_man11) < 1e-7:
        print("      ✓ P1 = Pinf + Pstar (approximate diffuse 合成, 差分块 "
              "κ=1e6 对角 + ARMA 块 Lyapunov) 数值等价成立")
    else:
        print("      ✗ 需检查: P1 合成方式或 diffuse 块初值")
    # v_t/F_t 对照
    res11 = mod11.filter(params11)
    v11_sm = np.asarray(res11.forecasts_error)[0]
    F11_sm = np.asarray(res11.forecasts_error_cov)[0, 0]
    print(f"      v_t diff = {np.max(np.abs(v11_sm - innov11)):.2e}  "
          f"F_t diff = {np.max(np.abs(F11_sm - Fv11)):.2e}")

    # --- [B 补充] burn 应用方式: llf(burn=1) vs llf(burn=0) 差 = 首期贡献 ---
    mod_b0 = SARIMAX(y, order=(1, 1, 1), enforce_stationarity=True,
                     enforce_invertibility=True, concentrate_scale=False,
                     loglikelihood_burn=0)
    ll_b0 = mod_b0.loglike(params11)
    ll_b1 = mod11.loglike(params11)
    # 手工全贡献和
    ll_full, contribs_full, _, _ = manual_kf(y, Z11, T11, R11, Q11, H11,
                                             a1_11, P1_11, burn=0)
    print(f"  burn 应用: llf(burn=0) − llf(burn=1) = {ll_b0 - ll_b1:.10f}")
    print(f"             手工首期贡献 = {contribs_full[0]:.10f}  "
          f"diff = {abs((ll_b0 - ll_b1) - contribs_full[0]):.2e}")
    assert abs(ll_b0 - ll_full) < 1e-7  # 同 KF 尾差量级 (diffuse κ 放大)
    print("      ✓ burn = 贡献求和起点 (非丢弃观测, 滤波全程运行)")


def probe_concentrate():
    """[D 补充] concentrate_scale 口径 (SPEC K4: σ² 非集中化 vs statsmodels 默认)"""
    from statsmodels.tsa.statespace.sarimax import SARIMAX
    y = gen_arma11(T=120)
    print("=" * 74)
    print("[D2] concentrate_scale 口径 (K4)")
    mod_c = SARIMAX(y, order=(1, 0, 1), concentrate_scale=True)
    mod_f = SARIMAX(y, order=(1, 0, 1), concentrate_scale=False)
    res_c = mod_c.fit(disp=False)
    res_f = mod_f.fit(disp=False)
    # ⚠️ 集中化 k_params = p+q (无 σ², params = [φ, θ]); 非集中化 params =
    # [φ, θ, σ²] — res_c.params[-1] 是 θ 不是 σ² (打印陷阱)
    print(f"  集中化 (默认): k_params={len(res_c.params)} llf = "
          f"{res_c.llf:.6f}  params = {res_c.params}")
    print(f"  非集中化:      k_params={len(res_f.params)} llf = "
          f"{res_f.llf:.6f}  params = {res_f.params}")
    print(f"  llf diff = {abs(res_c.llf - res_f.llf):.2e} (优化收敛差)")
    print(f"  params(φ,θ) diff = "
          f"{np.max(np.abs(res_c.params[:2] - res_f.params[:2])):.2e}")
    print("  → statsmodels fit 默认 concentrate_scale=True (σ² 集中出,"
          " params 末端无 σ²); 本库 K4 非集中化 (σ² 为优化参数) —"
          " 对拍需同口径 (llf 主锚, 交叉验证以 K4 声明)")


def main():
    probe_matrices()
    probe_burn()
    probe_transform()
    probe_llf()
    probe_concentrate()
    print("=" * 74)
    print("全部探明项冻结: 见文件头注释 + SPEC §3.2 Step 1 裁决")


if __name__ == "__main__":
    main()
