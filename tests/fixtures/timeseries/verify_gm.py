# -*- coding: utf-8 -*-
# verify_gm.py - arch 8.0.0 ARCHInMean 对照基准生成 (Phase 7C M0, spec §2.3)
#
# 用途: 固定合成 GARCH-M 数据 + arch ARCHInMean(form='vol'/'var'/'log',
#       rescale=False) 三 form 基准 (params/llf/robust SE), 供
#       test_garch_m_model.cpp 硬编码断言 (容差 1e-8~1e-10, spec §1.3)
# 环境: conda open-webui (arch 8.0.0): C:/Users/Peng/.conda/envs/open-webui/python.exe
#
# 关键对齐 (GM 幻觉点, zivot 同款头注惯例):
#   - rescale=False (spec §2.3 Step5: 关闭重标定, llf 差 n·log(s²) 陷阱)
#   - dist=Normal, 无 lags/x (纯常数均值 + in-mean)
#   - robust SE = Bollerslev-Wooldridge QMLE (GM4 对照 std_errors)
#
# 运行: python verify_gm.py  (打印基准 + 落盘 gm_arch_baselines.json
#       + 生成 tests/unit/timeseries/gm_baseline.inc)
import json
import os

import numpy as np
from arch.univariate import ARCHInMean, GARCH, Normal

HERE = os.path.dirname(os.path.abspath(__file__))
INC_PATH = os.path.normpath(os.path.join(HERE, "..", "..", "unit", "timeseries",
                                         "gm_baseline.inc"))


def make_gm_data() -> np.ndarray:
    """固定合成数据 (RandomState(42), T=1000, vol 形真值 λ=0.50)

    DGP: r_t = 0.05 + 0.50·√h_t + ε_t;  h_t = 0.05 + 0.10·ε²_{t-1} + 0.85·h_{t-1}
    (三 form 用同一数据 — 估计的是不同设定, 非恢复各自 DGP)
    """
    rng = np.random.RandomState(42)
    T = 1000
    mu, lam = 0.05, 0.50
    om, al, be = 0.05, 0.10, 0.85
    h = np.empty(T)
    r = np.empty(T)
    h_prev = om / (1.0 - al - be)  # 无条件方差起步
    e_prev = 0.0
    for t in range(T):
        h_t = om + al * e_prev**2 + be * h_prev
        e_t = np.sqrt(h_t) * rng.standard_normal()
        r[t] = mu + lam * np.sqrt(h_t) + e_t
        h[t] = h_t
        e_prev, h_prev = e_t, h_t
    # 预烧 100 弱化初值影响
    return r[100:]


def main() -> None:
    y = make_gm_data()
    np.savetxt(os.path.join(HERE, "gm_smoke_data.csv"), y, delimiter=",",
               header="r", comments="")
    out = {}
    for form in ("vol", "var", "log"):
        mod = ARCHInMean(y, volatility=GARCH(p=1, q=1),
                         distribution=Normal(), rescale=False, form=form)
        res = mod.fit(disp="off")  # cov_type 默认 'robust' → std_errors = BW QMLE
        out[form] = {
            "params": {k: float(v) for k, v in
                       zip(res.params.index, res.params.values)},
            "llf": float(res.loglikelihood),
            "nobs": int(res.nobs),
            "robust_se": {k: float(v) for k, v in
                          zip(res.params.index, res.std_err)},
            "converged": int(res.convergence_flag) == 0,
        }
        print(f"[{form}] llf={out[form]['llf']:.6f} params={out[form]['params']}")
        print(f"      robust_se={out[form]['robust_se']}")
    path = os.path.join(HERE, "gm_arch_baselines.json")
    json.dump(out, open(path, "w"), indent=1)
    print("baselines ->", path)
    gen_inc(y, out)


def gen_inc(y: np.ndarray, out: dict) -> None:
    """生成 tests/unit/timeseries/gm_baseline.inc (数据数组 + 三 form 基线)

    容差注记 (7B test_garch_model.cpp 头注先例, 参数可达精度):
      - params/llf: 1e-4 名义 (实测 1e-5~1e-6; scipy SLSQP 解析梯度 vs
        自研数值梯度的落点差, spec 名义 1e-8~1e-10 仅同优化器可达)
      - robust SE: 5e-2 名义 (实测 ≤2e-2; 数值 Hessian/OPG 步长舍入噪声)
    """
    vals = [repr(float(v)) for v in y]
    rows = ["    " + ", ".join(vals[i:i + 4]) + ","
            for i in range(0, len(vals), 4)]

    def case(form, name):
        c = out[form]
        p, se = c["params"], c["robust_se"]
        return "\n".join([
            f"constexpr GmCase {name}{{  // form='{form}'",  # noqa: F541
            f"    {repr(p['Const'])}, {repr(p['kappa'])}, {repr(p['omega'])},",
            f"    {repr(p['alpha[1]'])}, {repr(p['beta[1]'])},",
            f"    {repr(c['llf'])},",
            f"    {repr(se['Const'])}, {repr(se['kappa'])}}};",
        ])

    inc = "\n".join([
        "// 自动生成: tests/fixtures/timeseries/verify_gm.py (arch 8.0.0)",
        "//   数据: RandomState(42), T=900 (vol 形 DGP λ=0.50, 预烧 100)",
        "//   基准: ARCHInMean(rescale=False, Normal, GARCH(1,1)) 三 form",
        "// 勿手改 — 重新生成请运行该脚本",
        "#pragma once",
        '#include "cpphub/core/types.hpp"',
        "namespace cpphub { inline namespace v1 { namespace timeseries {",
        "namespace garch { namespace gm_baseline {",
        "",
        "constexpr Size T = 900;",
        "",
        "constexpr Real Y_GM[T] = {",
        *rows,
        "};",
        "",
        "// (mu, lambda, omega, alpha, beta, llf, se_mu, se_lambda)",
        "struct GmCase { Real mu, lambda, omega, alpha, beta, llf, se_mu,",
        "                se_lambda; };",
        case("vol", "ARCH_VOL"),
        case("var", "ARCH_VAR"),
        case("log", "ARCH_LOG"),
        "",
        "}}}}}  // namespace cpphub::v1::timeseries::garch::gm_baseline",
        "",
    ])
    with open(INC_PATH, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(inc)
    print("inc ->", INC_PATH)


if __name__ == "__main__":
    main()
