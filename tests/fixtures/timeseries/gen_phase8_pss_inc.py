# gen_phase8_pss_inc.py — PSS 响应面系数表 .inc 生成器
# Phase 8 v1.8 M1 Step 1 (PHASE8_SPEC.md v1.1 §2.4 Step 2; probe_pss_cv.py 布局冻结)
#
# 来源: 本地 statsmodels 0.14.6 tsa/ardl/pss_critical_values.py (BSD-3)
#   四张表各 100 键 (k ∈ 1..10, case ∈ 1..5, I1 ∈ {False,True}):
#     large_p   — 4 系数 (3 阶多项式, 升幂), stat ≤ stat_star 时使用
#     small_p   — 3 系数 (2 阶多项式, 升幂), stat > stat_star 时使用
#     stat_star — 1 标量 (阈值, 原始 stat 尺度比较)
#     crit_vals — 4 档临界值 {90%, 95%, 99%, 99.9%} (表常量直查, 非多项式求值)
#   转录总量: 100×4 + 100×3 + 100 + 100×4 = 1200 数字
#
# 求值语义 (probe_pss_cv.py [B] 冻结, 数值自证 ≤1e-16):
#   p = 1 − Φ(c₀ + c₁·log(stat) + c₂·log(stat)² + c₃·log(stat)³)
#
# 输出: include/cpphub/timeseries/ardl/pss_response_surface.inc
#   扁平索引: idx = (k−1)·10 + (case−1)·2 + (i1 ? 1 : 0)  — 与 PSS10 裁决一致,
#   接口的 k 即 force 变量数 (PSS 语义; statsmodels bounds_test 调用层偏移 +1 不镜像)
#
# 重新生成: & C:\Users\Peng\.conda\envs\open-webui\python.exe gen_phase8_pss_inc.py
import sys
from pathlib import Path

from statsmodels.tsa.ardl import pss_critical_values as pcv

OUT = Path(__file__).resolve().parents[3] / "include" / "cpphub" / "timeseries" / "ardl" / "pss_response_surface.inc"


def idx(k: int, case: int, i1: bool) -> int:
    return (k - 1) * 10 + (case - 1) * 2 + (1 if i1 else 0)


def fmt(v: float) -> str:
    return repr(float(v))  # 最短往返表示, 保证逐位一致


def collect(table: dict, width: int) -> list:
    rows = [None] * 100
    for (k, case, i1), vals in table.items():
        assert len(vals) == width, f"宽度异常 key={(k, case, i1)}: {vals}"
        rows[idx(k, case, i1)] = vals
    assert all(r is not None for r in rows), "键空间不完整 (100 键)"
    return rows


def emit_rows(out, name: str, rows: list, per_row: int) -> None:
    out.write(f"inline constexpr Real {name}[][{per_row}] = {{\n")
    for r in rows:
        out.write("    {" + ", ".join(fmt(v) for v in r) + "},\n")
    out.write("};\n\n")


def main() -> None:
    large = collect(pcv.large_p, 4)
    small = collect(pcv.small_p, 3)
    star = [None] * 100
    for (k, case, i1), v in pcv.stat_star.items():
        star[idx(k, case, i1)] = v
    assert all(s is not None for s in star)
    crit = collect(pcv.crit_vals, 4)

    lines = []
    lines.append("// ============================================================================")
    lines.append("// pss_response_surface.inc — PSS 边限检验响应面系数表 (自动生成, 勿手改)")
    lines.append("//")
    lines.append("// Phase 8 v1.8 M1 (PHASE8_SPEC.md v1.1 §2.4; 布局/求值语义 probe_pss_cv.py 冻结)")
    lines.append("//")
    lines.append("// 生成脚本: tests/fixtures/timeseries/gen_phase8_pss_inc.py")
    lines.append("// 来源: statsmodels 0.14.6 tsa/ardl/pss_critical_values.py (BSD-3)")
    lines.append("//   逐键逐值机器转录 (非手抄); 全量 diff 由 probe_pss_cv.py 对照验收")
    lines.append("//")
    lines.append("// 键空间: k ∈ 1..10 (force 变量数, PSS 语义 — PSS10: statsmodels bounds_test")
    lines.append("//   调用层用 nvar = 1+k 偏移查表, 本库不镜像) × case ∈ 1..5 × I1 ∈ {0,1}")
    lines.append("// 扁平索引: idx = (k−1)·10 + (case−1)·2 + i1")
    lines.append("//")
    lines.append("// 求值语义 (probe [B] 数值自证 ≤1e-16):")
    lines.append("//   临界值 = PSS_CRIT_VALS[idx] 直查 ({90%, 95%, 99%, 99.9%} 四档, PSS9:")
    lines.append("//     PSS 原表为 {10,5,2.5,1}%, statsmodels 无 2.5% 档而增 0.1% 档)")
    lines.append("//   p 值 = 1 − Φ(Σ_{i=0..n} cᵢ·log(stat)ⁱ); 系数按 stat > PSS_STAT_STAR[idx]")
    lines.append("//     ? PSS_SMALL_P (n=2) : PSS_LARGE_P (n=3) 切换 (阈值用原始 stat 比较)")
    lines.append("//")
    lines.append("// 锚点 (PSS 2001 Table CI(iii) k=1 5% ≈ 4.94/5.73, Statalist Kripfganz 旁证):")
    lines.append("//   crit_vals[(1,3,False)][1] = 4.9228967, crit_vals[(1,3,True)][1] = 5.7217442")
    lines.append("//   — pss_critical_values.hpp static_assert 双锚 (编译期)")
    lines.append("//")
    lines.append("// 本文件被 pss_critical_values.hpp 于 namespace detail 内直接 #include")
    lines.append("// ============================================================================")
    lines.append("#pragma once")
    lines.append("")
    header = "\n".join(lines)

    with OUT.open("w", newline="\n") as f:
        f.write(header + "\n\n")
        emit_rows(f, "PSS_LARGE_P", large, 4)
        emit_rows(f, "PSS_SMALL_P", small, 3)
        f.write("inline constexpr Real PSS_STAT_STAR[] = {\n")
        for s in star:
            f.write(f"    {fmt(s)},\n")
        f.write("};\n\n")
        emit_rows(f, "PSS_CRIT_VALS", crit, 4)

    n = sum(len(r) for r in large) + sum(len(r) for r in small) + len(star) + sum(len(r) for r in crit)
    print(f"OK: {OUT}")
    print(f"  键 100 × 4 表; 数字总数 = {n} (预期 1200: 400+300+100+400)")
    print(f"  双锚自检: (1,3,False)[1] = {fmt(crit[idx(1,3,False)][1])} (预期 4.9228967)")
    print(f"            (1,3,True)[1]  = {fmt(crit[idx(1,3,True)][1])} (预期 5.7217442)")
    assert fmt(crit[idx(1, 3, False)][1]) == "4.9228967"
    assert fmt(crit[idx(1, 3, True)][1]) == "5.7217442"
    assert n == 1200
    print("  断言全部通过")


if __name__ == "__main__":
    sys.exit(main())
