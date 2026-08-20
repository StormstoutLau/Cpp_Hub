# -*- coding: utf-8 -*-
# gen_phase7c_granger_baseline.py - M1 Granger 测试基准 .inc 生成 (Phase 7C)
#
# 数据源 (verify_granger.py 生成, 固化):
#   granger_smoke_data.csv        (gen_granger_fixture.py, T=250)
#   coint_smoke_data.csv          (gen_coint_fixture.py, T=250, I(1) 对)
#   granger_statsmodels_baselines.txt (statsmodels 0.14.4:
#     gr_* = grangercausalitytests 4 统计量;
#     grty_* = OLS f_test 增广 Wald (TY);
#     grhac_* = 显式 NW 三明治 Wald + statsmodels HAC f_test 交叉)
#
# 输出: tests/unit/timeseries/granger_baseline.inc
#   命名: GR_<TAG>_P<p> (10 值) / GRTY_<TAG>_P<k>_D<d> (6 值)
#         / GRHAC_<TAG>_P<p>_L<L> (6 值) / X, Y, Z, Y1, Y2 (夹具列)

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "..", "unit", "timeseries"))


def fmt(v):
    v = float(v)
    if v != v:  # NaN → 合法 constexpr 表达式
        return "NQ"
    return repr(v)


def array(name, vals):
    body = ",\n    ".join(fmt(v) for v in vals)
    return "inline constexpr double %s[] = {\n    %s,\n};\n" % (name, body)


def load(path):
    d = {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2 and parts[1].strip():
            d[parts[0]] = [float(v) for v in parts[1].split()]
    return d


def read_csv_cols(name, header):
    cols = []
    with open(os.path.join(HERE, name)) as fh:
        for row in csv.reader(fh):
            if row and row[0] != header:
                cols.append([float(v) for v in row])
    return list(map(list, zip(*cols)))  # 行 → 列


def main():
    base = load(os.path.join(HERE, "granger_statsmodels_baselines.txt"))
    gcols = read_csv_cols("granger_smoke_data.csv", "x")
    ccols = read_csv_cols("coint_smoke_data.csv", "y1")
    gx, gy, gz = gcols[0], gcols[1], gcols[2]
    cy = ccols

    out = []
    out.append("// =============================================================================\n")
    out.append("// granger_baseline.inc - M1 Granger 因果检验测试基准 (Phase 7C v1.7)\n")
    out.append("//\n")
    out.append("// 生成: tests/fixtures/timeseries/gen_phase7c_granger_baseline.py (勿手改)\n")
    out.append("// 源: granger_statsmodels_baselines.txt (statsmodels 0.14.4)\n")
    out.append("//   gr_*   = grangercausalitytests (4 统计量, 1e-10)\n")
    out.append("//   grty_* = TY 增广 Wald (OLS f_test, df=k, 1e-8)\n")
    out.append("//   grhac_* = NW 三明治 Wald (cov_hac_simple 无修正, 1e-8)\n")
    out.append("// 数据: granger_smoke_data.csv (x→y 真因果) + coint_smoke_data.csv (I(1) 对)\n")
    out.append("// =============================================================================\n\n")
    out.append("#include <limits>\n\n")
    out.append("namespace cpphub {\ninline namespace v1 {\nnamespace timeseries {\n")
    out.append("namespace granger_baseline {\n\n")
    out.append("// NaN 编码 (§1.4-5)\n")
    out.append("inline constexpr double NQ =\n")
    out.append("    std::numeric_limits<double>::quiet_NaN();\n\n")

    out.append(array("X", gx))
    out.append(array("Y", gy))
    out.append(array("Z", gz))
    out.append(array("Y1", cy[0]))
    out.append(array("Y2", cy[1]))

    for tag, vals in base.items():
        out.append(array(tag.upper(), vals))

    out.append("}  // namespace granger_baseline\n")
    out.append("}  // namespace timeseries\n")
    out.append("}  // inline namespace v1\n")
    out.append("}  // namespace cpphub\n")

    path = os.path.join(OUT, "granger_baseline.inc")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(out))
    print("wrote", path, "(%d arrays)" % len(base))


if __name__ == "__main__":
    main()
