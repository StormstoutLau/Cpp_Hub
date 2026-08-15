# probe_arch_unitroot.py - arch 8.0.0 unitroot 源码约定探查 (M2 排幻觉)
# 用法: python probe_arch_unitroot.py
# 输出: Schwert 规则位置 / PP Z(tau) / KPSS _autolag / VR / MacKinnon CV 关键源码行
import inspect
import arch.unitroot.unitroot as u
import arch.unitroot.critical_values.dickey_fuller as df_cv
import arch.unitroot.critical_values.dfgls as dfgls_cv

SRC = inspect.getsourcefile(u).replace("\\", "/")
lines = open(SRC, encoding="utf-8").read().splitlines()


def show(lo, hi, label):
    print(f"===== {label} ({SRC.split('/')[-1]}:{lo}-{hi}) =====")
    for i in range(lo - 1, min(hi, len(lines))):
        print(f"{i+1}: {lines[i]}")


print("### Schwert 规则出现位置 ###")
for i, l in enumerate(lines):
    if "12.0" in l and "0.25" in l:
        print(f"L{i+1}: {l.strip()}")

print("\n### ADF 默认参数与 lag 选择 ###")
for i, l in enumerate(lines):
    if "class ADF" in l:
        show(i + 1, i + 80, "ADF class")
        break

print("\n### PP 类 (Z(tau) 修正) ###")
for i, l in enumerate(lines):
    if "class PhillipsPerron" in l:
        show(i + 1, i + 115, "PhillipsPerron class")
        break
