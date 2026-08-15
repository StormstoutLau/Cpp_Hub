# probe_arch_unitroot2.py - 精确区段源码核对 (M2 排幻觉, 第二轮)
import inspect
import arch.unitroot.unitroot as u

SRC = inspect.getsourcefile(u).replace("\\", "/")
lines = open(SRC, encoding="utf-8").read().splitlines()


def show(lo, hi, label):
    print(f"===== {label} (unitroot.py:{lo}-{hi}) =====")
    for i in range(lo - 1, min(hi, len(lines))):
        print(f"{i+1}: {lines[i]}")


# 找到所有关键函数/类的行号
for i, l in enumerate(lines):
    s = l.strip()
    if s.startswith("def ") or s.startswith("class "):
        print(f"L{i+1}: {s[:80]}")
