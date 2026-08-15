# probe_arch_unitroot5.py - _select_best_ic / _autolag_ols (AIC/BIC 精确公式)
import inspect
import arch.unitroot.unitroot as u

SRC = inspect.getsourcefile(u).replace("\\", "/")
lines = open(SRC, encoding="utf-8").read().splitlines()


def show(lo, hi, label):
    print(f"===== {label} (unitroot.py:{lo}-{hi}) =====")
    for i in range(lo - 1, min(hi, len(lines))):
        print(f"{i+1}: {lines[i]}")


show(147, 204, "_select_best_ic")
show(299, 371, "_autolag_ols")
