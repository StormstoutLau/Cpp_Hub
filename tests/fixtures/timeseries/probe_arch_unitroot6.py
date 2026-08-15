# probe_arch_unitroot6.py - ADF._compute_statistic / KPSS._autolag / DFGLS._compute_statistic 源码
import inspect

import arch.unitroot.unitroot as u

SRC = inspect.getsourcefile(u).replace("\\", "/")
lines = open(SRC, encoding="utf-8").read().splitlines()


def show(lo, hi, label):
    print(f"===== {label} (unitroot.py:{lo}-{hi}) =====")
    for i in range(lo - 1, min(hi, len(lines))):
        print(f"{i+1}: {lines[i]}")


# mackinnonp / mackinnoncrit / kpss_crit 源码
import inspect as _ins
from arch.unitroot import critical_values as _cv_pkg  # noqa: F401

import arch.unitroot.unitroot as uu

m = _ins.getsource(uu.mackinnonp)
print("===== mackinnonp =====")
print(m[:4200])
k = _ins.getsource(uu.kpss_crit)
print("===== kpss_crit =====")
print(k)
