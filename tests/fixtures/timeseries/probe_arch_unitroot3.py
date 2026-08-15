# probe_arch_unitroot3.py - 精确区段源码核对 (M2 排幻觉, 第三轮)
# dump: _df_select_lags / ADF._compute_statistic / DFGLS._compute_statistic /
#       PP._compute_statistic / KPSS._compute_statistic + _autolag /
#       VR._compute_statistic / mackinnonp / mackinnoncrit / kpss_crit
import inspect
import arch.unitroot.unitroot as u

SRC = inspect.getsourcefile(u).replace("\\", "/")
lines = open(SRC, encoding="utf-8").read().splitlines()


def show(lo, hi, label):
    print(f"===== {label} (unitroot.py:{lo}-{hi}) =====")
    for i in range(lo - 1, min(hi, len(lines))):
        print(f"{i+1}: {lines[i]}")


show(372, 447, "_df_select_lags")
show(773, 819, "ADF._select_lag + _compute_statistic")
show(936, 997, "DFGLS._compute_statistic")
show(1121, 1187, "PP._compute_statistic")
show(1309, 1377, "KPSS._compute_statistic + _autolag")
show(1704, 1767, "VR._compute_statistic")
show(1768, 1928, "mackinnonp + mackinnoncrit")
show(1929, 1966, "kpss_crit")
