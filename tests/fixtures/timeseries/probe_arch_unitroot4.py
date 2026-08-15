# probe_arch_unitroot4.py - _estimate_df_regression + cov_nw + ADF/VR 默认参数
import inspect
import arch.unitroot.unitroot as u

SRC = inspect.getsourcefile(u).replace("\\", "/")
lines = open(SRC, encoding="utf-8").read().splitlines()


def show(lo, hi, label):
    print(f"===== {label} (unitroot.py:{lo}-{hi}) =====")
    for i in range(lo - 1, min(hi, len(lines))):
        print(f"{i+1}: {lines[i]}")


show(1, 60, "imports")
show(454, 491, "_estimate_df_regression")
show(754, 772, "ADF.__init__")
show(903, 935, "DFGLS.__init__")
show(1095, 1110, "PP.__init__")
show(1278, 1298, "KPSS.__init__")
show(1644, 1703, "VarianceRatio.__init__")

# cov_nw 定义
import arch.utility.cov as cov
print("===== cov_nw =====")
print(inspect.getsource(cov.cov_nw))
