# -*- coding: utf-8 -*-
"""probe_jidt_te.py — H2/风险⑭ 关闭探针: 本库 TE 组装层 vs JIDT 运行时对照.

目标 (PHASE9_SPEC §10.1 H2 = 风险⑭):
  1. 单位行为: JIDT 连续 TE 计算器输出单位实证 (调研 A47 "bits 默认" 的真伪)
  2. k/l/delay 组装语义: 本库 TransferEntropyConfig{source_history, target_history,
     delay} 的嵌入 [X_{t-d},...,X_{t-kd}, Y_t, Y_{t-d},...,Y_{t-ld}] 与 JIDT
     TransferEntropyCalculatorKraskov 的 (k, k_tau, l, l_tau, delay) 嵌入映射
  3. 残余偏差量化: 同数据同参数同 K 下 |te_cpp - te_jidt|

JIDT 1.6.1 属性语义 (源码实证, tests/fixtures/tmp/jidt 发行包):
  - initialise(k, k_tau, l, l_tau, delay) 五参重载:
      k = 目标历史长度, k_tau = 目标嵌入延迟, l = 源历史长度,
      l_tau = 源嵌入延迟, delay = 源-目标延迟 (源锚定 t-delay)
  - 裸属性 "k" = Kraskov K (经 CMI 计算器); "k_HISTORY"/"l_HISTORY" 为历史长度
  - "NOISE_LEVEL_TO_ADD" 默认 1e-8 (加噪破坏平手) → 置 0 以确定性对齐本库 jitter=0
  - "DYN_CORR_EXCL" (Theiler 窗) 默认 0 = 与本库"仅排除自身"语义一致
  - 结构映射: delay=1 时本库嵌入与 JIDT (k=kHist, k_tau=1, l=lHist, l_tau=1, delay=1)
    逐点同构; delay>1 时 JIDT 目标历史锚定 t-1 (硬编码), 本库锚定 t-delay →
    结构不同构, 只做行为记录 (见 CONFIGS_DELAY2)

前置 (一次性):
  - pip install jpype1
  - JRE: https://api.adoptium.net/v3/binary/latest/21/ga/windows/x64/jre/hotspot/normal/eclipse
    解压至 tmp/jidt/jdk/ (或设 JIDT_JAVA_HOME 指向任意 JDK/JRE)
  - jar: https://github.com/jlizier/jidt/releases/download/v1.6.1/infodynamics-dist-1.6.1.zip
    解压至 tmp/jidt/dist/infodynamics.jar (或设 JIDT_JAR)
  - C++ 驱动 (scratch): tmp/h2_te_driver.exe — 构建命令见 DRIVER_BUILD_HINT
"""
import os
import subprocess
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
JAR = os.environ.get("JIDT_JAR", os.path.join(REPO, "tmp", "jidt", "dist", "infodynamics.jar"))
JAVA_HOME = os.environ.get("JIDT_JAVA_HOME", "")
DRIVER = os.path.join(REPO, "tmp", "h2_te_driver.exe")
DATA_DIR = os.path.join(REPO, "tmp", "h2")

DRIVER_BUILD_HINT = (
    'cl /nologo /O2 /EHsc /std:c++20 /utf-8 /Iinclude '
    '/Ithird_party\\random123\\include tmp\\h2_te_driver.cpp /Fe:tmp\\h2_te_driver.exe'
)

PHI, A, N = 0.5, 0.3, 2000
GC_HALF_NATS = 0.5 * np.log(1.0 + A * A / (1.0 - PHI * PHI))  # 0.056664 nats

# (k=source_history, l=target_history, delay, knn) — delay=1 为结构精确可对照区
CONFIGS_DELAY1 = [
    (1, 1, 1, 3),  # 组装基线 (knn=3 与本库默认一致)
    (2, 1, 1, 3),  # 源历史深度
    (1, 2, 1, 3),  # 目标历史深度
    (3, 2, 1, 3),  # 双边历史
    (2, 2, 1, 3),
    (1, 1, 1, 4),  # knn=4 (JIDT kraskov 默认 K)
]
# delay=2: JIDT 目标历史锚定 t-1 与本库锚定 t-delay 结构不同构 — 仅行为记录
CONFIG_DELAY2 = (1, 1, 2, 3)
SEEDS = [3000, 3001]


def driver(args):
    out = subprocess.run([DRIVER] + args, capture_output=True, text=True, check=True)
    return out.stdout.strip()


def driver_gen(seed, path):
    txt = driver(["gen", str(N), str(PHI), str(A), str(seed), path])
    assert txt.startswith("gen_ok"), txt


def driver_te(path, k, l, d, knn):
    txt = driver(["te", path, str(k), str(l), str(d), str(knn), "1"])
    kv = dict(p.split("=", 1) for p in txt.split())
    return float(kv["te"]), float(kv["bits"]), int(kv["n_eff"])


def zscore(v):
    v = np.asarray(v, dtype=float)
    return (v - v.mean()) / v.std(ddof=1)


def main():
    if not os.path.exists(JAR):
        sys.exit(f"JIDT jar 缺失: {JAR} (见文件头前置说明)")
    if not os.path.exists(DRIVER):
        sys.exit(f"C++ 驱动缺失: {DRIVER}\n构建: {DRIVER_BUILD_HINT}")
    if JAVA_HOME:
        os.environ["JAVA_HOME"] = JAVA_HOME
    elif not os.environ.get("JAVA_HOME"):
        jdk = os.path.join(REPO, "tmp", "jidt", "jdk")
        if os.path.isdir(jdk):
            subs = [os.path.join(jdk, d) for d in os.listdir(jdk)]
            subs = [p for p in subs if os.path.isdir(p)]
            if len(subs) == 1:
                os.environ["JAVA_HOME"] = subs[0]

    import jpype

    jpype.startJVM(classpath=[JAR])
    JDouble = jpype.JArray(jpype.JDouble)
    TECalc = jpype.JClass(
        "infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov")

    def jidt_te(xs, ys, k, l, d, knn):
        """JIDT Kraskov TE, 与本库 (k,l,d) 组装对齐.

        delay=1: initialise(目标历史=本库l, k_tau=1, 源历史=本库k, l_tau=1, delay=1)
        与本库嵌入逐点同构 (见文件头结构映射); 裸 "k" 属性 = Kraskov K.
        """
        c = TECalc()
        c.setProperty("k", str(knn))                # 转发至 CMI 计算器 = Kraskov K
        c.setProperty("NOISE_LEVEL_TO_ADD", "0")    # 关默认 1e-8 加噪 (确定性)
        c.setProperty("DYN_CORR_EXCL", "0")         # Theiler=0 = 本库排除自身语义
        c.initialise(l, 1, k, 1, 1 if d == 1 else d)  # (kHist, k_tau, lHist, l_tau, delay)
        c.setObservations(JDouble(xs.tolist()), JDouble(ys.tolist()))
        return float(c.computeAverageLocalOfObservations())

    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"解析锚: GC/2 = {GC_HALF_NATS:.6f} nats (phi={PHI}, a={A}); "
          f"JIDT 1.6.1 Kraskov, NOISE=0, DYN_CORR_EXCL=0")
    print(f"{'config(k,l,d,K)':>18} {'seed':>5} {'cpp_te':>10} {'jidt_te':>10} "
          f"{'|diff|':>10} {'cpp/GC2':>8} {'jidt/GC2':>9}")

    worst = 0.0
    rows = []
    for seed in SEEDS:
        path = os.path.join(DATA_DIR, f"h2_var_{seed}.txt")
        driver_gen(seed, path)
        raw = np.loadtxt(path)
        xs, ys = zscore(raw[:, 0]), zscore(raw[:, 1])
        for (k, l, d, knn) in CONFIGS_DELAY1:
            cpp_te, cpp_bits, n_eff = driver_te(path, k, l, d, knn)
            j_te = jidt_te(xs, ys, k, l, d, knn)
            diff = abs(cpp_te - j_te)
            worst = max(worst, diff)
            rows.append((f"({k},{l},{d},K{knn})", seed, cpp_te, j_te, diff))
            print(f"{f'({k},{l},{d},K{knn})':>18} {seed:>5} {cpp_te:>10.6f} {j_te:>10.6f} "
                  f"{diff:>10.6f} {cpp_te / GC_HALF_NATS:>8.3f} {j_te / GC_HALF_NATS:>9.3f}")
        # delay=2 行为记录: JIDT 源锚定 t-2 但目标历史锚定 t-1 (结构不同构)
        k, l, d, knn = CONFIG_DELAY2
        cpp_te, _, n_eff = driver_te(path, k, l, d, knn)
        j_te = jidt_te(xs, ys, k, l, d, knn)
        print(f"{f'({k},{l},{d},K{knn})':>18} {seed:>5} {cpp_te:>10.6f} {j_te:>10.6f} "
              f"{abs(cpp_te - j_te):>10.6f}  <- delay=2 结构不同构, 仅记录")

    print(f"\n[VERDICT] delay=1 结构精确区: max|cpp - jidt| = {worst:.6f} nat "
          f"({len(rows)} 组)")
    jpype.shutdownJVM()


if __name__ == "__main__":
    main()
