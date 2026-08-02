# Cpp_Hub

> A header-first C++20 quantitative pricing library — BS / Heston / PDE / Tree / MC / AAD Greeks / VaR / SVI, with Python bindings and optional CUDA MC.
>
> 一个以头文件为主的 C++20 量化定价库：解析解 / 蒙特卡洛 / PDE / 树形 / AAD Greeks / VaR / SVI 波动率曲面，附带 Python 绑定与可选 CUDA MC 内核。
>
> v1.4.0 起新增高频计量经济学模块 (HFE), 对标 R `highfrequency` 1.0.3 (Boudt, Kleen, Sjørup 2022, JSS doi:10.18637/jss.v104.i08), 提供 C++ 工程级 RV/BPV/RSV/BNS 跳跃检验实现. v1.4.1 新增 Realized Kernel (BNS 2008 ECTA) 微结构噪声稳健估计, 含 12 种核函数 + 噪声方差 + 最优带宽选择.

---

## Status

| Stage | Scope | Tests |
|---|---|---|
| **Phase 1–3** | Core / PayOff / BS / Heston / MC / Sobol QMC / PDE / Tree / AAD-Lite Greeks / Pathwise & LR Greeks / Greeks Factory / VaR & Backtesting / Calibration / SVI | 286 / 286 ✅ |
| **Phase 4 LITE** | SSVI 跨期限 / Python 绑定 (nanobind) / GPU MC (CUDA, 主控站可选) | 320 / 320 ✅ |
| **Phase 5 v1.4.0** | HFE: TAQ Reader / Realized Measures (RV/RVol/RQ/BPV/RSV) / BNS Jump Test | 18 / 18 ✅ (R baseline 对标, spec 矩阵 15 + R baseline exact 3) |
| **Phase 5 v1.4.1** | HFE: Realized Kernel (BNS 2008) / 12 核函数 / 噪声方差 ω² / 最优带宽 H* | 14 / 14 ✅ (R baseline 对标, 11 核函数 + 3 R baseline) |
| **Cross-platform** | MSVC (Win10) + GCC (Ubuntu NEX / GTR-Pro) 三平台一致 | 1300 / 1300 位精确一致 ✅ (v1.4.0+v1.4.1 三平台验证通过) |

跨平台浮点确定性由 `-ffp-contract=off` (GCC) + `/fp:precise` (MSVC) 保证 IEEE-754 严格模式。

---

## Features

- **Models**: GBM, Heston (Heston 1993 / Albrecher 2007 "Little Trap" with branch-cut fix)
- **Pricing engines**: Analytic (BSM, Heston CF), Monte Carlo (path + Sobol QMC + Brownian Bridge), PDE (Crank–Nicolson PSOR), Tree (Leisen–Reimer binomial, trinomial)
- **Greeks**: Analytic, Pathwise, Likelihood Ratio (LR), AAD (autodiff reverse mode), Factory auto-dispatch
- **Risk**: MC VaR / ES, Kupiec POF backtesting
- **Volatility surface**: SVI / SSVI (Gatheral-Jacquier, no-arb butterfly / calendar conditions)
- **Calibration**: LM / Nelder–Mead / DE optimizers
- **High-frequency econometrics (v1.4.0+)**: TAQ CSV reader + time aggregation, Realized Variance / Volatility / Quarticity, Bipower Variation (jump-robust), Realized Semivariance (±), BNS Jump Test (BN-S 2006, with Truncated Power Quarticity), 多资产 Realized Covariance — 对标 R `highfrequency` 1.0.3, 容差 1e-12
- **High-frequency econometrics (v1.4.1+)**: Realized Kernel (BNS 2008 ECTA) — 12 种核函数 (Rectangular/Bartlett/Second/Epanechnikov/Cubic/Fifth/Sixth/Seventh/Eighth/Parzen/TukeyHanning/ModifiedTukeyHanning), 噪声方差 ω²=RV/(2n) (H-L 2006), 最优带宽 H*=c·ξ^(4/5)·(ω²/IV)^(2/5)·n^(3/5) (BNS 2008 eq.51), 严格对标 R `rKernelCov` (KK() + kernelEstimator() 源码实测, 容差 1e-12)
- **Python bindings**: `import cpphub` via nanobind (BSM / Heston / MC / Greeks / VaR)
- **GPU MC** (optional): CUDA Philox4x64-10 RNG, bit-exact match with CPU `cpphub::core::rng`

---

## Project Layout

```
Cpp_Hub/
├── include/cpphub/         # Header-only library
│   ├── core/               # config, constants, error, linalg, math, rng, simd, types
│   ├── models/             # gbm, svi
│   ├── pricing/            # analytic, monte_carlo, pde, tree
│   ├── risk/               # greeks (aad/pathwise/lr), var
│   ├── calibration/        # calibrator, objective, optimizer
│   ├── hfecon/             # v1.4.0+ 高频计量经济学 (对标 R highfrequency 1.0.3)
│   │   ├── data/           # TAQ reader (CSV, 时间聚合)
│   │   ├── measures/       # RV / RVol / RQ / BPV / RSV / rCov + (v1.4.1) 12 核函数 + Realized Kernel
│   │   ├── tests/          # BNS Jump Test (BN-S 2006, TPQ)
│   │   ├── noise/          # (v1.4.1) 噪声方差 ω² + 最优带宽 H* (BNS 2008)
│   │   ├── models/         # (v1.4.2) HAR, HEAVY
│   │   └── liquidity/      # (v1.4.3) 流动性度量
│   └── performance/gpu/    # gpu_mc (CUDA optional)
├── src/                    # Non-header sources (CUDA kernels + CPU stubs)
├── tests/                  # Unit + validation tests (GoogleTest)
│   └── fixtures/hfe/       # R baseline 生成脚本 + baselines.json (CI gate)
├── benchmarks/             # Performance benchmarks
├── python/                 # Python bindings (nanobind) + pytest
├── third_party/autodiff/   # Vendored autodiff headers (forward/reverse mode)
├── docs/                   # Architecture / Phase specs / ADR / Audit checklist
└── CMakeLists.txt
```

---

## Build

### Requirements

- C++20 compiler (MSVC 19.x, GCC 13+, Clang 14+)
- CMake ≥ 3.25
- Python ≥ 3.9 (optional, for Python bindings)
- CUDA ≥ 12 (optional, for GPU MC)

### Configure & Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release

# Run tests
ctest --test-dir build --build-config Release --output-on-failure
```

### Python bindings

```bash
pip install -e python/
pytest python/tests/
```

### Optional: GPU MC

```bash
cmake -S . -B build_cuda -DCMAKE_BUILD_TYPE=Release -DCPPHUB_ENABLE_CUDA=ON
cmake --build build_cuda --config Release
```

Without CUDA, a CPU stub is compiled automatically and 15 GPU MC tests still run on CPU.

---

## Numerical Benchmarks

| Module | Target | Achieved |
|---|---|---|
| BSM analytic vs benchmark | 1e-12 | ✅ |
| Heston CF vs Schoutens table | 1e-8 | ✅ (after branch-cut fix) |
| MC convergence order | -0.5 ± 0.05 | ✅ |
| Sobol QMC variance reduction | ≥ 10x | ✅ |
| AAD Greeks vs analytic | 1e-10 | ✅ |
| Pathwise / LR / Analytic / AAD Greeks cross-check | 1e-6 | ✅ |
| Cross-platform float determinism | bit-exact | ✅ (MSVC + 2× GCC) |
| HFE Realized Measures vs R `highfrequency` 1.0.3 | 1e-12 | ✅ (RV/RVol/RQ/BPV/RSV, 9 case baseline) |
| HFE BNS Jump Test vs R `BNSjumpTest` | 1e-10 | ✅ (Z statistic + p-value, no-jump + jump 双向) |
| HFE Realized Kernel vs R `rKernelCov` | 1e-12 | ✅ (BNS 2008, 12 核函数, 17 case B1/B2 baseline) |
| HFE 核函数解析值 vs R `KK()` 源码 | 1e-15 | ✅ (12 核函数, R 源码 realizedMeasures.cpp L16-74) |

---

## Documentation

Authoritative docs live under `docs/`:

- `docs/architecture/ARCHITECTURE.md` — layered architecture, module contracts
- `docs/phases/phase{1-4}/PHASEx_SPEC.md` — phase execution specs
- `docs/decisions/ADR_INDEX.md` — architecture decision records
- `docs/audit/AUDIT_CHECKLIST.md` — weighted audit checklist
- `BUILD_PLAN.md` — 12-week / 4-phase roadmap
- `TRACEABILITY_REPORT.md` — source traceability (textbook / open-source, zero hallucination)

---

## Methodology

This library is built test-driven (TDD) with a strict **form-first → numerical falsification → retrospective formalization** loop:

1. Write the theorem / formula first
2. Implement and verify numerically against textbook tables (Joshi, Duffy, Schoutens)
3. Cross-validate across compilers and platforms; bit-exact determinism required
4. Log deviations from literature as "discoveries" (potential research output)

Cross-platform verification is treated as a **standard-compliance test**, not just a "compiles & runs" test — MSVC's permissiveness masks several standard-compliance issues that GCC strictly rejects.

---

## License

Apache License 2.0 — see [LICENSE](LICENSE).

Third-party headers in `third_party/autodiff/` are vendored from [autodiff](https://github.com/autodiff/autodiff) (MIT License).
