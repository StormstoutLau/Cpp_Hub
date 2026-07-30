# Cpp_Hub

> A header-first C++20 quantitative pricing library — BS / Heston / PDE / Tree / MC / AAD Greeks / VaR / SVI, with Python bindings and optional CUDA MC.
>
> 一个以头文件为主的 C++20 量化定价库：解析解 / 蒙特卡洛 / PDE / 树形 / AAD Greeks / VaR / SVI 波动率曲面，附带 Python 绑定与可选 CUDA MC 内核。

---

## Status

| Stage | Scope | Tests |
|---|---|---|
| **Phase 1–3** | Core / PayOff / BS / Heston / MC / Sobol QMC / PDE / Tree / AAD-Lite Greeks / Pathwise & LR Greeks / Greeks Factory / VaR & Backtesting / Calibration / SVI | 286 / 286 ✅ |
| **Phase 4 LITE** | SSVI 跨期限 / Python 绑定 (nanobind) / GPU MC (CUDA, 主控站可选) | 320 / 320 ✅ |
| **Cross-platform** | MSVC (Win10) + GCC (Ubuntu NEX / GTR-Pro) 三平台一致 | 286 / 286 位精确一致 |

跨平台浮点确定性由 `-ffp-contract=off` (GCC) + `/fp:precise` (MSVC) 保证 IEEE-754 严格模式。

---

## Features

- **Models**: GBM, Heston (Heston 1993 / Albrecher 2007 "Little Trap" with branch-cut fix)
- **Pricing engines**: Analytic (BSM, Heston CF), Monte Carlo (path + Sobol QMC + Brownian Bridge), PDE (Crank–Nicolson PSOR), Tree (Leisen–Reimer binomial, trinomial)
- **Greeks**: Analytic, Pathwise, Likelihood Ratio (LR), AAD (autodiff reverse mode), Factory auto-dispatch
- **Risk**: MC VaR / ES, Kupiec POF backtesting
- **Volatility surface**: SVI / SSVI (Gatheral-Jacquier, no-arb butterfly / calendar conditions)
- **Calibration**: LM / Nelder–Mead / DE optimizers
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
│   └── performance/gpu/    # gpu_mc (CUDA optional)
├── src/                    # Non-header sources (CUDA kernels + CPU stubs)
├── tests/                  # Unit + validation tests (GoogleTest)
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
