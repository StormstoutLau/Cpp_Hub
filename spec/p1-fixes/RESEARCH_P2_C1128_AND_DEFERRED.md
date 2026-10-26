# P2 遗留专项调研: C1128 /bigobj + 三个数值契约延迟项

> **调研日期**: 2026-10-26 · **方法**: 官方文档/社区案例检索（WebSearch + WebFetch 逐页）+ 本地 MSVC 实证（逐 TU 编译复现 + `/bigobj` 探针构建）+ mpmath 60 位高精度锚
> **上游**: [EXTERNAL_REVIEW_2026-10-20.md](../../docs/audit/EXTERNAL_REVIEW_2026-10-20.md) §二；[RESEARCH_P2_RISKS.md](./RESEARCH_P2_RISKS.md) §4.5/§4.6/§5.6；[DEVELOPMENT_LOG.md](../../docs/DEVELOPMENT_LOG.md)「P2 数值契约批」
> **证据分级**（沿用 [ASSERTION_EVIDENCE_FRAMEWORK.md](../../docs/ASSERTION_EVIDENCE_FRAMEWORK.md)）: **A** = 一手来源逐字/本地原始输出；**B** = 本地实测或推理链；**C** = 判断类（仅 rationale）

## 0. 结论速览

| 项 | 复核结论 | 推荐方案 | 变更面 | 优先级 |
|---|---|---|---|---|
| C1128 `/bigobj` | **5 个直接 TU**（原记录"4 TU"漏 `test_tool_registry`）；**2 个真正依赖目标**（`test_c_api` / `test_c_api_tools`）；仅 Debug 触发 | MSVC 全局 `/bigobj`（一处）或 5 目标定向 `/bigobj` | 根 `CMakeLists.txt` 1 行（全局）或 helper + 5 处调用 | **P0** |
| `midas_diagnostics:285` | `1 - regularized_lower_gamma` = χ² 上尾 `1-P`，深尾归零同 RISK-036 | 迁移到 `chi2_sf(stat, df)` | 1 文件 1 行 + 注释 | P1 |
| SMM 深尾 | `1 - [2Φ(CD)-1]^m` 在 `CD=10, m=4` 输出 **0**（真值 6.1e-23）；`CD=8` 相对误差 ~7% | `detail::smm_joint_p` = `-expm1(m·log1p(-q))`, `q=2Φ(-CD)` | 1 头文件 + 1 测试 | **P2** |
| 关联缺陷 D-1 | 单 horizon 与 `z2_p_value` 走 `2-2Φ`，同样深尾归零（`vr_two_sided_p:153`） | 改 `2*normal_cdf(-|stat|)` | 1 行 + 1 测试 | P2a（建议同批） |
| Chow-Denning note | 正相关时独立 SMM 近似 p 偏大 → 检验保守（欠拒绝） | note 增"正相关 → p 偏保守"说明 | registry note 1 行 + **golden skeleton REGEN** | P3 |

---

## 1. C1128 `/bigobj` 专项

### 1.1 事实复核（A，本地原始输出）

Debug 构建逐 TU 复现（MSVC 19.x / VS18 v180，`build` 目录）：

| TU | 无 `/bigobj` | 加 `/bigobj`（探针 `build_probe_bigobj`） |
|---|---|---|
| `src/c_api/c_api.cpp`（`cpphub_c_api`） | `c_api.cpp(1,1): error C1128` | `cpphub_c_api.lib` 生成 |
| `tests/unit/timeseries/test_integration_phase7c.cpp` | `test_integration_phase7c.cpp(1,1): error C1128` | `.exe` 生成 |
| `tests/unit/timeseries/test_integration_phase8.cpp` | `test_integration_phase8.cpp(1,1): error C1128` | `.exe` 生成 |
| `tests/unit/timeseries/test_vecm_model.cpp` | `test_vecm_model.cpp(1,1): error C1128` | `.exe` 生成 |
| `tests/unit/c_api/test_tool_registry.cpp` | `test_tool_registry.cpp(1,1): error C1128` | `.exe` 生成 |
| `tests/unit/c_api/test_c_api.cpp` | 依赖阻断（`cpphub_c_api` 失败） | `.exe` 生成；**删掉自身 `/bigobj` 重编仍过** |
| `tests/unit/c_api/test_c_api_tools.cpp` | 依赖阻断 | `.exe` 生成；**删掉自身 `/bigobj` 重编仍过** |

完整错误文本（本地，中文 locale）：

```
error C1128: 节数超过对象文件格式限制: 请使用 /bigobj 进行编译
```

**对既有记录的更正（B）**: 原记录"4 个巨型 TU + 3 个 c_api 依赖目标"不准确。实测 `test_tool_registry` **不链接** `cpphub_c_api`，它自身就是第 5 个直接触发 C1128 的 TU；真正因依赖 `cpphub_c_api` 而失败的是 `test_c_api` / `test_c_api_tools`。该更正应同步回 `DEVELOPMENT_LOG.md` 与 `project_memory.md` 的批条目（见 §6）。

Debug 旗标（`build/tests/*.vcxproj` + `CMakeCache.txt`，A）：

- `CMAKE_CXX_FLAGS_DEBUG = /Ob0 /Od /RTC1`
- VS 生成器 Debug: `<Optimization>Disabled</Optimization>`、`<BasicRuntimeChecks>EnableFastChecks</BasicRuntimeChecks>`、`<RuntimeLibrary>MultiThreadedDebugDLL</RuntimeLibrary>`、`<DebugInformationFormat>ProgramDatabase</DebugInformationFormat>`（即 `/Zi`，非 `/ZI`）
- Release: `/O2 /Ob2 /DNDEBUG`，无调试节 → **同一 TU 在 Release 不触发**（CI 主线 `CMAKE_BUILD_TYPE: Release` 因此全绿）

### 1.2 官方文档与社区证据（A）

**C1128 官方说明** — [Fatal Error C1128](https://learn.microsoft.com/en-us/cpp/error-messages/compiler-errors-1/fatal-error-c1128?view=msvc-170)：

- "An .obj file exceeded the number of allowable sections, a COFF object file format limitation."
- "Reaching this section limitation can be the result of using /Gy and a debug build; /Gy causes functions to go into their own COMDAT sections. In a debug build, there is a debug info section for each COMDAT function."
- "C1128 can also be caused when there are too many inline functions."
- 修法三选一："divide your source file into multiple source code files, compile without /Gy, or compile with /bigobj."
- x64 放大机制："x64 will have at least 4 sections associated with each function compiled /Gy or inlined from templates or class-inline: code, pdata, and debug info, and possibly xdata."

**`/bigobj` 官方说明** — [/bigobj](https://learn.microsoft.com/en-us/cpp/build/reference/bigobj-increase-number-of-sections-in-dot-obj-file?view=msvc-170)：

- 默认上限 **65,279**（"almost 2^16"）addressable sections；`/bigobj` 扩到 **4,294,967,296（2^32）**。
- "machine-generated code, or code that makes heavy use of template libraries, may require .obj files that can hold more sections."
- "`/bigobj` is enabled by default on Universal Windows Platform (UWP) projects"（微软自家把重模板/生成代码场景默认打开）。

**与调试信息兼容性** — [Supported Code Changes (C++)](https://learn.microsoft.com/en-us/visualstudio/debugger/supported-code-changes-cpp?view=vs-2022)：

- "Projects that use the MSVC v120 toolset and the C/C++ `/bigobj` switch. Edit and Continue with `/bigobj` is only supported in the MSVC v140 and later build tools."
- 本仓库用 VS18（v180）+ `/Zi`，`/bigobj` 无兼容性障碍。

**社区/上游先例（A）**：

| 来源 | 做法 | 与本仓库关系 |
|---|---|---|
| 本仓 `third_party/eigen/cmake/EigenTesting.cmake:77` | `ei_add_target_property(${targetname} COMPILE_FLAGS "/bigobj")` | Eigen 自家测试在 MSVC 下**定向**加 `/bigobj`；与本库同型（重模板测试） |
| [LadybugDB CMakeLists.txt](https://raw.githubusercontent.com/LadybugDB/ladybug/master/CMakeLists.txt) | `add_compile_options("/bigobj")`，注释："Large test objects ... exceed the default 64k section limit; /bigobj raises it (fatal error C1128 otherwise)" | 现代 CMake 全局做法先例 |
| [NITRC NIREP #5589](https://www.nitrc.org/tracker/?aid=5589&atid=1644&func=detail&group_id=435) | `ADD_DEFINITIONS(/bigobj)` / `CMAKE_CXX_FLAGS += /bigobj` | 旧式全局做法；说明该问题长期普遍存在 |
| [Ceres Conan patch](http://raw.githubusercontent.com/conan-io/conan-center-index/master/recipes/ceres-solver/all/patches/2.2.0-0002-fix-mingw-build.patch) | MinGW 用 `-Wa,-mbig-obj`（等价物） | 本仓库不产 MinGW，**不适用**；GCC/ELF 无 65k 节限制 |

**本仓库既有教训（B）**：`$<$<CONFIG:Debug>:...>` 等 genex 在 `add_compile_options` 中被 VS 生成器**静默丢弃**（D8016 回归链实测，见 DEVELOPMENT_LOG RISK-029/030/031 批）。因此本项**不要**用 `$<$<CONFIG:Debug>:/bigobj>` 只在 Debug 生效；应全配置启用。

### 1.3 方案对比（C）

| 方案 | pros | cons | 结论 |
|---|---|---|---|
| **A. MSVC CXX 全局 `/bigobj`** | 1 行；覆盖当前/未来所有 TU；官方推荐；社区主流 | 所有 MSVC C++ 目标（含 gtest）目标文件格式变化；面最大 | **推荐** |
| B. 5 个目标定向 `/bigobj` | 变更面最小；Eigen 同型做法 | 需人工枚举；**已实测原清单漏第 5 个**，测试增长后会反复漏；新重 TU 复发 | 可选（若要求最小面） |
| C. 拆分巨型 TU | 根治节数；不依赖编译选项 | 测试文件是完整集成场景，拆分会破坏组织；`c_api.cpp` 拆分涉及静态初始化/ODR；改动大、回归风险高 | 不推荐（长期可评估） |
| D. `compile without /Gy` | 官方三选一 | 本仓 Debug 未启用 `/Gy`（vcxproj 无 `FunctionLevelLinking`）；Release `/O2` 隐含 `/Gy` 但无调试节；不是本因 | 无效 |
| E. `/Z7` 替代 `/Zi` | 可能减少节 | 改变调试信息格式/增量链接；不保证解决；无必要 | 不推荐 |

### 1.4 推荐（B）

**方案 A：根 `CMakeLists.txt` 的 MSVC 块给 CXX 全配置加 `/bigobj`**：

```cmake
if(MSVC)
    add_compile_options(
        $<$<COMPILE_LANGUAGE:CXX>:/utf-8>
        $<$<COMPILE_LANGUAGE:CXX>:/W3>
        # C1128 根治: COFF .obj 默认 65,279 节上限; Debug (/Zi) 下每个 COMDAT
        # 函数额外带调试节, 本库 5 个巨型 TU (c_api + VECM + 2 集成 +
        # tool_registry) 实测超限。 /bigobj 仅把节容量扩到 2^32, 不改代码
        # 生成/ABI; VS18 (v140+) 与 /Zi 及 /ZI 兼容。 全配置启用 (非仅 Debug):
        # genex CONFIG 在 VS 生成器下有静默丢弃先例, 且配置间格式一致无收益;
        # 用 COMPILE_LANGUAGE genex 与现有 /utf-8 同款, 避免传给 nvcc。
        $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()
```

要点：

- **不传入 INTERFACE 库**：`add_compile_options` 只影响本工程当前目录及之后 `add_subdirectory` 创建的目标，不会写进 `cpphub_core` 的导出接口，消费者不受影响。
- **不在 GCC 分支加**：ELF 无 COFF 65k 限制；不要用 `-Wa,-mbig-obj`（仅 MinGW 需要）。
- **CUDA 兼容**：`/bigobj` 必须包在 `$<COMPILE_LANGUAGE:CXX>` 内，避免 nvcc 收到未知选项。**genex 落点已实测**（最小 CMake 工程 + VS18 生成器）：Debug/Release/MinSizeRel/RelWithDebInfo 四个配置的 `AdditionalOptions` 均为 `%(AdditionalOptions) /utf-8 /bigobj`。D8016 事故记录的是 `$<CONFIG:>`/嵌套 genex 被丢弃；本方案只用简单 `$<COMPILE_LANGUAGE:CXX>`，但仍按既有教训在实施后检查 vcxproj 落点。
- **配置一致性**：Release 也带 `/bigobj`，保证同一 TU 在不同配置下目标文件格式一致；对运行期/链接产物零影响。
- 如坚持最小面，使用 helper：
  ```cmake
  function(cpphub_msvc_bigobj target)
      if(MSVC)
          target_compile_options(${target} PRIVATE /bigobj)
      endif()
  endfunction()
  ```
  调用点：`src/c_api/CMakeLists.txt:9` 后 + `tests/CMakeLists.txt` 的 `test_tool_registry`(:215) / `test_vecm_model`(:735) / `test_integration_phase7c`(:751) / `test_integration_phase8`(:813) 后。

**可选 CI 加固（建议但非必需）**：现有 CI 只有 Release，C1128 永远看不见。若要把 Debug 编译纳入门禁，新增一个 **Windows Debug build-only** job（不跑全量 ctest，避免小时级）：

```yaml
  windows-debug-guard:
    name: Debug compile guard (MSVC)
    runs-on: windows-2022
    steps:
      - uses: actions/checkout@v4
        with:
          submodules: recursive
      - uses: ilammy/msvc-dev-cmd@v1
        with:
          arch: x64
      - name: Configure (Debug, C API ON)
        shell: bash
        run: cmake -S . -B build_debug -G "Visual Studio 17 2022" -A x64 -DCPPHUB_ENABLE_C_API=ON
      - name: Build heavy TUs (Debug)
        shell: bash
        run: cmake --build build_debug --config Debug --parallel 4 --target cpphub_c_api test_tool_registry test_c_api test_c_api_tools test_vecm_model test_integration_phase7c test_integration_phase8
```

成本估算（本地探针实测）：configure ~9s + 5 个重 TU + 2 依赖目标 Debug 全编 ~4–6 min；加 checkout/setup 约 8–10 min。收益：把"C1128 类"与 D 族 Debug-only 缺陷（RISK-029/030/031 批的 `Size` 下溢族）纳入门禁。

### 1.5 验证（B，已执行）

1. 复现：`cmake --build build --config Debug --target <5 个 TU 之一>` → C1128（已逐一捕获）。
2. 探针：`cmake -S . -B build_probe_bigobj -G "Visual Studio 18 2026" -DCPPHUB_ENABLE_C_API=ON -DCMAKE_CXX_FLAGS_DEBUG="/bigobj /Ob0 /Od /RTC1"` → 7 个目标 Debug 全部生成（探针目录为临时产物，验完删除）。
3. 依赖目标自证：在探针中删除 `test_c_api` / `test_c_api_tools` 自身 vcxproj 的 `/bigobj` 后重编 → 仍通过（证明它们**不需要**自身 `/bigobj`，只需 `cpphub_c_api` 成功）。
4. 实施后验收：`cmake --build build --config Debug --target cpphub_c_api test_tool_registry test_c_api test_c_api_tools test_vecm_model test_integration_phase7c test_integration_phase8` → 0 C1128；Release CI 维持 9/9。

### 1.6 陷阱（B）

- **数量更正**：5 直接 TU + 2 依赖目标（不是 4 + 3）。
- **`test_tool_registry` 不是依赖受害者**：它只链 `cpphub_core + cpphub_timeseries_mat`，自身即超限。
- **`/bigobj` 不是"掩盖 bug"**：C1128 是 COFF 格式上限，不是代码错误；官方修复手段即 `/bigobj`/拆文件/去 `/Gy`。
- **不要只对 Debug 生效**：`$<CONFIG:>` genex 在 VS 生成器下有静默丢弃先例（本仓 D8016），且配置间格式不一致无收益；`$<COMPILE_LANGUAGE:CXX>` 已单独验证落点。
- **Release CI 盲区**：不新增 Debug job 的话，后续新重 TU 仍可能在开发者本地先炸。
- **目标文件格式**：`/bigobj` 产物由同一 MSVC 链接器消费，无 ABI/运行期差异；若下游有只认标准 COFF 的旧工具，才需退回定向方案。

---

## 2. `midas_diagnostics.hpp:285` 迁移 `chi2_sf`

### 2.1 事实与变更（A）

现状 [midas_diagnostics.hpp](../../include/cpphub/timeseries/midas/midas_diagnostics.hpp#L282-L287)：

```cpp
if (out.df == 0) {
    out.p = std::numeric_limits<Real>::quiet_NaN();
} else {
    out.p = 1.0 - regularized_lower_gamma(
                       static_cast<Real>(out.df) / 2.0, out.stat / 2.0);
}
```

`1 - P` 正是 RISK-036 已判定并修复的模式（[math.hpp](../../include/cpphub/core/math.hpp#L150-L165) 的上尾重写 + `chi2_sf`）。改为：

```cpp
} else {
    out.p = chi2_sf(out.stat, static_cast<Real>(out.df));
}
```

同时把文件头 `:17` 的公式注释 `p = 1 − χ²cdf(stat, df)` 改为 `p = chi2_sf(stat, df)`（χ² 上尾）。

### 2.2 语义与风险（B）

- 中尾（`stat=2.108, df=2`）：`x/2 = 1.054 < a+1 = 2`，`chi2_sf` 内部仍走 `1-P` 支路，**数值逐位一致**（不是"数 ULP 差"）；现有 R 基准锚 `stat=2.1080208, p=0.3485372`（容差 1e-4）不变。
- 深尾：旧路径在 `p < ~1e-16` 归零；新路径保留到 underflow。
- 边界：`chi2_sf(x<=0) = 1.0`；旧 `regularized_lower_gamma` 对负 `x` **抛异常**。HAH 统计量是半正定二次型 `h0'(I−PΔP')h0`，理论非负；若浮点舍入出极小负值，新路径返回 1（不拒绝）而非让诊断整体抛错——更稳健。
- `df == 0` 的 NaN 守卫保留（`chi2_sf` 要求 `k>0`）。

### 2.3 测试（B）

- 既有 `MidasModel.HahTestMatchesR`（[test_midas_model.cpp:246-252](../../tests/unit/timeseries/test_midas_model.cpp#L246-L252)）覆盖 `stat/df/p` 与 R 基准；迁移后应保持绿。
- 深尾能力由 `chi2_sf` 的 4 个新锚（RISK-036 批）保证，不重复造锚。
- 结论：**等价迁移 + 回归锚保护**；不新增行为测试。

---

## 3. SMM 联合 p 值深尾精度

### 3.1 现状与缺陷（A/B）

[variance_ratio_test.hpp:246-251](../../include/cpphub/timeseries/unit_root/variance_ratio_test.hpp#L246-L251)：

```cpp
const Real one_sided = 2.0 * detail::normal_cdf(cd_abs) - 1.0;
Real powm = 1.0;
for (Size i = 0; i < k_list.size(); ++i) powm *= one_sided;
r.chow_denning_p_value = 1.0 - powm;
```

两处抵消：`2Φ(CD)-1`（上尾信息被 1 吃掉）+ `1-powm`（`p` 极小时归零）。本地 mpmath 60 位对照（B）：

| m | CD | 真值 `1-(1-q)^m` | 现实现（double） | 稳健式（double） |
|---|---|---|---|---|
| 1 | 8 | 1.2441921148543568247e-15 | 1.3322676295501878e-15（+7%） | 1.2441921148543639e-15 |
| 4 | 8 | 4.9767684594174180107e-15 | 5.329070518200751e-15（+7%） | 4.976768459417446e-15 |
| 4 | 10 | 6.0958824193284208528e-23 | **0（归零）** | 6.095882419328476e-23 |
| 4 | 3 | 0.010755529521374500943 | 0.010755529521374596 | 0.010755529521374515 |

### 3.2 文献与社区（A）

- [cppreference expm1](https://kurzy.kpi.fei.tuke.sk/c-reference/en/c/numeric/math/expm1.html)：expm1 "is more accurate than the expression exp(arg)-1.0 if arg is close to zero"；Notes 明示 `(1+x)^n - 1 = expm1(n * log1p(x))` 用于金融复利——与本式同构。
- [cppreference log1p](https://kurzy.kpi.fei.tuke.sk/c-reference/en/c/numeric/math/log1p.html)：log1p "is more precise than log(1+arg) if arg is close to zero"。
- [Boost.Math Complements](https://www.boost.org/doc/libs/1_84_0/libs/math/doc/html/math_toolkit/stat_tut/overview/complements.html)：`1 - p` 在 `p≈1` 时"cancellation error will cause you to lose accuracy, perhaps totally"。
- [Goldberg 1991 §2.4 Cancellation](https://docs.oracle.com/cd/E19422-01/819-3693/ncg_goldberg.html)；scipy `sf`："also defined as 1-cdf, but sf is sometimes more accurate"（[RESEARCH_P2_RISKS §5.3](./RESEARCH_P2_RISKS.md#53-社区案例a)）。
- 公式来源：Chow & Denning (1993), DOI [10.1016/0304-4076(93)90051-6](https://doi.org/10.1016/0304-4076(93)90051-6)；R `vrtest::Chow.Denning` 源码代数反演（[RESEARCH_P2_RISKS §4.2](./RESEARCH_P2_RISKS.md#42-文献与实现a)）。

### 3.3 推荐实现（B）

在 `variance_ratio_test.hpp` 的 `detail` 增加单一来源：

```cpp
/// SMM(m,∞) 联合 p 值: p = 1 - [2Φ(|CD|)-1]^m = 1 - (1-q)^m,
/// q = 2Φ(-|CD|) = erfc(|CD|/√2) (直接上尾, 不经过 1-Φ)。
/// 深尾用 -expm1(m·log1p(-q)): cd=10,m=4 真值 6.1e-23, 朴素式在 double 下归零。
inline Real smm_joint_p(Real cd_abs, Size m) {
    if (m == 0) return 1.0;
    const Real q = 2.0 * normal_cdf(-cd_abs);  // = erfc(|CD|/√2)
    if (q >= 1.0) return 1.0;                  // cd = 0
    return -std::expm1(static_cast<Real>(m) * std::log1p(-q));
}
```

调用改成：

```cpp
r.chow_denning_p_value = detail::smm_joint_p(cd_abs, k_list.size());
```

边界：`cd=0 → q=1 → 1.0`；`q` 下溢为 0 → `-expm1(0) = 0`；`m` 为 horizon 个数（非 `max(k)`，§4.6 陷阱已守）。

### 3.4 关联缺陷 D-1（建议同批）

单 horizon 路径与 `z1/z2/clm` p 值走 [vr_two_sided_p:152-154](../../include/cpphub/timeseries/unit_root/variance_ratio_test.hpp#L152-L154)：

```cpp
inline Real vr_two_sided_p(Real stat) noexcept {
    return 2.0 - 2.0 * normal_cdf(std::fabs(stat));
}
```

`2 - 2Φ(10)` 在 double 下 = **0**（真值 1.5239706048321052132e-23）；`|z|=8` 时相对误差 ~7%。建议同批改为：

```cpp
inline Real vr_two_sided_p(Real stat) noexcept {
    return 2.0 * normal_cdf(-std::fabs(stat));  // = erfc(|stat|/√2)
}
```

理由：MCP `cpphub_variance_ratio_test` 同时返回 `z1/z2/clm` 的 p 值；只修多 horizon 联合 p 会留下同型不一致。若严格限定原延迟项范围，可将 D-1 单独排期，但需在日志显式登记。

### 3.5 测试锚（B）

新增 `VarianceRatioTest.SmmJointPDeepTail`（TDD 先红：引用 `detail::smm_joint_p` 会 C3861 编译红）：

| 输入 | 期望（mpmath 60dps） | 容差 |
|---|---|---|
| `(0, 4)` | 1.0 | 精确 |
| `(3, 4)` | 0.010755529521374500943 | 相对 1e-12 |
| `(8, 1)` | 1.2441921148543568247e-15 | 相对 1e-12 |
| `(8, 4)` | 4.9767684594174180107e-15 | 相对 1e-12 |
| `(10, 4)` | 6.0958824193284208528e-23 | 相对 1e-12；并断言 `> 1e-24`（朴素式为 0） |

D-1 同法加 `detail::vr_two_sided_p(10.0)` 锚（真值 1.5239706048321052132e-23，旧实现 0 → 红）。

既有 `MultiChowDenning*` 用例的 `p_ref` 是 `1-(2Φ-1)^2`（中尾，容差 1e-8），迁移后仍绿；可选择顺手把 oracle 改成稳健式，但非必需。

Golden 影响：中尾 p 变化 ≤ 1 ULP，`DEFAULT_TOL rtol=1e-9` 内，**无需更新数值**。

---

## 4. Chow-Denning note 保守性说明

### 4.1 事实（A/B）

- 公式是 Šidák / SMM(m,∞) 独立近似；Chow & Denning 原文有限样本用 SMM(m,T)；各 `Z2(k)` 因 horizon 重叠而**正相关**（[RESEARCH_P2_RISKS §4.2](./RESEARCH_P2_RISKS.md#42-文献与实现a)）。
- 正相关使 `P(max|Z2| ≤ x)` 大于独立情形（极端同相关 ρ=1 时 `Φ(x) > Φ(x)^m`），故按独立式算出的 p **偏大** → 实际拒绝率低于名义水平 = **保守（欠拒绝）**，不是虚高显著性。
- R `vrtest` 不返回 p 值（只给临界值），因此该 note 是库内唯一的语义提醒。

### 4.2 变更（B）

[cpphub_tool_registry.hpp:543-547](../../include/cpphub/c_api/cpphub_tool_registry.hpp#L543-L547) 的 note：

```cpp
{"多 horizon 联合 SMM(m,inf) (U16: 多重检验校正; horizon 间正相关 → p 偏保守)",
 horizon_note}
```

同时在 [variance_ratio_test.hpp:205-209](../../include/cpphub/timeseries/unit_root/variance_ratio_test.hpp#L205-L209) 的注释补一句"正相关时保守"（代码注释不影响契约）。

### 4.3 Golden 影响（A）

`python/tests/fixtures/mcp_golden.json` 的 `skeleton` 是**含 note 原文的 byte-exact 字符串**（`_normalize` 只把数字替换为 `<num>`，不改文本）。改 note 必改 skeleton：

1. 重建 nanobind `_core`（与 P2 数值契约批同法）；
2. `CPPHUB_MCP_REGEN_GOLDEN=1 pytest python/tests/test_mcp_tools.py -k write_golden`；
3. 检查 golden diff **只包含该 note 文本**（其余数字/结构不应变）；
4. 取消环境变量后跑全量 `pytest python/tests/test_mcp_tools.py` + `test_version.py`（基线 111 passed, 1 skipped）。

不推荐手改长转义 skeleton（易引入 JSON 转义错误；上一批已发生 golden 尾逗号事故）。

### 4.4 C++ 侧测试（B）

在 `ToolRegistry.VarianceRatioSmoke`（[test_tool_registry.cpp:245-257](../../tests/unit/c_api/test_tool_registry.cpp#L245-L257)）增加断言：找到 `method == "chow_denning_joint"` 的行，其 `notes` 任一元素含 `"保守"`。先加断言 → 红；再改 note → 绿。这同时防止 note 被后续重构悄悄删除。

---

## 5. 实施顺序与验收

> 遵循 TDD 与"最小改动"；P0 可独立提交，P1–P3 可同批。

| 阶段 | 变更 | 红证据 | 绿证据 |
|---|---|---|---|
| **P0** C1128 | 根 CMake 加 `/bigobj` | 5 TU 各自 `error C1128`（已捕获） | 7 目标 Debug 全生成；Release CI 9/9 |
| **P1** midas | `chi2_sf` 迁移 | 等价迁移，无行为红；以 R 锚回归保护 | `MidasModel.HahTestMatchesR` 绿 |
| **P2** SMM | `detail::smm_joint_p` + 5 锚 | 新测试引用未实现函数 → C3861 编译红 | 5 锚相对 1e-12；`cd=10` 不再归零 |
| **P2a** D-1 | `vr_two_sided_p` 上尾化 + 锚 | `vr_two_sided_p(10)` 期望 1.52e-23，旧值 0 → FAIL | 锚绿；既有 p 用例不回退 |
| **P3** note | registry note + C++ 断言 + golden REGEN | C++ note 断言红；pytest skeleton 红 | C++ 绿；pytest 111 passed 1 skipped |

**范围套件验收**：

```powershell
# C++ (Release 或 Debug 均可; C1128 专项需 Debug)
ctest --test-dir build -C Debug -R "VarianceRatioTest|MidasModel" --output-on-failure
cmake --build build --config Debug --target test_tool_registry test_c_api test_c_api_tools
# Python
& "C:\Users\Peng\AppData\Local\Programs\Python\Python312\python.exe" -m pytest python/tests -q
```

**全量权威门**：CI Release 全量 + 可选新增 Debug build-only job。

---

## 6. 风险与遗留

1. **`/bigobj` 全局化**：只影响 COFF 节容量，无 ABI/运行期差异。若担心第三方目标（gtest）格式变化，可选定向方案；但需接受"清单漏项"风险（本专项已实证一次）。
2. **CI Debug job 成本**：约 8–10 min/次；若不加，C1128 类问题只能靠本地 Debug 构建发现。
3. **文档更正**：`DEVELOPMENT_LOG.md` 与 `project_memory.md` 中"4 个巨型 TU + 3 个依赖目标"应改为"5 个直接 TU + 2 个依赖目标（`test_tool_registry` 为第 5 个直接成员）"。建议在 P0 提交说明中一并更正。
4. **同型扫描（本轮未修）**：`pss_critical_values.hpp:70` 的 `1.0 - normal_cdf(poly)` 也属深尾 `1-Φ` 模式；PSS 边限 p 值在极端统计量时精度受损。建议登记为后续批（不在本专项三项内）。
5. **golden REGEN 依赖重建扩展**：若 CI/本地未重建 `_core`，note 修改不会反映到 pytest，导致 skeleton 断言假绿/假红；按 §4.3 顺序执行。
6. **D-1 范围**：若决定不同批修 `vr_two_sided_p`，必须在日志中显式写"单 horizon 深尾仍归零"，避免审计误判为遗漏。

---

## 7. 参考

### 官方文档（A）

- [Fatal Error C1128](https://learn.microsoft.com/en-us/cpp/error-messages/compiler-errors-1/fatal-error-c1128?view=msvc-170)
- [/bigobj (Increase Number of Sections in .Obj file)](https://learn.microsoft.com/en-us/cpp/build/reference/bigobj-increase-number-of-sections-in-dot-obj-file?view=msvc-170)
- [Supported Code Changes (C++) — /bigobj & Edit and Continue](https://learn.microsoft.com/en-us/visualstudio/debugger/supported-code-changes-cpp?view=vs-2022)
- [cppreference expm1](https://kurzy.kpi.fei.tuke.sk/c-reference/en/c/numeric/math/expm1.html) / [log1p](https://kurzy.kpi.fei.tuke.sk/c-reference/en/c/numeric/math/log1p.html)（镜像；原站 en.cppreference.com 当日抓取失败）
- [Boost.Math Complements](https://www.boost.org/doc/libs/1_84_0/libs/math/doc/html/math_toolkit/stat_tut/overview/complements.html)

### 社区/上游先例（A）

- [LadybugDB CMakeLists — 全局 /bigobj + C1128 注释](https://raw.githubusercontent.com/LadybugDB/ladybug/master/CMakeLists.txt)
- [NITRC NIREP #5589 — /bigobj 全局修复案例](https://www.nitrc.org/tracker/?aid=5589&atid=1644&func=detail&group_id=435)
- [Ceres Conan MinGW `-Wa,-mbig-obj` 补丁](http://raw.githubusercontent.com/conan-io/conan-center-index/master/recipes/ceres-solver/all/patches/2.2.0-0002-fix-mingw-build.patch)
- 本仓 `third_party/eigen/cmake/EigenTesting.cmake:77`（Eigen 测试定向 `/bigobj`）

### 本仓既有文档（A/B）

- [RESEARCH_P2_RISKS.md](./RESEARCH_P2_RISKS.md) §4.2/§4.5/§4.6（Chow-Denning 公式、golden skeleton 不变、expm1 建议、保守性 note）
- [RESEARCH_P2_RISKS.md](./RESEARCH_P2_RISKS.md) §5.2/§5.3/§5.6（χ² 上尾、Boost/Goldberg/scipy、midas 迁移）
- [DEVELOPMENT_LOG.md](../../docs/DEVELOPMENT_LOG.md)「P2 数值契约批: RISK-035/036」（RISK-036 实施与 C1128 既有性裁定）
