// =============================================================================
// conduction_metrics.hpp - 信息论事前度量 (v2.0+ 预留)
//
// Phase 7A 附录 C: v2.0+ 预留接口空文件
//
// Scope 声明 (PHASE7A_FALSIFICATION_SPEC.md §Scope):
//   - Phase 7A 严格聚焦事后证伪统计量 (post-estimation)
//   - 信息论事前度量 (L0 层) 属于 v2.0+ scope
//   - 当前版本不实现, 仅预留接口位置以保持文件结构完整
//
// 预留接口 (v2.0+ 实现):
//   - ConductionMatrix { MatrixXD M; Real tau; Size r_eff; }
//   - compute_conduction_strength(model, param_grid) -> Real tau
//   - fisher_mapping(M) -> MatrixXD I(theta)
//   - effective_rank(M, threshold) -> Size
//
// 定位 (双层诊断体系):
//   L0 事前 (信息论度量, v2.0+):
//     - tau (传导强度): 模型选择层, 跨模型可比
//     - r_eff (有效秩): 独立脆弱性通道数
//     - v_max (最脆弱方向): 参数失效定位
//
//   L1/L2/L3 事后 (v1.6 Phase 7A 已实现):
//     - 计量: JB/LB/BG/BP/White (residual_diagnostics.hpp)
//     - 波动率: 标准化残差/z^2 LB (volatility_diagnostics.hpp)
//     - 高频: HAR 残差 LB (hfecon_diagnostics.hpp)
//     - 风险: DQ/Berkowitz/ES (risk_diagnostics.hpp)
//     - 定价: IV 拟合优度 (pricing_diagnostics.hpp)
//     - Greeks: 跨方法一致性 (greeks_consistency.hpp)
//
// 参考:
//   - INFORMATION_THEORY_METRICS_RESEARCH.md (v2.0+ 调研报告)
//   - ADR-015 方案 B (Eigen3 隔离原则, v2.0+ 将依赖 linalg_dynamic.hpp)
// =============================================================================
#pragma once

#include "cpphub/core/types.hpp"
// v2.0+ 启用: #include "cpphub/core/linalg_dynamic.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// v2.0+ 信息论事前度量 (见 INFORMATION_THEORY_METRICS_RESEARCH.md)
// 当前版本 (v1.6) 不实现, 仅预留接口位置
//
// 计划接口:
//   struct ConductionMatrix {
//       MatrixXD M;       // 传导矩阵
//       Real tau;         // 传导强度
//       Size r_eff;       // 有效秩
//   };
//   Real compute_conduction_strength(model, param_grid);
//   MatrixXD fisher_mapping(const MatrixXD& M);
//   Size effective_rank(const MatrixXD& M, Real threshold);

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
