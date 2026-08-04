// SOURCE: PHASE6_IMPLEMENTATION_PLAN §2.3 - 双层 linalg 动态尺寸层 (ADR-013)
// 版本: Eigen 3.4.0 (header-only, MPL2.0, C++14 兼容稳定版)
// 用途: econometrics/ 模块 (回归/SVD/QR/稀疏)
// 隔离: pricing/models/risk 模块不依赖此文件, 保持 core/linalg.hpp 固定尺寸
// 排幻觉点: LLT 适用范围窄于 LU, 仅用于 SPD 矩阵 (Eigen 官方基准测试文档)
#pragma once

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <stdexcept>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace linalg::dynamic {

// =============================================================================
// MatrixXD - 动态尺寸矩阵 (封装 Eigen::MatrixXd)
// 设计: 不暴露 Eigen 类型到公共 API (除 eigen() 显式转换)
// =============================================================================
class MatrixXD {
public:
    MatrixXD() = default;
    explicit MatrixXD(Size rows, Size cols) : data_(rows, cols) {}

    // 从 Eigen 构造 (显式, 避免隐式泄漏)
    explicit MatrixXD(const Eigen::MatrixXd& m) : data_(m) {}

    // 元素访问
    Real& operator()(Size i, Size j) { return data_(i, j); }
    Real operator()(Size i, Size j) const { return data_(i, j); }

    // 维度
    Size rows() const { return static_cast<Size>(data_.rows()); }
    Size cols() const { return static_cast<Size>(data_.cols()); }

    // 数据指针 (与 C API 互操作)
    Real* data() { return data_.data(); }
    const Real* data() const { return data_.data(); }

    // 显式转换 (边界明确, 不隐式泄漏 Eigen)
    const Eigen::MatrixXd& eigen() const { return data_; }
    Eigen::MatrixXd& eigen() { return data_; }

    // 矩阵运算
    MatrixXD operator*(const MatrixXD& rhs) const {
        return MatrixXD(data_ * rhs.data_);
    }
    MatrixXD operator+(const MatrixXD& rhs) const {
        return MatrixXD(data_ + rhs.data_);
    }
    MatrixXD operator-(const MatrixXD& rhs) const {
        return MatrixXD(data_ - rhs.data_);
    }
    MatrixXD transpose() const {
        return MatrixXD(data_.transpose());
    }

private:
    Eigen::MatrixXd data_;
};

// =============================================================================
// VectorXD - 动态尺寸向量 (封装 Eigen::VectorXd)
// =============================================================================
class VectorXD {
public:
    VectorXD() = default;
    explicit VectorXD(Size n) : data_(n) {}
    explicit VectorXD(const Eigen::VectorXd& v) : data_(v) {}

    Real& operator()(Size i) { return data_(i); }
    Real operator()(Size i) const { return data_(i); }

    Size size() const { return static_cast<Size>(data_.size()); }

    Real* data() { return data_.data(); }
    const Real* data() const { return data_.data(); }

    const Eigen::VectorXd& eigen() const { return data_; }
    Eigen::VectorXd& eigen() { return data_; }

    Real dot(const VectorXD& rhs) const {
        return data_.dot(rhs.data_);
    }

    VectorXD operator+(const VectorXD& rhs) const {
        return VectorXD(data_ + rhs.data_);
    }
    VectorXD operator-(const VectorXD& rhs) const {
        return VectorXD(data_ - rhs.data_);
    }

private:
    Eigen::VectorXd data_;
};

// =============================================================================
// 计量核心操作 (委托 Eigen + LAPACK)
// =============================================================================

// SVD 最小二乘求解 Ax = b (超定/适定系统)
// 复杂度: O(min(mn², m²n))
inline VectorXD svd_solve(const MatrixXD& A, const VectorXD& b) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A.eigen(), Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd x = svd.solve(b.eigen());
    return VectorXD(x);
}

// 对称正定矩阵求逆 (X'X)^{-1} via LLT
// 注: LLT 仅适用于 SPD 矩阵, 复杂度 n³/3, 比 LU 快约 2x
//     非 SPD 矩阵应使用 solve() (LU) 或 svd_solve() (SVD)
// 排幻觉点: LLT 适用范围窄于 LU, 非 SPD 矩阵会抛异常
inline MatrixXD inverse_symmetric(const MatrixXD& A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("inverse_symmetric: matrix must be square");
    }
    Eigen::LLT<Eigen::MatrixXd> llt(A.eigen());
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("inverse_symmetric: matrix is not positive definite (LLT failed)");
    }
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(), A.cols());
    Eigen::MatrixXd inv = llt.solve(I);
    return MatrixXD(inv);
}

// Cholesky 分解 (下三角 L, A = L * L')
inline MatrixXD cholesky_dynamic(const MatrixXD& A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("cholesky_dynamic: matrix must be square");
    }
    Eigen::LLT<Eigen::MatrixXd> llt(A.eigen());
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("cholesky_dynamic: matrix is not positive definite");
    }
    // LLT 存储 L, matrixL() 返回下三角视图
    Eigen::MatrixXd L = llt.matrixL();
    return MatrixXD(L);
}

// Ax = b 求解 via LU (通用, 适用于非 SPD 矩阵)
// 复杂度: 2n³/3
// 注: PartialPivLU 不报告 info() (与 LLT/SVD 不同), 奇异矩阵会返回 inf/nan
//     如需检测奇异性, 用 svd_solve() 检查奇异值
inline VectorXD solve(const MatrixXD& A, const VectorXD& b) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("solve: matrix must be square for LU solve");
    }
    if (A.rows() != b.size()) {
        throw std::invalid_argument("solve: dimension mismatch between A and b");
    }
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(A.eigen());
    Eigen::VectorXd x = lu.solve(b.eigen());
    return VectorXD(x);
}

// QR 分解, 返回 Q (正交矩阵)
inline MatrixXD qr_decompose(const MatrixXD& A) {
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A.eigen());
    Eigen::MatrixXd Q = qr.householderQ();
    return MatrixXD(Q);
}

// 完整 SVD 分解 A = U * diag(S) * V'
inline void svd_full(const MatrixXD& A, MatrixXD& U, VectorXD& S, MatrixXD& V) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A.eigen(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = MatrixXD(svd.matrixU());
    S = VectorXD(svd.singularValues());
    V = MatrixXD(svd.matrixV());
}

}  // namespace linalg::dynamic
}  // namespace v1
}  // namespace cpphub
