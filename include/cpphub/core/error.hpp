// SOURCE: PHASE1_SPEC §2.4 - 错误码、异常层级
#pragma once
#include <stdexcept>
#include <string>

namespace cpphub {
inline namespace v1 {

enum class ErrorCode {
    Success = 0,
    InvalidArgument,
    OutOfRange,
    ConvergenceFailure,
    NotConverged,
    NumericalError,
    NotImplemented,
};

class CppHubException : public std::runtime_error {
public:
    explicit CppHubException(const std::string& msg, ErrorCode code = ErrorCode::NumericalError)
        : std::runtime_error(msg), code_(code) {}
    ErrorCode code() const noexcept { return code_; }
private:
    ErrorCode code_;
};

inline const char* error_message(ErrorCode code) noexcept {
    switch (code) {
        case ErrorCode::Success: return "success";
        case ErrorCode::InvalidArgument: return "invalid argument";
        case ErrorCode::OutOfRange: return "out of range";
        case ErrorCode::ConvergenceFailure: return "convergence failure";
        case ErrorCode::NotConverged: return "not converged";
        case ErrorCode::NumericalError: return "numerical error";
        case ErrorCode::NotImplemented: return "not implemented";
        default: return "unknown error";
    }
}

}  // namespace v1
}  // namespace cpphub
