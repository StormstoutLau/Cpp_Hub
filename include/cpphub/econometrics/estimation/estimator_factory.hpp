// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.3 任务 2.5 - EstimatorFactory (ADR-003 风格)
// 估计器工厂 + 静态注册模板, 运行时按名称/枚举创建估计器实例
//
// 教材锚点:
//   - ADR-003 (Factory + 静态注册模板, Joshi Ch.10)
//   - ADR-002 (Bridge + Virtual Constructor, clone() 语义)
//   - Gamma et al. "Design Patterns" Abstract Factory
//
// 设计要点:
//   1. 单例工厂 (Meyers Singleton, C++11 线程安全局部变量初始化)
//   2. Creator 函数对象: std::function<std::unique_ptr<Estimator>()>
//   3. EstimatorRegistrar<T> 静态注册模板 (链接即注册, 无需修改工厂代码)
//   4. 同时支持按字符串名称和按 MLEFamily 枚举创建 (MLEFamily 自动转为字符串)
//
// 排幻觉点 F1 (clone 语义 vs 共享状态):
//   工厂每次 create() 都返回全新的 unique_ptr 实例, 不缓存也不共享内部状态.
//   原因: 估计器持有迭代历史 (n_iter_/converged_)/协方差类型/容差等可变配置,
//   共享实例会导致多调用方相互污染. ADR-002 的 clone() 保证派生类完整复制.
//
// 排幻觉点 F2 (未注册名称抛异常, 不返回 nullptr):
//   未注册的 name 抛 std::invalid_argument, 不返回 nullptr.
//   原因: 返回 nullptr 会将错误检测推迟到调用方解引用时 (UB 风险),
//   违反 fail-fast 原则. RAII 要求构造即有效.
//
// 排幻觉点 F3 (工厂不持有引用, unique_ptr 转移所有权):
//   create() 返回 std::unique_ptr<Estimator>, 所有权转移给调用方.
//   工厂仅保存 Creator 函数 (无状态 lambda), 不持有任何 Estimator 实例.
//   原因: 工厂持有实例会导致生命周期管理复杂化 (何时销毁? 多调用方共享?),
//   违反"单一所有权"原则 (ADR-002 unique_ptr 单一所有权约定).
//
// 排幻觉点 F4 (静态注册顺序, Meyers Singleton 规避):
//   C++11 保证函数内静态变量初始化线程安全且仅一次 (Meyers Singleton).
//   传统全局静态变量初始化顺序未定义 (SIOF, Static Initialization Order Fiasco),
//   Meyers Singleton 规避 SIOF: 首次调用 instance() 时才初始化.
//   注册器 (EstimatorRegistrar<T>) 在文件作用域静态实例化,
//   其构造函数调用 instance() 触发工厂首次初始化, 顺序正确.
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <functional>   // std::function
#include <memory>       // std::unique_ptr, std::make_unique
#include <stdexcept>    // std::invalid_argument, std::runtime_error
#include <string>       // std::string
#include <unordered_map>  // std::unordered_map
#include <utility>      // std::move
#include <vector>       // std::vector

#include "cpphub/core/types.hpp"  // Size
#include "cpphub/econometrics/core/estimator_base.hpp"
#include "cpphub/econometrics/estimation/mle.hpp"
#include "cpphub/econometrics/estimation/ols.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// =============================================================================
// EstimatorFactory - 估计器工厂单例 (ADR-003)
//
// 用法:
//   // 1. 按名称创建
//   auto ols = EstimatorFactory::instance().create("OLS");
//   auto mle_logit = EstimatorFactory::instance().create("MLE.Logistic");
//
//   // 2. 按 MLEFamily 枚举创建 (便捷接口)
//   auto mle_pois = EstimatorFactory::instance().createMLE(MLEFamily::Poisson);
//
//   // 3. 自定义配置后估计
//   mle_logit->setCovarianceType(CovarianceType::Hessian);
//   EstimationResult r = mle_logit->estimate(data);
//
// 注册新估计器:
//   // 在 .hpp 或 .cpp 文件中
//   static EstimatorRegistrar<MyEstimator> reg_my("MyEstimator");
//   // 链接即注册, 无需修改 EstimatorFactory 代码
// =============================================================================
class EstimatorFactory {
public:
    /// Creator 函数类型: 无参构造 Estimator 派生类
    using Creator = std::function<std::unique_ptr<Estimator>()>;

private:
    std::unordered_map<std::string, Creator> registry_;

    // 私有构造 (单例)
    EstimatorFactory() = default;

public:
    /// @brief 工厂单例访问 (Meyers Singleton, C++11 线程安全)
    /// @return 工厂引用
    static EstimatorFactory& instance() {
        static EstimatorFactory factory;
        return factory;
    }

    // 禁用拷贝/移动 (单例语义)
    EstimatorFactory(const EstimatorFactory&) = delete;
    EstimatorFactory& operator=(const EstimatorFactory&) = delete;
    EstimatorFactory(EstimatorFactory&&) = delete;
    EstimatorFactory& operator=(EstimatorFactory&&) = delete;

    /// @brief 注册估计器 (Creator 函数)
    /// @param name 估计器名称 (如 "OLS", "MLE.Logistic")
    /// @param creator 创建函数 (返回 unique_ptr<Estimator>)
    /// @note 线程安全: 仅在静态初始化期调用 (注册器构造时)
    void registerEstimator(const std::string& name, Creator creator) {
        registry_[name] = std::move(creator);
    }

    /// @brief 按名称创建估计器实例
    /// @param name 估计器名称 (区分大小写)
    /// @return 新的估计器实例 (unique_ptr, 所有权转移给调用方)
    /// @throws std::invalid_argument 若 name 未注册 (排幻觉点 F2)
    std::unique_ptr<Estimator> create(const std::string& name) const {
        const auto it = registry_.find(name);
        if (it == registry_.end()) {
            throw std::invalid_argument(
                "EstimatorFactory::create: unknown estimator name '" + name +
                "' (not registered)");
        }
        // 排幻觉点 F1: 每次 create 返回全新实例, 不缓存
        // 排幻觉点 F3: unique_ptr 转移所有权, 工厂不持有
        return it->second();
    }

    /// @brief 按 MLEFamily 枚举创建 MLE 估计器 (便捷接口)
    /// @param family MLE 分布族
    /// @return 新的 MLEEstimator 实例 (unique_ptr, 所有权转移给调用方)
    /// @throws std::invalid_argument 若 family 未注册
    std::unique_ptr<MLEEstimator> createMLE(MLEFamily family) const {
        const std::string name = "MLE." + to_string(family);
        auto est = create(name);
        // 类型转换: 工厂内部以 Estimator 基类存储, 此处安全向下转换
        // (注册时确保 Creator 返回 MLEEstimator)
        MLEEstimator* mle_ptr = dynamic_cast<MLEEstimator*>(est.get());
        if (mle_ptr == nullptr) {
            throw std::runtime_error(
                "EstimatorFactory::createMLE: internal error - creator for '" +
                name + "' did not return MLEEstimator");
        }
        est.release();  // 释放所有权, 由 MLEEstimator* 接管
        return std::unique_ptr<MLEEstimator>(mle_ptr);
    }

    /// @brief 列出已注册的估计器名称
    /// @return 名称向量
    std::vector<std::string> listNames() const {
        std::vector<std::string> names;
        names.reserve(registry_.size());
        for (const auto& kv : registry_) {
            names.push_back(kv.first);
        }
        return names;
    }

    /// @brief 检查名称是否已注册
    /// @param name 估计器名称
    /// @return true 若已注册
    bool isRegistered(const std::string& name) const {
        return registry_.find(name) != registry_.end();
    }

    /// @brief 已注册估计器数量
    /// @return 数量
    Size size() const { return static_cast<Size>(registry_.size()); }
};

// =============================================================================
// EstimatorRegistrar - 静态注册模板 (ADR-003)
//
// 用法:
//   static EstimatorRegistrar<OLSEstimator> reg_ols("OLS");
//   static EstimatorRegistrar<MLEEstimator> reg_mle_gauss("MLE.Gaussian");
//
// 注: 对于需要特定构造参数的估计器 (如 MLEEstimator 需要 MLEFamily),
//     可特化 Creator lambda 而非直接用 EstimatorRegistrar<T>.
//     本文件末尾提供 register_mle_family() 辅助函数.
// =============================================================================
template <typename EstT>
class EstimatorRegistrar {
public:
    /// @brief 注册无参构造的估计器
    /// @param name 估计器名称
    explicit EstimatorRegistrar(const std::string& name) {
        EstimatorFactory::instance().registerEstimator(
            name, []() -> std::unique_ptr<Estimator> {
                return std::make_unique<EstT>();
            });
    }

    /// @brief 注册带构造参数的估计器 (Creator 自定义)
    /// @param name 估计器名称
    /// @param creator 创建函数
    EstimatorRegistrar(const std::string& name,
                        std::function<std::unique_ptr<Estimator>()> creator) {
        EstimatorFactory::instance().registerEstimator(name, std::move(creator));
    }
};

// =============================================================================
// 内置估计器注册 (链接即注册)
//
// 排幻觉点 F4: 静态注册器在文件作用域实例化, 其构造函数调用
//   EstimatorFactory::instance() 触发 Meyers Singleton 首次初始化,
//   规避 SIOF (Static Initialization Order Fiasco).
//
// 注册清单:
//   - "OLS"               -> OLSEstimator (默认 Classical 协方差)
//   - "MLE.Gaussian"      -> MLEEstimator(Gaussian)
//   - "MLE.Logistic"      -> MLEEstimator(Logistic)
//   - "MLE.Bernoulli"     -> MLEEstimator(Bernoulli)  [同 Logistic, 别名]
//   - "MLE.Probit"        -> MLEEstimator(Probit)
//   - "MLE.Poisson"       -> MLEEstimator(Poisson)
//   - "MLE.NegativeBinomial" -> MLEEstimator(NegativeBinomial)
// =============================================================================
namespace detail {

// OLS 注册: 默认 Classical 协方差, 调用方可 setCovarianceType 切换 HC/HAC
inline EstimatorRegistrar<OLSEstimator> register_ols("OLS");

// MLE 族注册: 每种 MLEFamily 单独注册, Creator 捕获 family 枚举
// 注: lambda 返回 unique_ptr<Estimator> 但实际为 MLEEstimator,
//     createMLE() 中 dynamic_cast 安全向下转换
inline EstimatorRegistrar<MLEEstimator> register_mle_gaussian(
    "MLE.Gaussian",
    []() -> std::unique_ptr<Estimator> {
        return std::make_unique<MLEEstimator>(MLEFamily::Gaussian);
    });

inline EstimatorRegistrar<MLEEstimator> register_mle_logistic(
    "MLE.Logistic",
    []() -> std::unique_ptr<Estimator> {
        return std::make_unique<MLEEstimator>(MLEFamily::Logistic);
    });

inline EstimatorRegistrar<MLEEstimator> register_mle_bernoulli(
    "MLE.Bernoulli",
    []() -> std::unique_ptr<Estimator> {
        return std::make_unique<MLEEstimator>(MLEFamily::Bernoulli);
    });

inline EstimatorRegistrar<MLEEstimator> register_mle_probit(
    "MLE.Probit",
    []() -> std::unique_ptr<Estimator> {
        return std::make_unique<MLEEstimator>(MLEFamily::Probit);
    });

inline EstimatorRegistrar<MLEEstimator> register_mle_poisson(
    "MLE.Poisson",
    []() -> std::unique_ptr<Estimator> {
        return std::make_unique<MLEEstimator>(MLEFamily::Poisson);
    });

inline EstimatorRegistrar<MLEEstimator> register_mle_nb(
    "MLE.NegativeBinomial",
    []() -> std::unique_ptr<Estimator> {
        return std::make_unique<MLEEstimator>(MLEFamily::NegativeBinomial);
    });

}  // namespace detail

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
