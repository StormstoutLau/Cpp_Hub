// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.1 - econometrics/core 基础类型 (任务 1.2-1.5, ADR-002/ADR-013)
// 内容: CovarianceType/数据容器/估计结果/Estimator 抽象基类
#include <gtest/gtest.h>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "cpphub/econometrics/core/covariance_type.hpp"
#include "cpphub/econometrics/core/data_types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"
#include "cpphub/econometrics/core/estimator_base.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::linalg::dynamic::MatrixXD;
using cpphub::v1::linalg::dynamic::VectorXD;

namespace {

// 用于 Estimator 抽象基类测试的 Mock: 仅实现纯虚 (estimate/name/clone)
// 其余方法使用基类默认行为 (isParametric=true, estimatorClass=Parametric)
class MockEstimator final : public Estimator {
public:
    explicit MockEstimator(CovarianceType ct = CovarianceType::Classical) { cov_type_ = ct; }

    EstimationResult estimate(const EconData&) override {
        EstimationResult r;
        r.cov_type = covarianceType();
        return r;
    }

    std::string name() const override { return "mock_estimator"; }

    std::unique_ptr<Estimator> clone() const override {
        return std::make_unique<MockEstimator>(covarianceType());
    }
};

}  // namespace

// =============================================================================
// §1 CovarianceType (任务 1.2)
// =============================================================================

TEST(EconometricsCoreTest, CovarianceType_ToString) {
    EXPECT_EQ(to_string(CovarianceType::Classical), "Classical");
    EXPECT_EQ(to_string(CovarianceType::HC0), "HC0");
    EXPECT_EQ(to_string(CovarianceType::HC1), "HC1");
    EXPECT_EQ(to_string(CovarianceType::HC2), "HC2");
    EXPECT_EQ(to_string(CovarianceType::HC3), "HC3");
    EXPECT_EQ(to_string(CovarianceType::HC4), "HC4");
    EXPECT_EQ(to_string(CovarianceType::HC5), "HC5");
    EXPECT_EQ(to_string(CovarianceType::HAC_Bartlett), "HAC_Bartlett");
    EXPECT_EQ(to_string(CovarianceType::HAC_QuadraticSpectral), "HAC_QuadraticSpectral");
    EXPECT_EQ(to_string(CovarianceType::HAC_Parzen), "HAC_Parzen");
    EXPECT_EQ(to_string(CovarianceType::HAC_TukeyHanning), "HAC_TukeyHanning");
    EXPECT_EQ(to_string(CovarianceType::Cluster_OneWay), "Cluster_OneWay");
    EXPECT_EQ(to_string(CovarianceType::Cluster_TwoWay), "Cluster_TwoWay");
    EXPECT_EQ(to_string(CovarianceType::OPG), "OPG");
    EXPECT_EQ(to_string(CovarianceType::Hessian), "Hessian");
    EXPECT_EQ(to_string(CovarianceType::Sandwich), "Sandwich");
    EXPECT_EQ(to_string(CovarianceType::Bootstrap), "Bootstrap");
    EXPECT_EQ(to_string(CovarianceType::Custom), "Custom");
}

// =============================================================================
// §2 数据容器 (任务 1.3)
// =============================================================================

TEST(EconometricsCoreTest, CrossSectionData_Construct) {
    CrossSectionData data;
    data.X = MatrixXD(4, 2);
    data.y = VectorXD(4);
    data.X(0, 0) = 1.0;
    data.y(0) = 2.5;
    data.x_names = {"const", "income"};
    data.y_name = "consumption";

    EXPECT_EQ(data.X.rows(), 4);
    EXPECT_EQ(data.X.cols(), 2);
    EXPECT_EQ(data.y.size(), 4);
    EXPECT_DOUBLE_EQ(data.X(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(data.y(0), 2.5);
    ASSERT_EQ(data.x_names.size(), 2);
    EXPECT_EQ(data.x_names[1], "income");
    EXPECT_EQ(data.y_name, "consumption");
}

TEST(EconometricsCoreTest, PanelData_Construct) {
    PanelData data;
    data.X = MatrixXD(6, 2);
    data.y = VectorXD(6);
    data.entity_id = {1, 1, 1, 2, 2, 2};
    data.time_id = {1, 2, 3, 1, 2, 3};
    data.balanced = true;

    EXPECT_EQ(data.entity_id.size(), 6);
    EXPECT_EQ(data.time_id.size(), 6);
    EXPECT_EQ(data.entity_id[3], 2);
    EXPECT_EQ(data.time_id[2], 3);
    EXPECT_TRUE(data.balanced);
    EXPECT_EQ(data.X.rows(), 6);
    EXPECT_EQ(data.y.size(), 6);
}

TEST(EconometricsCoreTest, TimeSeriesData_Construct) {
    TimeSeriesData data;
    data.y = VectorXD(5);
    data.X = MatrixXD(5, 1);
    data.timestamps = {1.0, 2.0, 3.0, 4.0, 5.0};
    data.y(1) = 10.0;

    EXPECT_EQ(data.timestamps.size(), 5);
    EXPECT_DOUBLE_EQ(data.timestamps[0], 1.0);
    EXPECT_DOUBLE_EQ(data.y(1), 10.0);
    EXPECT_EQ(data.X.rows(), 5);
}

TEST(EconometricsCoreTest, MakeCrossSection_Helper) {
    MatrixXD X(3, 2);
    VectorXD y(3);
    auto data = make_cross_section(X, y, {"a", "b"}, "target");

    EXPECT_EQ(data.X.rows(), 3);
    EXPECT_EQ(data.y.size(), 3);
    ASSERT_EQ(data.x_names.size(), 2);
    EXPECT_EQ(data.x_names[0], "a");
    EXPECT_EQ(data.y_name, "target");
}

TEST(EconometricsCoreTest, MakePanel_Helper) {
    MatrixXD X(4, 2);
    VectorXD y(4);
    auto data = make_panel(X, y, {1, 1, 2, 2}, {1, 2, 1, 2}, {"k1", "k2"}, "y", true);

    EXPECT_EQ(data.entity_id[0], 1);
    EXPECT_EQ(data.time_id[3], 2);
    EXPECT_TRUE(data.balanced);
    EXPECT_EQ(data.x_names[1], "k2");
}

TEST(EconometricsCoreTest, MakeTimeSeries_Helper) {
    MatrixXD X(3, 1);
    VectorXD y(3);
    auto data = make_time_series(y, X, {0.0, 1.0, 2.0}, {"f"}, "z");

    EXPECT_EQ(data.timestamps.size(), 3);
    EXPECT_EQ(data.y.size(), 3);
    EXPECT_EQ(data.x_names[0], "f");
    EXPECT_EQ(data.y_name, "z");
}

TEST(EconometricsCoreTest, EconData_Variant_Access) {
    CrossSectionData cs;
    cs.X = MatrixXD(2, 1);
    cs.y = VectorXD(2);

    PanelData pd;
    pd.X = MatrixXD(2, 1);
    pd.y = VectorXD(2);

    EconData e1 = cs;
    EconData e2 = pd;

    EXPECT_TRUE(std::holds_alternative<CrossSectionData>(e1));
    EXPECT_TRUE(std::holds_alternative<PanelData>(e2));
    EXPECT_EQ(std::get<CrossSectionData>(e1).X.rows(), 2);
    EXPECT_EQ(std::get<PanelData>(e2).X.rows(), 2);
    EXPECT_EQ(std::get<PanelData>(e2).entity_id.size(), 0);
}

TEST(EconometricsCoreTest, MakeVariantHelpers) {
    MatrixXD X(2, 1);
    VectorXD y(2);
    EconData e = make_cross_section(X, y, {"x"}, "y");
    EXPECT_TRUE(std::holds_alternative<CrossSectionData>(e));
}

// =============================================================================
// §3 估计结果 (任务 1.4)
// =============================================================================

TEST(EconometricsCoreTest, EstimationResult_Construct) {
    EstimationResult r;
    r.coefficients = VectorXD(2);
    r.coefficients(0) = 1.0;
    r.std_errors = VectorXD(2);
    r.t_statistics = VectorXD(2);
    r.p_values = VectorXD(2);
    r.vcov = MatrixXD(2, 2);
    r.log_likelihood = -10.0;
    r.r_squared = 0.9;
    r.adj_r_squared = 0.89;
    r.n_obs = 100;
    r.n_params = 2;
    r.df_residual = 98;
    r.cov_type = CovarianceType::HC1;

    EXPECT_EQ(r.coefficients.size(), 2);
    EXPECT_DOUBLE_EQ(r.coefficients(0), 1.0);
    EXPECT_DOUBLE_EQ(r.log_likelihood, -10.0);
    EXPECT_DOUBLE_EQ(r.r_squared, 0.9);
    EXPECT_DOUBLE_EQ(r.adj_r_squared, 0.89);
    EXPECT_EQ(r.n_obs, 100);
    EXPECT_EQ(r.n_params, 2);
    EXPECT_EQ(r.df_residual, 98);
    EXPECT_EQ(r.cov_type, CovarianceType::HC1);
    EXPECT_EQ(r.vcov.rows(), 2);
}

TEST(EconometricsCoreTest, InferenceResult_Construct) {
    InferenceResult r;
    r.coefficients = VectorXD(1);
    r.vcov = MatrixXD(1, 1);
    r.wald_statistic = 12.5;
    r.lr_statistic = 3.2;
    r.lm_statistic = 0.1;
    r.wald_pvalue = 0.02;
    r.lr_pvalue = 0.07;
    r.lm_pvalue = 0.75;
    r.df_test = 1;

    EXPECT_DOUBLE_EQ(r.wald_statistic, 12.5);
    EXPECT_DOUBLE_EQ(r.lr_statistic, 3.2);
    EXPECT_DOUBLE_EQ(r.lm_statistic, 0.1);
    EXPECT_DOUBLE_EQ(r.wald_pvalue, 0.02);
    EXPECT_DOUBLE_EQ(r.lr_pvalue, 0.07);
    EXPECT_DOUBLE_EQ(r.lm_pvalue, 0.75);
    EXPECT_EQ(r.df_test, 1);
}

TEST(EconometricsCoreTest, BootstrapResult_Construct) {
    BootstrapResult r;
    r.coef_mean = VectorXD(2);
    r.coef_std = VectorXD(2);
    r.coef_mean(0) = 3.5;
    r.coef_vcov = MatrixXD(2, 2);
    r.bootstrap_samples = {VectorXD(2), VectorXD(2), VectorXD(2)};
    r.lower_ci = -0.5;
    r.upper_ci = 7.5;
    r.n_replicates = 3;
    r.n_failed = 0;

    EXPECT_DOUBLE_EQ(r.coef_mean(0), 3.5);
    EXPECT_EQ(r.bootstrap_samples.size(), 3);
    EXPECT_DOUBLE_EQ(r.lower_ci, -0.5);
    EXPECT_DOUBLE_EQ(r.upper_ci, 7.5);
    EXPECT_EQ(r.n_replicates, 3);
    EXPECT_EQ(r.n_failed, 0);
    EXPECT_EQ(r.coef_vcov.rows(), 2);
}

// =============================================================================
// §4 Estimator 抽象基类 (任务 1.5, ADR-002)
// =============================================================================

TEST(EconometricsCoreTest, EstimatorClass_Enum) {
    EXPECT_EQ(static_cast<int>(EstimatorClass::Parametric), 0);
    EXPECT_EQ(static_cast<int>(EstimatorClass::Semiparametric), 1);
    EXPECT_EQ(static_cast<int>(EstimatorClass::Nonparametric), 2);
    EXPECT_EQ(static_cast<int>(EstimatorClass::MachineLearning), 3);
}

TEST(EconometricsCoreTest, Estimator_Derived_Mock) {
    MockEstimator est(CovarianceType::HC0);
    CrossSectionData data;
    data.X = MatrixXD(3, 1);
    data.y = VectorXD(3);

    EconData e = data;
    auto res = est.estimate(e);
    EXPECT_EQ(res.cov_type, CovarianceType::HC0);
    EXPECT_EQ(est.name(), "mock_estimator");
}

TEST(EconometricsCoreTest, Estimator_Clone) {
    MockEstimator est(CovarianceType::HC2);
    std::unique_ptr<Estimator> c = est.clone();

    EXPECT_NE(c.get(), nullptr);
    EXPECT_NE(c.get(), &est);
    EXPECT_EQ(c->name(), "mock_estimator");
    EXPECT_EQ(c->covarianceType(), CovarianceType::HC2);

    // 多态调用: clone 出的对象仍是 MockEstimator
    EXPECT_FALSE(c->isParametric() == false);
}

TEST(EconometricsCoreTest, Estimator_SetCovarianceType) {
    MockEstimator est;
    EXPECT_EQ(est.covarianceType(), CovarianceType::Classical);
    est.setCovarianceType(CovarianceType::HAC_Bartlett);
    EXPECT_EQ(est.covarianceType(), CovarianceType::HAC_Bartlett);
}

TEST(EconometricsCoreTest, Estimator_IsParametric_Default) {
    MockEstimator est;
    EXPECT_TRUE(est.isParametric());
    EXPECT_FALSE(est.isSemiparametric());
    EXPECT_FALSE(est.isNonparametric());
}

TEST(EconometricsCoreTest, Estimator_EstimatorClass_Default) {
    MockEstimator est;
    EXPECT_EQ(est.estimatorClass(), EstimatorClass::Parametric);
}