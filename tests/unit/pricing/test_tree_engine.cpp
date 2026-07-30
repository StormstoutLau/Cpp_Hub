#include <gtest/gtest.h>
#include "cpphub/pricing/tree/binomial.hpp"
#include "cpphub/pricing/tree/trinomial.hpp"
#include "cpphub/pricing/tree/tree_engine.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
#include "cpphub/core/math.hpp"
#include <cmath>

using namespace cpphub;

inline double bsm_call(double S0, double K, double T, double r, double q, double sigma) {
    double d1 = (std::log(S0/K) + (r - q + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
    double d2 = d1 - sigma*std::sqrt(T);
    return S0*std::exp(-q*T)*normal_cdf(d1) - K*std::exp(-r*T)*normal_cdf(d2);
}

inline double bsm_put(double S0, double K, double T, double r, double q, double sigma) {
    double d1 = (std::log(S0/K) + (r - q + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
    double d2 = d1 - sigma*std::sqrt(T);
    return K*std::exp(-r*T)*normal_cdf(-d2) - S0*std::exp(-q*T)*normal_cdf(-d1);
}

TEST(tree_BinomialCRR, EuropeanCallMatchesBSM) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    BinomialTreeEngine bte(p, BinomialType::CRR);
    CallPayOff call(100);
    Real price = bte.price_european(call);
    Real bsm = bsm_call(100, 100, 1.0, 0.05, 0.0, 0.2);
    EXPECT_NEAR(price, bsm, bsm * 0.001);
}

TEST(tree_BinomialCRR, EuropeanPutMatchesBSM) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    BinomialTreeEngine bte(p, BinomialType::CRR);
    PutPayOff put(100);
    Real price = bte.price_european(put);
    Real bsm = bsm_put(100, 100, 1.0, 0.05, 0.0, 0.2);
    EXPECT_NEAR(price, bsm, bsm * 0.001);
}

TEST(tree_BinomialCRR, AmericanPutGreaterThanEuropean) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    BinomialTreeEngine bte(p, BinomialType::CRR);
    PutPayOff put(100);
    Real euro = bte.price_european(put);
    Real amer = bte.price_american(put);
    EXPECT_GE(amer + 1e-12, euro);
}

TEST(tree_BinomialCRR, AmericanPutMatchesBenchmark) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    BinomialTreeEngine bte(p, BinomialType::CRR);
    PutPayOff put(100);
    Real price = bte.price_american(put);
    EXPECT_NEAR(price, 6.0909, 6.0909 * 0.005);
}

TEST(tree_BinomialLR, EuropeanCallMatchesBSMHighAccuracy) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.02, 0.2, 201};
    BinomialTreeEngine bte(p, BinomialType::LeisenReimer);
    CallPayOff call(100);
    Real price = bte.price_european(call);
    Real bsm = bsm_call(100, 100, 1.0, 0.05, 0.02, 0.2);
    EXPECT_NEAR(price, bsm, bsm * 1e-5);
}

TEST(tree_BinomialLR, EuropeanPutMatchesBSMHighAccuracy) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.02, 0.2, 201};
    BinomialTreeEngine bte(p, BinomialType::LeisenReimer);
    PutPayOff put(100);
    Real price = bte.price_european(put);
    Real bsm = bsm_put(100, 100, 1.0, 0.05, 0.02, 0.2);
    EXPECT_NEAR(price, bsm, bsm * 1e-5);
}

TEST(tree_BinomialLR, ConvergesFasterThanCRR) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.02, 0.2, 201};
    CallPayOff call(100);
    Real bsm = bsm_call(100, 100, 1.0, 0.05, 0.02, 0.2);

    BinomialTreeEngine crr(p, BinomialType::CRR);
    Real crr_price = crr.price_european(call);

    BinomialTreeEngine lr(p, BinomialType::LeisenReimer);
    Real lr_price = lr.price_european(call);

    Real crr_err = std::abs(crr_price - bsm);
    Real lr_err = std::abs(lr_price - bsm);
    EXPECT_LT(lr_err, crr_err / 10.0);
}

TEST(tree_BinomialLR, AmericanPutHighAccuracy) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 501};
    BinomialTreeEngine bte(p, BinomialType::LeisenReimer);
    PutPayOff put(100);
    Real price = bte.price_american(put);
    EXPECT_NEAR(price, 6.0909, 6.0909 * 0.001);
}

TEST(tree_BinomialTian, EuropeanCallMatchesBSM) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    BinomialTreeEngine bte(p, BinomialType::Tian);
    CallPayOff call(100);
    Real price = bte.price_european(call);
    Real bsm = bsm_call(100, 100, 1.0, 0.05, 0.0, 0.2);
    EXPECT_NEAR(price, bsm, bsm * 0.0005);
}

TEST(tree_BinomialTian, MoreStableThanCRRForDeepOTM) {
    BinomialParams p{100, 200, 1.0, 0.05, 0.0, 0.2, 1000};
    CallPayOff call(200);
    Real bsm = bsm_call(100, 200, 1.0, 0.05, 0.0, 0.2);

    BinomialTreeEngine crr(p, BinomialType::CRR);
    Real crr_price = crr.price_european(call);

    BinomialTreeEngine tian(p, BinomialType::Tian);
    Real tian_price = tian.price_european(call);

    Real crr_err = std::abs(crr_price - bsm);
    Real tian_err = std::abs(tian_price - bsm);
    EXPECT_LT(tian_err, crr_err);
}

TEST(tree_Trinomial, EuropeanCallMatchesBSM) {
    TrinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    TrinomialTreeEngine tte(p, TrinomialType::Explicit);
    CallPayOff call(100);
    Real price = tte.price_european(call);
    Real bsm = bsm_call(100, 100, 1.0, 0.05, 0.0, 0.2);
    EXPECT_NEAR(price, bsm, bsm * 0.001);
}

TEST(tree_Trinomial, AmericanPutGreaterThanEuropean) {
    TrinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    TrinomialTreeEngine tte(p, TrinomialType::Explicit);
    PutPayOff put(100);
    Real euro = tte.price_european(put);
    Real amer = tte.price_american(put);
    EXPECT_GE(amer + 1e-12, euro);
}

TEST(tree_Trinomial, ExplicitMatchesBinomialAtLargeN) {
    BinomialParams bp{100, 100, 1.0, 0.05, 0.0, 0.2, 2000};
    TrinomialParams tp{100, 100, 1.0, 0.05, 0.0, 0.2, 2000};

    BinomialTreeEngine bte(bp, BinomialType::CRR);
    TrinomialTreeEngine tte(tp, TrinomialType::Explicit);
    CallPayOff call(100);

    Real binom_price = bte.price_european(call);
    Real trinom_price = tte.price_european(call);
    Real ref = std::max(std::abs(binom_price), std::abs(trinom_price));
    EXPECT_NEAR(binom_price, trinom_price, ref * 0.0005);
}

TEST(tree_TreeEngine, BinomialVariantWorks) {
    BinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    BinomialTreeEngine bte(p, BinomialType::CRR);
    TreeEngine te(bte);
    CallPayOff call(100);
    Real price = te.price(call, OptionType::European);
    Real bsm = bsm_call(100, 100, 1.0, 0.05, 0.0, 0.2);
    EXPECT_GT(price, 0.0);
    EXPECT_NEAR(price, bsm, bsm * 0.001);

    Real delta = te.delta(call, OptionType::European);
    Real gamma = te.gamma(call, OptionType::European);
    Real theta = te.theta(call, OptionType::European);
    Real d1 = (std::log(100.0/100.0) + (0.05 + 0.5*0.04)*1.0) / (0.2*1.0);
    Real bsm_delta = std::exp(-0.0*1.0) * normal_cdf(d1);
    EXPECT_NEAR(delta, bsm_delta, 0.15);
    EXPECT_GT(gamma, 0.0);
    EXPECT_TRUE(theta <= 0.0 || std::abs(theta) < 0.1);
    EXPECT_EQ(te.name(), "BinomialTreeEngine");
}

TEST(tree_TreeEngine, TrinomialVariantWorks) {
    TrinomialParams p{100, 100, 1.0, 0.05, 0.0, 0.2, 1000};
    TrinomialTreeEngine tte(p, TrinomialType::Explicit);
    TreeEngine te(tte);
    PutPayOff put(100);
    Real price = te.price(put, OptionType::European);
    Real bsm = bsm_put(100, 100, 1.0, 0.05, 0.0, 0.2);
    EXPECT_GT(price, 0.0);
    EXPECT_NEAR(price, bsm, bsm * 0.001);

    Real delta = te.delta(put, OptionType::European);
    Real gamma = te.gamma(put, OptionType::European);
    Real d1 = (std::log(100.0/100.0) + (0.05 + 0.5*0.04)*1.0) / (0.2*1.0);
    Real bsm_delta = -std::exp(-0.0*1.0) * normal_cdf(-d1);
    EXPECT_NEAR(delta, bsm_delta, 0.15);
    EXPECT_GT(gamma, 0.0);
    EXPECT_EQ(te.name(), "TrinomialTreeEngine");
}
