#include <gtest/gtest.h>
#include "cpphub/instruments/payoff/payoff.hpp"
#include "cpphub/instruments/payoff/payoff_bridge.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
#include "cpphub/instruments/payoff/exotic.hpp"
#include "cpphub/instruments/payoff/factory.hpp"
#include "cpphub/core/error.hpp"
#include <map>
#include <string>

using namespace cpphub;

class PayOffFactoryTest : public ::testing::Test {
protected:
    void SetUp() override {
        auto& f = PayOffFactory::instance();
        f.register_payoff("Call", [](const std::map<std::string, double>& p) {
            double k = p.count("strike") ? p.at("strike") : 0.0;
            return std::make_unique<CallPayOff>(k);
        });
        f.register_payoff("Put", [](const std::map<std::string, double>& p) {
            double k = p.count("strike") ? p.at("strike") : 0.0;
            return std::make_unique<PutPayOff>(k);
        });
        f.register_payoff("DigitalCall", [](const std::map<std::string, double>& p) {
            double k = p.count("strike") ? p.at("strike") : 0.0;
            double pay = p.count("payment") ? p.at("payment") : 1.0;
            return std::make_unique<DigitalCallPayOff>(k, pay);
        });
        f.register_payoff("DoubleDigital", [](const std::map<std::string, double>& p) {
            double lo = p.count("lower") ? p.at("lower") : 0.0;
            double up = p.count("upper") ? p.at("upper") : 0.0;
            double pay = p.count("payment") ? p.at("payment") : 1.0;
            return std::make_unique<DoubleDigitalPayOff>(lo, up, pay);
        });
    }
};

TEST_F(PayOffFactoryTest, CreateCall) {
    auto p = PayOffFactory::instance().create("Call", {{"strike", 100.0}});
    EXPECT_NE(p, nullptr);
    EXPECT_EQ(p->name(), "Call");
    EXPECT_DOUBLE_EQ((*p)(110.0), 10.0);
    EXPECT_DOUBLE_EQ((*p)(90.0), 0.0);
}

TEST_F(PayOffFactoryTest, CreatePut) {
    auto p = PayOffFactory::instance().create("Put", {{"strike", 50.0}});
    EXPECT_DOUBLE_EQ((*p)(40.0), 10.0);
    EXPECT_DOUBLE_EQ((*p)(60.0), 0.0);
}

TEST_F(PayOffFactoryTest, CreateDigitalCall) {
    auto p = PayOffFactory::instance().create("DigitalCall",
        {{"strike", 100.0}, {"payment", 5.0}});
    EXPECT_DOUBLE_EQ((*p)(110.0), 5.0);
    EXPECT_DOUBLE_EQ((*p)(90.0), 0.0);
}

TEST_F(PayOffFactoryTest, CreateDoubleDigital) {
    auto p = PayOffFactory::instance().create("DoubleDigital",
        {{"lower", 90.0}, {"upper", 110.0}, {"payment", 2.0}});
    EXPECT_DOUBLE_EQ((*p)(100.0), 2.0);
    EXPECT_DOUBLE_EQ((*p)(85.0), 0.0);
    EXPECT_DOUBLE_EQ((*p)(115.0), 0.0);
    EXPECT_DOUBLE_EQ((*p)(90.0), 2.0);
    EXPECT_DOUBLE_EQ((*p)(110.0), 2.0);
}

TEST_F(PayOffFactoryTest, UnknownTypeThrows) {
    EXPECT_THROW(PayOffFactory::instance().create("NonExistent"),
                 CppHubException);
}

TEST_F(PayOffFactoryTest, IsRegistered) {
    EXPECT_TRUE(PayOffFactory::instance().is_registered("Call"));
    EXPECT_TRUE(PayOffFactory::instance().is_registered("Put"));
    EXPECT_FALSE(PayOffFactory::instance().is_registered("Foo"));
}

TEST_F(PayOffFactoryTest, List) {
    auto names = PayOffFactory::instance().list();
    EXPECT_GE(names.size(), 4u);
    bool has_call = false;
    for (const auto& n : names) if (n == "Call") has_call = true;
    EXPECT_TRUE(has_call);
}

TEST_F(PayOffFactoryTest, SingletonUnique) {
    auto& f1 = PayOffFactory::instance();
    auto& f2 = PayOffFactory::instance();
    EXPECT_EQ(&f1, &f2);
}

TEST_F(PayOffFactoryTest, BridgeFromFactory) {
    auto p = PayOffFactory::instance().create("Call", {{"strike", 100.0}});
    PayOffBridge b(std::move(p));
    EXPECT_DOUBLE_EQ(b(120.0), 20.0);
}
