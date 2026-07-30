#include <gtest/gtest.h>
#include "cpphub/instruments/payoff/payoff.hpp"
#include "cpphub/instruments/payoff/payoff_bridge.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"

using namespace cpphub;

TEST(PayOffBridge, DefaultConstruct) {
    PayOffBridge b;
    EXPECT_FALSE(static_cast<bool>(b));
}

TEST(PayOffBridge, ConstructFromPayOffRef) {
    CallPayOff call(100.0);
    PayOffBridge b(call);
    EXPECT_TRUE(static_cast<bool>(b));
    EXPECT_DOUBLE_EQ(b(105.0), 5.0);
}

TEST(PayOffBridge, ConstructFromUniquePtr) {
    auto ptr = std::make_unique<CallPayOff>(50.0);
    PayOffBridge b(std::move(ptr));
    EXPECT_DOUBLE_EQ(b(60.0), 10.0);
}

TEST(PayOffBridge, CopyConstruct) {
    CallPayOff call(100.0);
    PayOffBridge b1(call);
    PayOffBridge b2(b1);
    EXPECT_DOUBLE_EQ(b1(110.0), 10.0);
    EXPECT_DOUBLE_EQ(b2(110.0), 10.0);
}

TEST(PayOffBridge, MoveConstruct) {
    CallPayOff call(100.0);
    PayOffBridge b1(call);
    PayOffBridge b2(std::move(b1));
    EXPECT_FALSE(static_cast<bool>(b1));
    EXPECT_DOUBLE_EQ(b2(120.0), 20.0);
}

TEST(PayOffBridge, CopyAssign) {
    CallPayOff call(100.0);
    PayOffBridge b1(call);
    PayOffBridge b2;
    b2 = b1;
    EXPECT_DOUBLE_EQ(b2(110.0), 10.0);
}

TEST(PayOffBridge, MoveAssign) {
    CallPayOff call(100.0);
    PayOffBridge b1(call);
    PayOffBridge b2;
    b2 = std::move(b1);
    EXPECT_FALSE(static_cast<bool>(b1));
    EXPECT_DOUBLE_EQ(b2(110.0), 10.0);
}

TEST(PayOffBridge, SelfAssign) {
    CallPayOff call(100.0);
    PayOffBridge b(call);
    PayOffBridge& ref = b;
    b = ref;
    EXPECT_DOUBLE_EQ(b(110.0), 10.0);
}

TEST(PayOffBridge, PolymorphicCall) {
    PayOffBridge b1(CallPayOff(100.0));
    PayOffBridge b2(PutPayOff(100.0));
    EXPECT_DOUBLE_EQ(b1(110.0), 10.0);
    EXPECT_DOUBLE_EQ(b2(110.0), 0.0);
    EXPECT_DOUBLE_EQ(b2(90.0), 10.0);
}

TEST(PayOffBridge, GetAccessor) {
    CallPayOff call(100.0);
    PayOffBridge b(call);
    const PayOff& ref = b.get();
    EXPECT_EQ(ref.name(), "Call");
}

TEST(PayOffBridge, SwapViaMove) {
    CallPayOff call1(100.0);
    PutPayOff put(50.0);
    PayOffBridge b1(call1);
    PayOffBridge b2(put);
    PayOffBridge tmp(std::move(b1));
    b1 = std::move(b2);
    b2 = std::move(tmp);
    EXPECT_DOUBLE_EQ(b1(60.0), 0.0);
    EXPECT_DOUBLE_EQ(b2(60.0), 0.0);
    EXPECT_DOUBLE_EQ(b2(110.0), 10.0);
}

TEST(PayOffBridge, EmptyBridgeCallThrows) {
    PayOffBridge b;
    EXPECT_ANY_THROW(b(100.0));
}
