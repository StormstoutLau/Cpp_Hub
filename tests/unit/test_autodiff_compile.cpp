// autodiff 头文件 MSVC 编译测试 (不依赖 Eigen)
// 注意: forward/dual 和 reverse/var 不能在同一编译单元同时 using namespace autodiff (wrt 歧义)
// 实际使用时: ad_dual.hpp 只包含 forward, aad_greeks.hpp 只包含 reverse

#include <iostream>

// === 前向模式测试 ===
#if 0
#include <autodiff/forward/dual.hpp>
using namespace autodiff;

dual forward_test(dual x) {
    return x * x + 2.0 * x + 1.0;
}

void test_forward() {
    dual x = 1.0;
    dual y = forward_test(x);
    double dydx = derivative(forward_test, wrt(x), at(x));
    std::cout << "Forward: y=" << val(y) << " dy/dx=" << dydx << std::endl;
}
#endif

// === 反向模式测试 ===
#include <autodiff/reverse/var.hpp>
using namespace autodiff;

var reverse_test(var a, var b) {
    return a * a + b * b;
}

void test_reverse() {
    var a = 3.0;
    var b = 4.0;
    var c = reverse_test(a, b);
    auto [ca, cb] = derivatives(c, wrt(a, b));
    std::cout << "Reverse: c=" << val(c) << " dc/da=" << ca << " dc/db=" << cb << std::endl;
}

int main() {
    test_reverse();
    return 0;
}
