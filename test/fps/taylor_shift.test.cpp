#define PROBLEM "https://judge.yosupo.jp/problem/polynomial_taylor_shift"

#include <iostream>
#include <vector>
#include "../../src/fps/fps.hpp"

int main() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);
    int n;
    modint998244353 c;
    if (!(std::cin >> n >> c)) return 0;
    fps<modint998244353> f(n);
    for (auto& x : f) std::cin >> x;
    auto result = f.taylor_shift(c);
    for (int i = 0; i < n; i++) std::cout << result[i] << (i + 1 == n ? '\n' : ' ');
}
