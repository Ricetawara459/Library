#define PROBLEM "https://judge.yosupo.jp/problem/multipoint_evaluation"
#include <iostream>
#include <vector>
#include "../../src/fps/fps.hpp"
int main() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);
    int n, m;
    if (!(std::cin >> n >> m)) return 0;
    fps<modint998244353> f(n);
    std::vector<modint998244353> xs(m);
    for (auto& x : f) std::cin >> x;
    for (auto& x : xs) std::cin >> x;
    auto ans = f.multipoint_eval(xs);
    for (int i = 0; i < m; i++) std::cout << ans[i] << (i + 1 == m ? '\n' : ' ');
}
