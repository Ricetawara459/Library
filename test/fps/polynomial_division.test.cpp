#define PROBLEM "https://judge.yosupo.jp/problem/division_of_polynomials"
#include <iostream>
#include "../../src/fps/fps.hpp"
int main() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);
    int n, m;
    if (!(std::cin >> n >> m)) return 0;
    fps<modint998244353> a(n), b(m);
    for (auto& x : a) std::cin >> x;
    for (auto& x : b) std::cin >> x;
    auto [q, r] = a.divmod(b);
    std::cout << q.size() << ' ' << r.size() << '\n';
    for (int i = 0; i < int(q.size()); i++) std::cout << q[i] << (i + 1 == int(q.size()) ? '\n' : ' ');
    if (q.empty()) std::cout << '\n';
    for (int i = 0; i < int(r.size()); i++) std::cout << r[i] << (i + 1 == int(r.size()) ? '\n' : ' ');
    if (r.empty()) std::cout << '\n';
}
