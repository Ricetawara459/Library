#define PROBLEM "https://judge.yosupo.jp/problem/inv_of_formal_power_series"

#include <cassert>
#include <iostream>
#include <vector>

#include "../../src/math/modint.hpp"
#include "../../src/fps/fps.hpp"

using namespace std;
using mint = modint998244353;
using poly = fps<mint>;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n;
    if (cin >> n) {
        poly a(n);
        for (int i = 0; i < n; i++) {
            long long x;
            cin >> x;
            a[i] = mint(x);
        }
        auto b = a.inv(n);
        for (int i = 0; i < n; i++) {
            cout << b[i].val() << (i + 1 == n ? '\n' : ' ');
        }
        return 0;
    }

    poly a{1, 2, 3};
    poly b{4, 5};
    assert((a * b == poly{4, 13, 22, 15}));

    poly f{1, 1, 1, 1, 1, 1};
    auto inv = f.inv(6);
    assert(((f * inv).pre(6) == poly{1, 0, 0, 0, 0, 0}));

    poly g{0, 1, 2, 3, 4, 5};
    assert(g.exp(6).log(6) == g);

    poly p{1, 1, 0, 0, 0, 0};
    assert((p.pow(3, 6) == poly{1, 3, 3, 1, 0, 0}));

    poly sq{1, 2, 1, 0, 0};
    assert((sq.sqrt(5) == poly{1, 1, 0, 0, 0}));

    assert(a.eval(2) == mint(17));
    return 0;
}
