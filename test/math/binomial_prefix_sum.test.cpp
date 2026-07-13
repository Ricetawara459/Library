#include <cassert>
#include <iostream>
#include <vector>

#include "../../src/math/modint.hpp"
#include "../../src/math/factorial.hpp"
#include "../../src/math/binomial_prefix_sum.hpp"

using mint2 = modint998244353;

mint2 brute(factorial<mint2>& fac, int n, int r) {
    mint2 res = 0;
    for (int i = 0; i < r; i++) res += fac.comb(n, i);
    return res;
}

int main() {
    factorial<mint2> fac(100);
    binomial_prefix_sum<mint2> bp(5, 3);
    assert(bp.get() == mint2(16));
    assert(bp.comb() == fac.comb(5, 3));

    bp.inc_r();
    assert(bp.n() == 5 && bp.r() == 4);
    assert(bp.get() == brute(fac, 5, 4));
    assert(bp.comb() == fac.comb(5, 4));

    bp.dec_r();
    assert(bp.get() == brute(fac, 5, 3));
    assert(bp.comb() == fac.comb(5, 3));

    bp.inc_n();
    assert(bp.n() == 6 && bp.r() == 3);
    assert(bp.get() == brute(fac, 6, 3));
    assert(bp.comb() == fac.comb(6, 3));

    bp.dec_n();
    assert(bp.n() == 5 && bp.r() == 3);
    assert(bp.get() == brute(fac, 5, 3));
    assert(bp.comb() == fac.comb(5, 3));

    std::vector<std::pair<int, int>> path = {
        {6, 3}, {7, 3}, {7, 4}, {8, 4}, {8, 5},
        {9, 5}, {9, 4}, {10, 4}, {10, 6}, {12, 6},
        {11, 6}, {11, 2}, {3, 2}, {3, 0},
    };
    for (auto [n, r] : path) {
        bp.move(n, r);
        assert(bp.n() == n && bp.r() == r);
        assert(bp.get() == brute(fac, n, r));
        assert(bp.comb() == fac.comb(n, r));
    }

    int n, r;
    if (std::cin >> n >> r) {
        binomial_prefix_sum<mint2> solver(n, r);
        std::cout << solver.get().val() << '\n';
    }
}
