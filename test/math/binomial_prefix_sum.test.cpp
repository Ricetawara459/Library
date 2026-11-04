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

template <class Mint>
void test_boundaries(int limit) {
    std::vector<std::vector<Mint>> choose(limit + 1), prefix(limit + 1);
    std::vector<std::pair<int, int>> states;
    for (int n = 0; n <= limit; n++) {
        choose[n].assign(n + 1, 1);
        for (int r = 1; r < n; r++) choose[n][r] = choose[n - 1][r - 1] + choose[n - 1][r];
        prefix[n].assign(n + 2, 0);
        for (int r = 0; r <= n; r++) prefix[n][r + 1] = prefix[n][r] + choose[n][r];
        for (int r = 0; r <= n + 1; r++) states.emplace_back(n, r);
    }
    auto check = [&](const binomial_prefix_sum<Mint>& state) {
        int n = state.n(), r = state.r();
        assert(state.get() == prefix[n][r]);
        assert(state.comb() == (r <= n ? choose[n][r] : Mint(0)));
        if (r == n + 1) assert(state.get() == Mint(2).pow(n));
    };

    for (auto [n, r] : states) {
        binomial_prefix_sum<Mint> initial(n, r);
        check(initial);
        initial.reset(0, 1);
        check(initial);
        initial.reset(n, r);
        check(initial);
        initial.extend(limit);
        check(initial);
        if (r <= n) {
            auto state = initial;
            state.inc_r();
            check(state);
            state.dec_r();
            check(state);
        }
        if (r > 0) {
            auto state = initial;
            state.dec_r();
            check(state);
            state.inc_r();
            check(state);
        }
        if (n < limit) {
            auto state = initial;
            state.inc_n();
            check(state);
            state.dec_n();
            check(state);
        }
        if (n > 0 && r <= n) {
            auto state = initial;
            state.dec_n();
            check(state);
            state.inc_n();
            check(state);
        }
    }

    // 全ての状態間を移動し、n+1 の境界をまたぐ move も照合する。
    for (auto [n, r] : states) {
        if (n > 12) continue;
        binomial_prefix_sum<Mint> state(n, r);
        for (auto [next_n, next_r] : states) {
            if (next_n > 12) continue;
            state.move(n, r);
            check(state);
            state.move(next_n, next_r);
            check(state);
        }
    }

    binomial_prefix_sum<Mint> full(0, 1);
    for (int n = 0; n < limit; n++) {
        full.inc_n();
        check(full);
        full.inc_r();
        check(full);
    }
    for (int n = limit; n > 0; n--) {
        full.dec_r();
        check(full);
        full.dec_n();
        check(full);
    }

    binomial_prefix_sum<Mint> prepared;
    prepared.extend(limit);
    check(prepared);
    for (auto [n, r] : states) {
        prepared.move(n, r);
        check(prepared);
    }
}

int main() {
    test_boundaries<modint998244353>(32);
    test_boundaries<modint1000000007>(32);
    test_boundaries<static_modint<17>>(16);
    test_boundaries<static_modint<3>>(2);
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
        {3, 4}, {4, 4}, {3, 4}, {0, 1}, {0, 0},
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
