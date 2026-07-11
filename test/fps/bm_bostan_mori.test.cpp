#include <cassert>
#include <vector>

#include "../../src/fps/fps.hpp"
#include "../../src/math/modint.hpp"

using mint = modint998244353;

int main() {
    std::vector<mint> fib{0, 1, 1, 2, 3, 5, 8, 13, 21, 34};
    auto coef = berlekamp_massey(fib);
    assert((coef == std::vector<mint>{1, 1}));
    fps<mint> ba = fps<mint>(berlekamp_massey(fib));
    assert((ba == fps<mint>{1, 1}));
    assert(linear_recurrence_kth(fib, coef, 0) == mint(0));
    assert(linear_recurrence_kth(fib, coef, 1) == mint(1));
    assert(linear_recurrence_kth(fib, coef, 10) == mint(55));
    assert(bmbm(fib, 10) == mint(55));

    formal_power_series<mint> p{1};
    formal_power_series<mint> q{1, -1, -1};
    assert(bostan_mori(p, q, 0) == mint(1));
    assert(bostan_mori(p, q, 1) == mint(1));
    assert(bostan_mori(p, q, 9) == mint(55));

    std::vector<mint> seq;
    for (int i = 0; i < 20; i++) seq.push_back(mint(3).pow(i));
    assert((berlekamp_massey(seq) == std::vector<mint>{3}));
    assert(bmbm(seq, 50) == mint(3).pow(50));
}
