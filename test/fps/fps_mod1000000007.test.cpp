#include <cassert>
#include <vector>

#include "../../src/fps/fps.hpp"
#include "../../src/math/modint.hpp"

using mint = modint1000000007;
using poly = fps<mint>;

int main() {
    std::vector<mint> a{1000000006, 1000000005, 123456789};
    std::vector<mint> b{1000000004, 987654321};
    auto c = convolution(a, b);
    std::vector<mint> expected(a.size() + b.size() - 1);
    for (int i = 0; i < int(a.size()); i++) {
        for (int j = 0; j < int(b.size()); j++) {
            expected[i + j] += a[i] * b[j];
        }
    }
    assert(c == expected);

    poly f{1, 2, 3, 4, 5, 6};
    auto inv = f.inv(6);
    assert(((f * inv).pre(6) == poly{1, 0, 0, 0, 0, 0}));

    std::vector<mint> fib{0, 1, 1, 2, 3, 5, 8, 13, 21, 34};
    auto coef = berlekamp_massey(fib);
    assert((coef == std::vector<mint>{1, 1}));
    assert(linear_recurrence_kth(fib, coef, 10) == mint(55));
    assert(bmbm(fib, 10) == mint(55));

    poly p{0, 1};
    poly q{1, -1, -1};
    assert(bostan_mori(p, q, 10) == mint(55));
}
