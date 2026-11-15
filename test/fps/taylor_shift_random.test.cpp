#include <cassert>
#include <random>
#include <vector>
#include "../../src/fps/fps.hpp"

template <class Mint>
void test(int max_n) {
    std::mt19937 rng(459);
    for (int n = 0; n <= max_n; n++) {
        for (int trial = 0; trial < 8; trial++) {
            fps<Mint> f(n), expected(n);
            for (auto& x : f) x = int(rng() % 2001) - 1000;
            Mint c = trial == 0 ? Mint(0) : Mint(int(rng() % 21) - 10);
            std::vector<Mint> powers(n, 1), row(n, 0);
            for (int i = 1; i < n; i++) powers[i] = powers[i - 1] * c;
            if (n) row[0] = 1;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= i; j++) expected[j] += f[i] * row[j] * powers[i - j];
                for (int j = i + 1; j > 0 && j < n; j--) row[j] += row[j - 1];
            }
            auto shifted = f.taylor_shift(c);
            assert(shifted == expected);
            assert(shifted.taylor_shift(-c) == f);
            assert(shifted.eval(Mint(7)) == f.eval(Mint(7) + c));
        }
    }
}

int main() {
    test<modint998244353>(140);
    test<modint1000000007>(140);
    test<static_modint<17>>(17);
}
