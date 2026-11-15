#include <cassert>
#include <random>
#include <vector>
#include "../../src/fps/fps.hpp"

template <class Mint>
void test() {
    std::mt19937 rng(459);
    for (int tc = 0; tc < 300; tc++) {
        int n = tc % 100 == 0 ? 500 : rng() % 160;
        int m = tc % 100 == 0 ? 201 : rng() % 130;
        fps<Mint> f(n), g(m);
        for (auto& x : f) x = int(rng() % 19) - 9;
        for (auto& x : g) x = int(rng() % 19) - 9;
        if (tc % 3 == 0 && !f.empty()) f.back() = 0;
        if (tc % 4 == 0 && !g.empty()) g.back() = 0;
        auto normalized = g;
        normalized.shrink();
        if (normalized.empty()) { g.push_back(1); normalized = g; normalized.shrink(); }

        fps<Mint> remainder = f;
        remainder.shrink();
        fps<Mint> quotient;
        if (remainder.size() >= normalized.size()) {
            quotient.resize(remainder.size() - normalized.size() + 1);
            Mint inv = normalized.back().inv();
            while (!remainder.empty() && remainder.size() >= normalized.size()) {
                int shift = int(remainder.size() - normalized.size());
                Mint value = remainder.back() * inv;
                quotient[shift] = value;
                for (int i = 0; i < int(normalized.size()); i++) remainder[i + shift] -= value * normalized[i];
                remainder.shrink();
            }
            quotient.shrink();
        }
        auto qr = f.divmod(g);
        assert(qr.first == quotient && qr.second == remainder);
        assert(f / g == quotient && f % g == remainder);
        auto a = f, b = f;
        a /= g; b %= g;
        assert(a == quotient && b == remainder);
        auto reconstructed = quotient * normalized + remainder;
        reconstructed.shrink();
        auto original = f; original.shrink();
        assert(reconstructed == original);

        std::vector<Mint> xs(m);
        for (auto& x : xs) x = tc % 2 ? Mint(static_cast<long long>(rng())) : Mint(int(rng() % 13) - 6);
        auto actual = f.multipoint_eval(xs);
        assert(actual.size() == xs.size());
        for (int i = 0; i < m; i++) assert(actual[i] == f.eval(xs[i]));
    }
    fps<Mint> f{1, 2, 3};
    assert(f.multipoint_eval({}).empty());
    assert((fps<Mint>{}.multipoint_eval({0, 1, 1}) == std::vector<Mint>{0, 0, 0}));
    auto self = f;
    self /= self;
    assert((self == fps<Mint>{1}));
    self = f;
    self %= self;
    assert(self.empty());
    for (auto [n, m] : std::vector<std::pair<int, int>>{{513, 256}, {1025, 128}, {127, 1025}, {2048, 1024}, {123, 4096}}) {
        fps<Mint> poly(n);
        std::vector<Mint> xs(m);
        for (auto& x : poly) x = static_cast<long long>(rng());
        for (auto& x : xs) x = static_cast<long long>(rng());
        auto result = poly.multipoint_eval(xs);
        for (int i = 0; i < m; i++) assert(result[i] == poly.eval(xs[i]));
    }
}

int main() {
    test<modint998244353>();
    test<modint1000000007>();
}
