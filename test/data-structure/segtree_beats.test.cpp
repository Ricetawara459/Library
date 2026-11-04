#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <random>
#include <vector>

void expect(bool ok) {
    if (!ok) std::abort();
}

#include "../../src/data-structure/segtree_beats.hpp"
#include "../../snippets/segtree_beats_chmin_chmax_sum.cpp"

using namespace beats_chmin_chmax_sum;

namespace all_prod_test {
struct node { ll sum; int len; };
int calls = 0;
node merge(node a, node b) { calls++; return {a.sum + b.sum, a.len + b.len}; }
node empty() { return {0, 0}; }
bool apply(int f, node& x) { calls++; x.sum += ll(f) * x.len; return true; }
int compose(int f, int g) { return f + g; }
int identity() { return 0; }
using tree = segtree_beats<node, merge, empty, int, apply, compose, identity>;
}

int main() {
    for (int n : {0, 1, 5, 255, 65535}) {
        all_prod_test::tree seg(std::vector<all_prod_test::node>(n, {1, 1}));
        seg.apply(0, n, 7);
        if (n) seg.set(n - 1, {11, 1});
        all_prod_test::calls = 0;
        const auto& view = seg;
        expect(view.all_prod().sum == 8LL * n + (n ? 3 : 0));
        expect(all_prod_test::calls == 0);
    }

    std::mt19937 rng(123456789);
    for (int tc = 0; tc < 200; tc++) {
        int n = 1 + int(rng() % 40);
        std::vector<ll> a(n);
        for (ll& x : a) x = int(rng() % 101) - 50;
        auto b = a;
        range_chmin_chmax_add_sum seg(a);

        for (int q = 0; q < 300; q++) {
            int l = int(rng() % n), r = int(rng() % n) + 1;
            if (l > r) std::swap(l, r);
            int typ = int(rng() % 7);
            ll x = int(rng() % 101) - 50;
            if (typ == 0) {
                seg.range_add(l, r, x);
                for (int i = l; i < r; i++) b[i] += x;
            } else if (typ == 1) {
                seg.range_chmin(l, r, x);
                for (int i = l; i < r; i++) b[i] = std::min(b[i], x);
            } else if (typ == 2) {
                seg.range_chmax(l, r, x);
                for (int i = l; i < r; i++) b[i] = std::max(b[i], x);
            } else if (typ == 3) {
                seg.set(l, x);
                b[l] = x;
            } else if (typ == 4) {
                ll expected = std::accumulate(b.begin() + l, b.begin() + r, 0LL);
                expect(seg.range_sum(l, r) == expected);
            } else if (typ == 5) {
                expect(seg.get(l) == b[l]);
            } else {
                ll expected = std::accumulate(b.begin(), b.end(), 0LL);
                expect(seg.all_sum() == expected);
            }
        }
        for (int l = 0; l <= n; l++) {
            for (int r = l; r <= n; r++) {
                ll expected = std::accumulate(b.begin() + l, b.begin() + r, 0LL);
                expect(seg.range_sum(l, r) == expected);
            }
        }
    }

    {
        std::vector<ll> a(30, 1);
        range_chmin_chmax_add_sum seg(a);
        seg.range_add(5, 25, 2);
        expect(seg.seg.max_right(0, [](S x) { return x.sum <= 20; }) == 10);
        expect(seg.seg.min_left(30, [](S x) { return x.sum <= 20; }) == 20);
    }
}
