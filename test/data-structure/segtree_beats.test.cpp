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

int main() {
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
