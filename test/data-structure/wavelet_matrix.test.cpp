#include <algorithm>
#include <cassert>
#include <limits>
#include <random>
#include <utility>
#include <vector>
#include "../../src/data-structure/wavelet_matrix.hpp"

int main() {
    std::mt19937 rng(459);
    for (int n : {0, 1, 2, 31, 63, 64, 65, 127, 128, 129, 250}) {
        for (int trial = 0; trial < 8; trial++) {
            std::vector<long long> a(n);
            for (auto& x : a) x = int(rng() % 41) - 20;
            if (n >= 2) {
                a[0] = std::numeric_limits<long long>::min();
                a[1] = std::numeric_limits<long long>::max();
            }
            wavelet_matrix<long long> wm(a);
            assert(wm.size() == n);
            for (int i = 0; i < n; i++) assert(wm[i] == a[i] && wm.access(i) == a[i]);
            for (int q = 0; q < 500; q++) {
                int l = rng() % (n + 1), r = rng() % (n + 1);
                if (l > r) std::swap(l, r);
                long long x = int(rng() % 61) - 30, y = int(rng() % 61) - 30;
                if (x > y) std::swap(x, y);
                std::vector<long long> v(a.begin() + l, a.begin() + r);
                std::sort(v.begin(), v.end());
                int lo = std::lower_bound(v.begin(), v.end(), x) - v.begin();
                int hi = std::lower_bound(v.begin(), v.end(), y) - v.begin();
                assert(wm.range_freq(l, r, x) == lo);
                assert(wm.range_freq(l, r, x, y) == hi - lo);
                assert(wm.count(l, r, x) == std::count(v.begin(), v.end(), x));
                assert(wm.rank(x, r) == std::count(a.begin(), a.begin() + r, x));
                auto prev = wm.prev_value(l, r, x);
                auto next = wm.next_value(l, r, x);
                assert(bool(prev) == (lo > 0));
                assert(bool(next) == (lo < int(v.size())));
                if (prev) assert(*prev == v[lo - 1]);
                if (next) assert(*next == v[lo]);
                if (!v.empty()) {
                    int k = rng() % v.size();
                    assert(wm.kth_smallest(l, r, k) == v[k]);
                    assert(wm.kth_largest(l, r, k) == v[v.size() - 1 - k]);
                }
            }
            wm.build({});
            assert(wm.size() == 0 && wm.range_freq(0, 0, 0) == 0);
            wm.build(a);
            for (int i = 0; i < n; i++) assert(wm.access(i) == a[i]);
        }
    }
    using P = std::pair<int, int>;
    wavelet_matrix<P> pairs(std::vector<P>{{1, 2}, {-1, 4}, {1, 2}, {1, 1}});
    assert(pairs.kth_smallest(0, 4, 1) == P(1, 1));
    assert(pairs.count(0, 4, P(1, 2)) == 2);
    assert(pairs.prev_value(0, 4, P(1, 2)) == P(1, 1));
}
