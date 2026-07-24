#define PROBLEM "https://judge.yosupo.jp/problem/static_range_frequency"

#include <algorithm>
#include <iostream>
#include <vector>

#include "../../src/misc/compress.hpp"

int main() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);

    int n, q;
    std::cin >> n >> q;
    std::vector<long long> a(n);
    compress<long long> cmp;
    for (auto& x : a) {
        std::cin >> x;
        cmp.add(x);
    }
    cmp.build();

    std::vector<std::vector<int>> positions(cmp.size());
    for (int i = 0; i < n; i++) positions[cmp(a[i])].push_back(i);

    while (q--) {
        int l, r;
        long long x;
        std::cin >> l >> r >> x;

        int id = cmp.get_id(x);
        if (id == cmp.size() || cmp[id] != x) {
            std::cout << 0 << '\n';
            continue;
        }
        auto first = std::lower_bound(
            positions[id].begin(), positions[id].end(), l);
        auto last = std::lower_bound(
            positions[id].begin(), positions[id].end(), r);
        std::cout << last - first << '\n';
    }
    return 0;
}
