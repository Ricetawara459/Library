#define PROBLEM "https://judge.yosupo.jp/problem/range_kth_smallest"

#include <iostream>
#include <vector>
#include "../../src/data-structure/wavelet_matrix.hpp"

int main() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);
    int n, q;
    if (!(std::cin >> n >> q)) return 0;
    std::vector<int> a(n);
    for (int& x : a) std::cin >> x;
    wavelet_matrix<int> wm(a);
    while (q--) {
        int l, r, k;
        std::cin >> l >> r >> k;
        std::cout << wm.kth_smallest(l, r, k) << '\n';
    }
}
