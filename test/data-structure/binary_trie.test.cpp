#define PROBLEM "https://judge.yosupo.jp/problem/set_xor_min"

#include <iostream>

#include "../../src/data-structure/binary_trie.hpp"

int main() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);

    int q;
    std::cin >> q;
    binary_trie<int, 29> bt;
    while (q--) {
        int t, x;
        std::cin >> t >> x;
        if (t == 0) {
            if (!bt.contains(x)) bt.insert(x);
        } else if (t == 1) {
            if (bt.contains(x)) bt.erase(x);
        } else {
            std::cout << bt.min_xor(x) << '\n';
        }
    }
    return 0;
}
