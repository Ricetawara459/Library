#define PROBLEM "https://atcoder.jp/contests/abc403/tasks/abc403_g"

#include <iostream>

#include "../../src/data-structure/binary_trie.hpp"

int main() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);

    int q;
    std::cin >> q;

    binary_trie<int, 29> bt;
    bt.reserve(q * 31 + 1);

    long long answer = 0;
    while (q--) {
        long long y;
        std::cin >> y;
        int x = int((y + answer) % 1000000000) + 1;
        bt.insert(x);
        answer = bt.odd_position_sum();
        std::cout << answer << '\n';
    }
    return 0;
}
