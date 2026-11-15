#include <algorithm>
#include <cassert>
#include <limits>
#include <random>
#include <vector>
#include "../../src/data-structure/binary_trie.hpp"

int main() {
    std::mt19937 rng(459);
    for (int tc = 0; tc < 100; tc++) {
        binary_trie<int, 7> trie;
        std::vector<int> a;
        for (int q = 0; q < 1000; q++) {
            int x = rng() % 256, op = rng() % 4;
            if (op == 0 || a.empty()) {
                trie.insert(x);
                a.push_back(x);
            } else if (op == 1) {
                int i = rng() % a.size();
                trie.erase(a[i]);
                a.erase(a.begin() + i);
            } else if (op == 2) {
                trie.xor_all(x);
                for (int& v : a) v ^= x;
            } else {
                trie.xor_all(x);
                for (int& v : a) v ^= x;
            }
            std::sort(a.begin(), a.end());
            assert(trie.size() == int(a.size()));
            assert(trie.empty() == a.empty());
            assert(trie.count(x) == int(std::count(a.begin(), a.end(), x)));
            assert(trie.contains(x) == (std::count(a.begin(), a.end(), x) != 0));
            for (int upper : {-1, 0, x, 255, 256, 1000}) {
                assert(trie.count_less(upper) == int(std::lower_bound(a.begin(), a.end(), upper) - a.begin()));
            }
            if (!a.empty()) {
                assert(trie.min_element() == a.front());
                assert(trie.max_element() == a.back());
                int k = rng() % a.size();
                assert(trie.kth_smallest(k) == a[k]);
            }
        }
    }
    binary_trie<unsigned long long, 63> wide;
    auto high = 1ULL << 63;
    wide.insert(0);
    wide.insert(high);
    wide.insert(~0ULL);
    wide.xor_all(high);
    assert(wide.min_element() == 0);
    assert(wide.max_element() == high);
    assert(wide.count_less(~0ULL) == 3);
    wide.erase(0);
    assert(!wide.contains(0));
    binary_trie<int> signed_default;
    signed_default.insert(std::numeric_limits<int>::max());
    assert(signed_default.max_element() == std::numeric_limits<int>::max());
}
