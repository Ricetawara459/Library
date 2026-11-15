#include <algorithm>
#include <cassert>
#include <functional>
#include <random>
#include <set>
#include <string>
#include <vector>
#include "../../src/data-structure/ordered_set.hpp"

template <class Compare>
void test_random(Compare cmp) {
    ordered_set<int, Compare> tree(cmp);
    std::set<int, Compare> reference(cmp);
    std::mt19937 rng(459);
    for (int q = 0; q < 50000; q++) {
        int x = int(rng() % 401) - 200;
        if (rng() % 2) {
            bool a = tree.insert(x);
            bool b = reference.insert(x).second;
            assert(a == b);
        } else {
            bool a = tree.erase(x);
            bool b = reference.erase(x) != 0;
            assert(a == b);
        }
        assert(tree.size() == int(reference.size()));
        assert(tree.empty() == reference.empty());
        assert(tree.contains(x) == (reference.count(x) != 0));
        assert(tree.lower_bound(x) == std::distance(reference.begin(), reference.lower_bound(x)));
        assert(tree.upper_bound(x) == std::distance(reference.begin(), reference.upper_bound(x)));
        if (!reference.empty()) {
            int k = rng() % reference.size();
            auto it = reference.begin();
            std::advance(it, k);
            assert(tree.kth_smallest(k) == *it);
        }
    }
    auto copy = tree;
    tree.insert(1000);
    assert(!copy.contains(1000));
}

int main() {
    test_random(std::less<int>());
    test_random(std::greater<int>());
    ordered_set<std::string> words;
    assert(words.insert("banana"));
    assert(words.insert("apple"));
    assert(!words.insert("apple"));
    assert(words.kth_smallest(0) == "apple");
    assert(words.lower_bound("banana") == 1);
    assert(words.upper_bound("banana") == 2);
    assert(words.erase("apple"));
    assert(!words.erase("apple"));
}
