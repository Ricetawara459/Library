#include <cassert>
#include <functional>
#include <iterator>
#include <random>
#include <set>
#include <utility>
#include "../../src/data-structure/ordered_multiset.hpp"

template <class Compare>
void test_pairs(Compare cmp) {
    using P = std::pair<int, int>;
    ordered_multiset<P, Compare> ms(cmp);
    ordered_set<P, Compare> set(cmp);
    std::multiset<P, Compare> ref(cmp);
    std::set<P, Compare> ref_set(cmp);
    std::mt19937 rng(459);
    for (int q = 0; q < 20000; q++) {
        P x{int(rng() % 13) - 6, int(rng() % 13) - 6};
        if (rng() % 2) {
            ms.insert(x);
            ref.insert(x);
            bool inserted = set.insert(x);
            assert(inserted == ref_set.insert(x).second);
        } else {
            bool erased = ms.erase(x);
            auto it = ref.find(x);
            bool expected = it != ref.end();
            if (expected) ref.erase(it);
            assert(erased == expected);
            bool set_erased = set.erase(x);
            assert(set_erased == (ref_set.erase(x) != 0));
        }
        assert(ms.size() == int(ref.size()));
        assert(ms.empty() == ref.empty());
        assert(ms.count(x) == int(ref.count(x)));
        assert(ms.contains(x) == (ref.count(x) != 0));
        assert(ms.lower_bound(x) == std::distance(ref.begin(), ref.lower_bound(x)));
        assert(ms.upper_bound(x) == std::distance(ref.begin(), ref.upper_bound(x)));
        assert(set.lower_bound(x) == std::distance(ref_set.begin(), ref_set.lower_bound(x)));
        if (!ref.empty()) {
            int k = rng() % ref.size();
            auto it = ref.begin();
            std::advance(it, k);
            assert(ms.kth_smallest(k) == *it);
        }
        if (!ref_set.empty()) {
            int k = rng() % ref_set.size();
            auto it = ref_set.begin();
            std::advance(it, k);
            assert(set.kth_smallest(k) == *it);
        }
    }
    auto copy = ms;
    ms.insert({100, 200});
    assert(!copy.contains({100, 200}));
}

int main() {
    using P = std::pair<int, int>;
    test_pairs(std::less<P>());
    test_pairs(std::greater<P>());
    ordered_multiset<int> ms;
    ms.insert(2); ms.insert(2); ms.insert(1);
    assert(ms.size() == 3 && ms.count(2) == 2);
    assert(ms.kth_smallest(0) == 1 && ms.kth_smallest(2) == 2);
    assert(ms.lower_bound(2) == 1 && ms.upper_bound(2) == 3);
    assert(ms.erase(2) && ms.count(2) == 1);
    assert(ms.erase(2) && !ms.erase(2));
}
