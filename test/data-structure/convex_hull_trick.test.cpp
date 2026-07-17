#include <cassert>
#include <algorithm>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "../../src/data-structure/convex_hull_trick.hpp"

long long eval(long long a, long long b, long long x) {
    return a * x + b;
}

int main() {
    std::vector<std::pair<long long, long long>> lines = {
        {2, 3}, {-1, 10}, {4, -5}, {2, 8}, {-3, 7},
    };

    cht_max<long long> mx;
    for (auto [a, b] : lines) mx.add(a, b);
    for (long long x = -20; x <= 20; x++) {
        long long want = std::numeric_limits<long long>::min();
        for (auto [a, b] : lines) want = std::max(want, eval(a, b, x));
        assert(mx.query(x) == want);
    }

    cht_min<long long> mn;
    for (auto [a, b] : lines) mn.add(a, b);
    for (long long x = -20; x <= 20; x++) {
        long long want = std::numeric_limits<long long>::max();
        for (auto [a, b] : lines) want = std::min(want, eval(a, b, x));
        assert(mn.query(x) == want);
    }

    int n;
    if (std::cin >> n) {
        cht_max<long long> q;
        for (int i = 0; i < n; i++) {
            long long a, b;
            std::cin >> a >> b;
            q.add(a, b);
        }
        int m;
        std::cin >> m;
        while (m--) {
            long long x;
            std::cin >> x;
            std::cout << q.query(x) << '\n';
        }
    }
}
