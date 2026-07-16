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
    convex_hull_trick<long long> cht;
    std::vector<std::pair<long long, long long>> lines = {
        {2, 3}, {-1, 10}, {4, -5}, {2, 8}, {-3, 7},
    };
    for (auto [a, b] : lines) cht.add_line(a, b);
    for (long long x = -20; x <= 20; x++) {
        long long want = std::numeric_limits<long long>::min();
        for (auto [a, b] : lines) want = std::max(want, eval(a, b, x));
        assert(cht.max(x) == want);
    }

    convex_hull_trick<long long> mn;
    for (auto [a, b] : lines) mn.add_min_line(a, b);
    for (long long x = -20; x <= 20; x++) {
        long long want = std::numeric_limits<long long>::max();
        for (auto [a, b] : lines) want = std::min(want, eval(a, b, x));
        assert(mn.min(x) == want);
    }

    int n;
    if (std::cin >> n) {
        convex_hull_trick<long long> q;
        for (int i = 0; i < n; i++) {
            long long a, b;
            std::cin >> a >> b;
            q.add_line(a, b);
        }
        int m;
        std::cin >> m;
        while (m--) {
            long long x;
            std::cin >> x;
            std::cout << q.max(x) << '\n';
        }
    }
}
