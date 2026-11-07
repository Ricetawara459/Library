#include <cassert>
#include <random>
#include <vector>
#include "../../src/data-structure/fenwick_tree.hpp"

int main() {
    std::mt19937 rng(459);
    for (int n = 0; n <= 80; n++) {
        fenwick_tree<long long> fw(n);
        std::vector<long long> a(n);
        for (int iteration = 0; iteration < 100; iteration++) {
            if (n) {
                int p = rng() % n;
                long long value = rng() % 8;
                fw.add(p, value - a[p]);
                a[p] = value;
            }
            long long total = 0;
            for (long long value : a) total += value;
            for (long long x = -1; x <= total + 1; x++) {
                int lower = 0, upper = 0;
                long long sum = 0;
                if (x > 0) {
                    while (lower < n && sum + a[lower] < x) sum += a[lower++];
                }
                sum = 0;
                if (x >= 0) {
                    while (upper < n && sum + a[upper] <= x) sum += a[upper++];
                }
                assert(fw.lower_bound(x) == lower);
                assert(fw.upper_bound(x) == upper);
            }
        }
    }
}
