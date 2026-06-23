#define PROBLEM "https://judge.yosupo.jp/problem/sum_of_floor_of_linear"

#include <iostream>
#include "../../src/math/math.hpp"

using namespace std;

void solve() {
    long long n, m, a, b;
    cin >> n >> m >> a >> b;
    cout << floor_sum(n, m, a, b) << "\n";
}

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    int t;
    if (cin >> t) {
        while (t--) solve();
    }
    return 0;
}