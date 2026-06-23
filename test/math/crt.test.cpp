#define PROBLEM "https://judge.yosupo.jp/problem/chinese_remainder_theorem"

#include <iostream>
#include <vector>
#include "../../src/math/math.hpp"

using namespace std;

void solve() {
    long long r1, m1, r2, m2;
    cin >> r1 >> m1 >> r2 >> m2;
    auto ans = crt({r1, r2}, {m1, m2});
    if (ans.second == 0) {
        cout << -1 << "\n";
    } else {
        cout << ans.first << " " << ans.second << "\n";
    }
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