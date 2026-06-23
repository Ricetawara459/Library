#define PROBLEM "https://judge.yosupo.jp/problem/sum_of_products"

#include <iostream>
#include <vector>
#include "../../src/math/modint.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n;
    if (!(cin >> n)) return 0;

    vector<mint998244353> a(n), b(n);
    for (int i = 0; i < n; i++) cin >> a[i];
    for (int i = 0; i < n; i++) cin >> b[i];

    mint998244353 ans = 0;
    for (int i = 0; i < n; i++) {
        ans += a[i] * b[i];
    }

    cout << ans << "\n";
    return 0;
}