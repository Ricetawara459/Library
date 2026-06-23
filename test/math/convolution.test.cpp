#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#include <iostream>
#include <vector>
#include "../../src/math/convolution.hpp"
#include "../../src/math/modint.hpp"

using mint=mint998244353;

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    vector<mint> a(n), b(m);
    for (int i = 0; i < n; i++) cin >> a[i];
    for (int i = 0; i < m; i++) cin >> b[i];

    auto c = convolution(a, b);

    for (int i = 0; i < n + m - 1; i++) {
        cout << c[i] << (i == n + m - 2 ? "" : " ");
    }
    cout << "\n";

    return 0;
}