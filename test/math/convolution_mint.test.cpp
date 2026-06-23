#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#include <iostream>
#include <vector>
#include "../../src/math/modint.hpp"       // ご自身の modint のパス
#include "../../src/math/convolution.hpp"

using namespace std;

// ご自身の環境の Mod 998244353 用の型名に書き換えてください
// 例: using mint = modint998244353;
using mint = modint998244353; 

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    vector<mint> a(n), b(m);
    for (int i = 0; i < n; i++) {
        long long x;
        cin >> x;
        a[i] = mint(x);
    }
    for (int i = 0; i < m; i++) {
        long long x;
        cin >> x;
        b[i] = mint(x);
    }

    // template版の convolution が自動で選択される
    auto c = convolution(a, b);

    for (int i = 0; i < n + m - 1; i++) {
        cout << c[i].val() << (i == n + m - 2 ? "" : " ");
    }
    cout << "\n";

    return 0;
}