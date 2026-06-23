#define PROBLEM "https://judge.yosupo.jp/problem/pow_of_matrix"

#include <iostream>
#include "../../src/math/modint.hpp"
#include "../../src/math/matrix.hpp"

using namespace std;

// Library Checker の指定 mod である 998244353 を使用
using mint = modint998244353;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n;
    long long k;
    if (!(cin >> n >> k)) return 0;

    matrix<mint> a(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            long long val;
            cin >> val;
            a[i][j] = val;
        }
    }

    // 行列累乗の実行
    auto ans = a.pow(k);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << ans[i][j].val() << (j == n - 1 ? "" : " ");
        }
        cout << "\n";
    }

    return 0;
}