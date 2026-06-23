#define PROBLEM "https://judge.yosupo.jp/problem/matrix_product"

#include <iostream>
#include <vector>
#include "../../src/math/modint.hpp"

using namespace std;
using mint = modint998244353;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, m, k;
    if (!(cin >> n >> m >> k)) return 0;

    vector<vector<mint>> a(n, vector<mint>(m));
    vector<vector<mint>> b(m, vector<mint>(k));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) cin >> a[i][j];
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) cin >> b[i][j];
    }

    vector<vector<mint>> c(n, vector<mint>(k, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int h = 0; h < k; h++) {
                c[i][h] += a[i][j] * b[j][h];
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            cout << c[i][j] << (j + 1 == k ? "" : " ");
        }
        cout << "\n";
    }

    return 0;
}