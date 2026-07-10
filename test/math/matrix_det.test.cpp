#define PROBLEM "https://judge.yosupo.jp/problem/matrix_det"

#include <iostream>

#include "../../src/math/matrix.hpp"
#include "../../src/math/modint.hpp"

using namespace std;
using mint = modint998244353;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n;
    cin >> n;
    matrix<mint> a(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            long long x;
            cin >> x;
            a[i][j] = x;
        }
    }
    cout << a.determinant().val() << '\n';
}
