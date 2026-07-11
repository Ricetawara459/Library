#define PROBLEM "https://judge.yosupo.jp/problem/kth_term_of_linearly_recurrent_sequence"

#include <iostream>
#include <vector>

#include "../../src/fps/fps.hpp"
#include "../../src/math/modint.hpp"

using namespace std;
using mint = modint998244353;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int d;
    long long k;
    cin >> d >> k;

    vector<mint> a(d), c(d);
    for (int i = 0; i < d; i++) cin >> a[i];
    for (int i = 0; i < d; i++) cin >> c[i];

    cout << linear_recurrence_kth(a, c, k).val() << '\n';
}
