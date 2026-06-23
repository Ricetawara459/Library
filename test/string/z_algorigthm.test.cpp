#define PROBLEM "https://judge.yosupo.jp/problem/zalgorithm"

#include <iostream>
#include <string>
#include <vector>
#include "../../src/string/z_algorithm.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    string s;
    if (!(cin >> s)) return 0;

    auto ans = z_algorithm(s);

    int n = int(ans.size());
    for (int i = 0; i < n; i++) {
        cout << ans[i] << (i + 1 == n ? "" : " ");
    }
    cout << "\n";

    return 0;
}