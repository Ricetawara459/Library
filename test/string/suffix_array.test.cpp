#define PROBLEM "https://judge.yosupo.jp/problem/suffixarray"

#include <iostream>
#include <string>
#include "../../src/string/suffix_array.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    string s;
    if (!(cin >> s)) return 0;

    auto sa = suffix_array(s);
    int n = int(sa.size());
    for (int i = 0; i < n; i++) {
        cout << sa[i] << (i + 1 == n ? "" : " ");
    }
    cout << "\n";

    return 0;
}