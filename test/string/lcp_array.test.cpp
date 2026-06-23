#define PROBLEM "https://judge.yosupo.jp/problem/number_of_substrings"

#include <iostream>
#include <string>
#include <numeric>
#include "../../src/string/suffix_array.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    string s;
    if (!(cin >> s)) return 0;

    long long n = s.size();
    auto sa = suffix_array(s);
    auto lcp = lcp_array(s, sa);

    long long ans = n * (n + 1) / 2;
    for (int x : lcp) {
        ans -= x;
    }

    cout << ans << "\n";

    return 0;
}