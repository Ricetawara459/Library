#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_14_B"

#include <iostream>
#include <string>
#include <vector>
#include "../../src/string/rolling_hash.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    string t, p;
    if (!(cin >> t >> p)) return 0;

    if (t.size() < p.size()) return 0;

    // テキスト T とパターン P のハッシュをそれぞれ構築
    rolling_hash rh_t(t);
    rolling_hash rh_p(p);

    unsigned long long target_hash = rh_p.get(0, p.size());
    int len_t = t.size();
    int len_p = p.size();

    // T の中から P と同じ長さの区間をスライドしながらハッシュを比較
    for (int i = 0; i <= len_t - len_p; i++) {
        if (rh_t.get(i, i + len_p) == target_hash) {
            cout << i << "\n";
        }
    }

    return 0;
}