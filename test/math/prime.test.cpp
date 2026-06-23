#define PROBLEM "https://judge.yosupo.jp/problem/factorize"

#include <iostream>
#include <vector>
#include "../../src/math/prime.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int q;
    if (!(cin >> q)) return 0;

    while (q--) {
        long long a;
        cin >> a;
        
        auto res = factorize(a);
        
        cout << res.size();
        for (long long p : res) {
            cout << " " << p;
        }
        cout << "\n";
    }

    return 0;
}