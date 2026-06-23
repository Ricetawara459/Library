#define PROBLEM "https://yukicoder.me/problems/no/186"

#include <iostream>
#include <vector>
#include "../../src/math/math.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    vector<long long> x(3), y(3);
    for (int i = 0; i < 3; i++) {
        cin >> x[i] >> y[i];
    }

    auto ans = crt(x, y);

    if (ans.second == 0) {
        // 解が存在しない場合
        cout << -1 << "\n";
    } else {
        long long res = ans.first;
        // 問題文の要求が「最小の正整数」なので、最小の非負整数が 0 の場合は、周期（LCM）を足して正にする
        if (res == 0) {
            res = ans.second;
        }
        cout << res << "\n";
    }

    return 0;
}