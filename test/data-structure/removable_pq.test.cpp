#define PROBLEM "https://judge.yosupo.jp/problem/double_ended_priority_queue"

#include <iostream>
#include <vector>
#include "../../src/data-structure/removable_pq.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, q;
    if (!(cin >> n >> q)) return 0;

    removable_pq<int, less<int>> max_pq;     // 最大値管理用
    removable_pq<int, greater<int>> min_pq;  // 最小値管理用

    for (int i = 0; i < n; i++) {
        int s;
        cin >> s;
        max_pq.push(s);
        min_pq.push(s);
    }

    for (int i = 0; i < q; i++) {
        int type;
        cin >> type;
        if (type == 0) {
            int x;
            cin >> x;
            max_pq.push(x);
            min_pq.push(x);
        } else if (type == 1) {
            // 最小値の取得と削除
            int val = min_pq.top();
            cout << val << "\n";
            max_pq.erase(val); // 最大値側から遅延削除
            min_pq.pop();
        } else if (type == 2) {
            // 最大値の取得と削除
            int val = max_pq.top();
            cout << val << "\n";
            min_pq.erase(val); // 最小値側から遅延削除
            max_pq.pop();
        }
    }

    return 0;
}