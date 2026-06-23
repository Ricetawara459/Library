#pragma once

#include <queue>
#include <vector>
#include <functional>

template <typename T, typename Compare = std::less<T>>
struct removable_pq {
  public:
    std::priority_queue<T, std::vector<T>, Compare> que, rm_que;
    int cnt = 0;
    T sum = 0;

    removable_pq() = default;

    bool empty() const { return cnt == 0; }
    int size() const { return cnt; }

    T top() {
        normalize();
        return que.top();
    }

    void push(T x) {
        que.push(x);
        cnt += 1;
        sum += x;
    }

    // 任意の要素 x を1個指定して削除
    void erase(T x) {
        rm_que.push(x);
        cnt -= 1;
        sum -= x;
    }

    // top の要素（最大値または最小値）を削除
    void pop() {
        normalize();
        cnt -= 1;
        sum -= que.top();
        que.pop();
    }

  private:
    void normalize() {
        while (!rm_que.empty() && rm_que.top() == que.top()) {
            rm_que.pop();
            que.pop();
        }
    }
};