#line 2 "src/data-structure/removable_pq.hpp"

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

    /// 現在残っている要素が空なら true。
    bool empty() const { return cnt == 0; }
    int size() const { return cnt; }

    /// 現在残っている要素の先頭を返す。削除予約済み要素は自動で捨てる。
    T top() {
        normalize();
        return que.top();
    }

    /// 値 x を追加する。
    void push(T x) {
        que.push(x);
        cnt += 1;
        sum += x;
    }

    // 任意の要素 x を1個指定して削除
    /// 値 x を 1 個削除予約する。存在する値に対して呼ぶ。
    void erase(T x) {
        rm_que.push(x);
        cnt -= 1;
        sum -= x;
    }

    // top の要素（最大値または最小値）を削除
    /// 現在の先頭要素を 1 個削除する。
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
