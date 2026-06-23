#line 2 "src/data-structure/partition_sum.hpp"

#include <vector>
#include <cassert>
#include <utility>
#include <functional>
#line 2 "src/data-structure/removable_pq.hpp"

#include <queue>
#line 6 "src/data-structure/removable_pq.hpp"

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
#line 8 "src/data-structure/partition_sum.hpp"

template <typename T, typename SumT = long long>
struct partition_sum {
  public:
    removable_pq<T, std::less<T>> que_l;          // 最大値ヒープ
    removable_pq<T, std::greater<T>> que_r;       // 最小値ヒープ

    void insert(T x) {
        if (que_l.empty() || que_l.top() < x) {
            que_r.push(x);
        } else {
            que_l.push(x);
        }
    }

    void erase(T x) {
        if (que_l.empty() || que_l.top() < x) {
            que_r.erase(x);
        } else {
            que_l.erase(x);
        }
    }

    // 小さい方から K 個を左 (que_l) に集め、{左の総和, 右の総和} を返す
    std::pair<SumT, SumT> get_sum(int K) {
        assert(0 <= K && K <= que_l.cnt + que_r.cnt);
        while (que_l.cnt < K) {
            que_l.push(que_r.top());
            que_r.pop();
        }
        while (que_l.cnt > K) {
            que_r.push(que_l.top());
            que_l.pop();
        }
        return {que_l.sum, que_r.sum};
    }
};
