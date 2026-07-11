#ifndef LIBRARY_DATA_STRUCTURE_PARTITION_SUM_HPP
#define LIBRARY_DATA_STRUCTURE_PARTITION_SUM_HPP

#line 1 "src/data-structure/partition_sum.hpp"

#include <vector>
#include <cassert>
#include <utility>
#include <functional>
#ifndef LIBRARY_DATA_STRUCTURE_REMOVABLE_PQ_HPP
#define LIBRARY_DATA_STRUCTURE_REMOVABLE_PQ_HPP

#line 1 "src/data-structure/removable_pq.hpp"

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

#endif  // LIBRARY_DATA_STRUCTURE_REMOVABLE_PQ_HPP


template <typename T, typename SumT = long long>
struct partition_sum {
  public:
    removable_pq<T, std::less<T>> que_l;          // 最大値ヒープ
    removable_pq<T, std::greater<T>> que_r;       // 最小値ヒープ

    /// 値 x を追加する。
    void insert(T x) {
        if (que_l.empty() || que_l.top() < x) {
            que_r.push(x);
        } else {
            que_l.push(x);
        }
    }

    /// 値 x を 1 個削除する。
    void erase(T x) {
        if (que_l.empty() || que_l.top() < x) {
            que_r.erase(x);
        } else {
            que_l.erase(x);
        }
    }

    // 小さい方から K 個を左 (que_l) に集め、{左の総和, 右の総和} を返す
    /// 小さい方から K 個の和と、残りの和を返す。
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

#endif  // LIBRARY_DATA_STRUCTURE_PARTITION_SUM_HPP
