#ifndef LIBRARY_DATA_STRUCTURE_PARTITION_SUM_HPP
#define LIBRARY_DATA_STRUCTURE_PARTITION_SUM_HPP

#pragma once

#include <vector>
#include <cassert>
#include <utility>
#include <functional>
#include "removable_pq.hpp"

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
