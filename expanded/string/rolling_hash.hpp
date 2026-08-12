#ifndef LIBRARY_STRING_ROLLING_HASH_HPP
#define LIBRARY_STRING_ROLLING_HASH_HPP

#line 1 "src/string/rolling_hash.hpp"


#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <type_traits>

struct rolling_hash {
    using u64 = unsigned long long;
    using u128 = unsigned __int128;
    static const u64 MOD = (1ULL << 61) - 1;

    // すべてのインスタンスで共通のランダム基数を取得
    static u64 get_base() {
        static u64 base = 0;
        if (base == 0) {
            std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
            base = rng() % (MOD - 2) + 2;
        }
        return base;
    }

    static inline std::vector<u64> pow_table{1};

    // 必要に応じて基数の冪乗テーブルを自動拡張
    static void extend_pow(int n) {
        if ((int)pow_table.size() > n) return;
        int old_sz = pow_table.size();
        pow_table.resize(n + 1);
        u64 b = get_base();
        for (int i = old_sz; i <= n; i++) {
            pow_table[i] = mul(pow_table[i - 1], b);
        }
    }

    // Mod 2^61 - 1 の高速な乗算アルゴリズム
    static u64 mul(u128 a, u128 b) {
        u128 ans = a * b;
        u64 res = (ans >> 61) + (ans & MOD);
        if (res >= MOD) res -= MOD;
        return res;
    }

    std::vector<u64> hash_table{0};

    rolling_hash() = default;
    template <typename T>
    rolling_hash(const T& s) { build(s); }

    /// 列 s の rolling hash を構築する。
    template <typename T>
    void build(const T& s) {
        int n = int(s.size());
        hash_table.resize(n + 1, 0);
        extend_pow(n);
        u64 b = get_base();
        for (int i = 0; i < n; i++) {
            u64 value;
            if constexpr (std::is_same_v<std::decay_t<decltype(s[i])>, char>) {
                value = static_cast<unsigned char>(s[i]);
            } else {
                value = static_cast<u64>(s[i]);
            }
            hash_table[i + 1] = mul(hash_table[i], b) + value;
            if (hash_table[i + 1] >= MOD) hash_table[i + 1] -= MOD;
        }
    }

    /// 半開区間 [l, r) のハッシュを返す。
    u64 get(int l, int r) const {
        u64 res = hash_table[r] + MOD - mul(hash_table[l], pow_table[r - l]);
        if (res >= MOD) res -= MOD;
        return res;
    }

    /// 現在の列の長さを返す。
    int size() const {
        return int(hash_table.size()) - 1;
    }

    /// 現在の列の右に文字列 s を連結する。計算量 O(|s|)。
    void concat_right(const std::string& s) {
        int old_size = size();
        int add_size = int(s.size());
        extend_pow(old_size + add_size);
        hash_table.reserve(old_size + add_size + 1);
        u64 b = get_base();
        for (char c : s) {
            u64 next = mul(hash_table.back(), b) + static_cast<unsigned char>(c);
            if (next >= MOD) next -= MOD;
            hash_table.push_back(next);
        }
    }

    /// 現在の列の右に文字 c を連結する。計算量 O(1)。
    void concat_right(char c) {
        extend_pow(size() + 1);
        u64 next = mul(hash_table.back(), get_base()) + static_cast<unsigned char>(c);
        if (next >= MOD) next -= MOD;
        hash_table.push_back(next);
    }

    /// 現在の列の左に文字列 s を連結する。計算量 O(size() + |s|)。
    void concat_left(const std::string& s) {
        if (s.empty()) return;
        int left_size = int(s.size());
        int old_size = size();
        extend_pow(left_size + old_size);

        std::vector<u64> next(left_size + old_size + 1, 0);
        u64 b = get_base();
        for (int i = 0; i < left_size; i++) {
            next[i + 1] = mul(next[i], b) + static_cast<unsigned char>(s[i]);
            if (next[i + 1] >= MOD) next[i + 1] -= MOD;
        }
        u64 left_hash = next[left_size];
        for (int i = 1; i <= old_size; i++) {
            next[left_size + i] = connect(left_hash, hash_table[i], i);
        }
        hash_table.swap(next);
    }

    /// 現在の列の左に文字 c を連結する。計算量 O(size())。
    void concat_left(char c) {
        int old_size = size();
        extend_pow(old_size + 1);
        u64 value = static_cast<unsigned char>(c);
        hash_table.resize(old_size + 2);
        for (int i = old_size; i >= 1; i--) {
            hash_table[i + 1] = connect(value, hash_table[i], i);
        }
        hash_table[1] = value;
    }

    /// ハッシュ h1 の後ろに長さ len2 の列のハッシュ h2 を連結する。
    static u64 connect(u64 h1, u64 h2, int len2) {
        extend_pow(len2);
        u64 res = mul(h1, pow_table[len2]) + h2;
        if (res >= MOD) res -= MOD;
        return res;
    }
};

#endif  // LIBRARY_STRING_ROLLING_HASH_HPP
