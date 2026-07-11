#ifndef LIBRARY_STRING_ROLLING_HASH_HPP
#define LIBRARY_STRING_ROLLING_HASH_HPP

#line 1 "src/string/rolling_hash.hpp"


#include <vector>
#include <string>
#include <chrono>
#include <random>

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

    std::vector<u64> hash_table;

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
            hash_table[i + 1] = mul(hash_table[i], b) + s[i];
            if (hash_table[i + 1] >= MOD) hash_table[i + 1] -= MOD;
        }
    }

    /// 半開区間 [l, r) のハッシュを返す。
    u64 get(int l, int r) const {
        u64 res = hash_table[r] + MOD - mul(hash_table[l], pow_table[r - l]);
        if (res >= MOD) res -= MOD;
        return res;
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
