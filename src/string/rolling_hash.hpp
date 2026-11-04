#ifndef LIBRARY_STRING_ROLLING_HASH_HPP
#define LIBRARY_STRING_ROLLING_HASH_HPP

#pragma once

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

    // 右側の累積ハッシュ。左追加分は left_hash_table で別に管理する。
    std::vector<u64> hash_table{0};

    rolling_hash() = default;
    template <typename T>
    rolling_hash(const T& s) { build(s); }

    /// 列 s の rolling hash を構築する。
    template <typename T>
    void build(const T& s) {
        int n = int(s.size());
        left_hash_table.assign(1, 0);
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

    /// 半開区間 [l, r) のハッシュを O(1) で返す。左右連結後も使える。
    u64 get(int l, int r) const {
        int left_size = int(left_hash_table.size()) - 1;
        u64 res;
        if (r <= left_size) {
            res = mul(left_hash_table[left_size - l], pow_table[r - l]) + MOD
                  - left_hash_table[left_size - r];
        } else if (l >= left_size) {
            res = hash_table[r - left_size] + MOD
                  - mul(hash_table[l - left_size], pow_table[r - l]);
        } else {
            res = mul(left_hash_table[left_size - l], pow_table[r - l])
                  + hash_table[r - left_size];
        }
        if (res >= MOD) res -= MOD;
        return res;
    }

    /// 現在の列の長さを返す。
    int size() const {
        return int(hash_table.size() + left_hash_table.size()) - 2;
    }

    /// 現在の列の右に文字列 s を連結する。償却 O(|s|)。
    void concat_right(const std::string& s) {
        int old_size = size();
        int add_size = int(s.size());
        extend_pow(old_size + add_size);
        u64 b = get_base();
        for (char c : s) {
            u64 next = mul(hash_table.back(), b) + static_cast<unsigned char>(c);
            if (next >= MOD) next -= MOD;
            hash_table.push_back(next);
        }
    }

    /// 現在の列の右に文字 c を連結する。償却 O(1)。
    void concat_right(char c) {
        extend_pow(size() + 1);
        u64 next = mul(hash_table.back(), get_base()) + static_cast<unsigned char>(c);
        if (next >= MOD) next -= MOD;
        hash_table.push_back(next);
    }

    /// 現在の列の左に文字列 s を連結する。償却 O(|s|)。
    void concat_left(const std::string& s) {
        for (auto it = s.rbegin(); it != s.rend(); ++it) concat_left(*it);
    }

    /// 現在の列の左に文字 c を連結する。償却 O(1)。
    void concat_left(char c) {
        extend_pow(size() + 1);
        u64 next = left_hash_table.back() + static_cast<unsigned char>(c);
        if (next >= MOD) next -= MOD;
        left_hash_table.push_back(mul(next, get_base_inverse()));
    }

    /// ハッシュ h1 の後ろに長さ len2 の列のハッシュ h2 を連結する。
    static u64 connect(u64 h1, u64 h2, int len2) {
        extend_pow(len2);
        u64 res = mul(h1, pow_table[len2]) + h2;
        if (res >= MOD) res -= MOD;
        return res;
    }

  private:
    // k 文字時点の左追加分のハッシュ / base^k を蓄積する。
    std::vector<u64> left_hash_table{0};

    static u64 get_base_inverse() {
        static const u64 inverse = [] {
            u64 res = 1, b = get_base();
            for (u64 e = MOD - 2; e > 0; e >>= 1) {
                if (e & 1) res = mul(res, b);
                b = mul(b, b);
            }
            return res;
        }();
        return inverse;
    }
};

#endif  // LIBRARY_STRING_ROLLING_HASH_HPP
