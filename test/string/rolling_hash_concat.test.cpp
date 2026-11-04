#include <algorithm>
#include <cassert>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../../src/string/rolling_hash.hpp"

using u64 = unsigned long long;

u64 naive_hash(const std::vector<u64>& values, int l, int r) {
    u64 hash = 0, base = rolling_hash::get_base();
    for (int i = l; i < r; i++) {
        hash = (static_cast<unsigned __int128>(hash) * base + values[i]) % rolling_hash::MOD;
    }
    return hash;
}

std::vector<u64> bytes(const std::string& s) {
    std::vector<u64> values(s.size());
    for (int i = 0; i < int(s.size()); i++) values[i] = static_cast<unsigned char>(s[i]);
    return values;
}

void check_all(const rolling_hash& hash, const std::vector<u64>& values) {
    int n = int(values.size());
    assert(hash.size() == n);
    for (int l = 0; l <= n; l++) {
        for (int r = l; r <= n; r++) assert(hash.get(l, r) == naive_hash(values, l, r));
    }
}

void check_samples(const rolling_hash& hash, const std::string& s, std::mt19937& rng) {
    int n = int(s.size());
    auto values = bytes(s);
    assert(hash.size() == n);
    assert(hash.get(0, n) == naive_hash(values, 0, n));
    assert(hash.get(0, 0) == 0 && hash.get(n, n) == 0);
    for (int q = 0; q < 20; q++) {
        int l = int(rng() % (n + 1)), r = int(rng() % (n + 1));
        if (l > r) std::swap(l, r);
        assert(hash.get(l, r) == naive_hash(values, l, r));
        assert(rolling_hash::connect(hash.get(0, l), hash.get(l, n), n - l) == hash.get(0, n));
    }
}

int main() {
    rolling_hash empty;
    check_all(empty, {});
    empty.concat_left(std::string());
    empty.concat_right(std::string());
    check_all(empty, {});

    rolling_hash mixed(std::string("de"));
    mixed.concat_left(std::string("bc"));
    mixed.concat_left('a');
    mixed.concat_right(std::string("fg"));
    check_all(mixed, bytes("abcdefg"));
    mixed.build(std::string("xyz"));
    check_all(mixed, bytes("xyz"));
    mixed.build(std::string());
    mixed.concat_left('a');
    check_all(mixed, bytes("a"));

    std::vector<u64> integers{0, 255, 256, rolling_hash::MOD - 2, 17};
    rolling_hash integer_hash(integers);
    integer_hash.concat_left('a');
    integers.insert(integers.begin(), 'a');
    integer_hash.concat_right('\0');
    integers.push_back(0);
    check_all(integer_hash, integers);

    std::mt19937 rng(459);
    for (int tc = 0; tc < 20; tc++) {
        rolling_hash hash;
        std::string s;
        for (int q = 0; q < 200; q++) {
            int type = int(rng() % 7);
            char c = static_cast<char>(rng() % 256);
            std::string part(int(rng() % 6), '\0');
            for (char& value : part) value = static_cast<char>(rng() % 256);
            if (type == 0) {
                hash.concat_left(c);
                s.insert(s.begin(), c);
            } else if (type == 1) {
                hash.concat_right(c);
                s.push_back(c);
            } else if (type == 2) {
                hash.concat_left(part);
                s = part + s;
            } else if (type == 3) {
                hash.concat_right(part);
                s += part;
            } else if (type == 4) {
                hash.build(part);
                s = part;
            } else if (type == 5) {
                auto copy = hash;
                copy.concat_left(c);
                check_samples(copy, c + s, rng);
                hash = copy;
                s = c + s;
            } else {
                auto copy = hash;
                rolling_hash moved(std::move(copy));
                check_samples(moved, s, rng);
                copy.build(part);
                check_samples(copy, part, rng);
            }
            check_samples(hash, s, rng);
        }
    }

    const int n = 65536;
    const std::string one = "a";
    std::vector<u64> repeated(n, 'a');
    u64 expected = naive_hash(repeated, 0, n);
    rolling_hash left_chars, left_strings, right_chars, right_strings;
    for (int i = 0; i < n; i++) {
        left_chars.concat_left('a');
        left_strings.concat_left(one);
        right_chars.concat_right('a');
        right_strings.concat_right(one);
    }
    for (const auto* hash : {&left_chars, &left_strings, &right_chars, &right_strings}) {
        assert(hash->size() == n);
        assert(hash->get(0, n) == expected);
        assert(hash->get(321, 654) == naive_hash(repeated, 321, 654));
    }
}
