#ifndef LIBRARY_STRING_TRIE_HPP
#define LIBRARY_STRING_TRIE_HPP

#line 1 "src/string/trie.hpp"


#include <array>
#include <cassert>
#include <string>
#include <vector>

/// 文字を [margin, margin + alphabet_size) として扱う Trie。重複する文字列も保持する。
template <int alphabet_size = 26, char margin = 'a'>
struct trie {
    static_assert(alphabet_size > 0);

  private:
    struct node {
        std::array<int, alphabet_size> next;
        int pass_count = 0;
        int terminal_count = 0;

        node() { next.fill(-1); }
    };

  public:
    /// 空のTrieを作る。
    trie() : _nodes(1) {}

    /// 頂点配列の容量を node_capacity まで予約する。
    void reserve(int node_capacity) {
        assert(node_capacity >= 0);
        _nodes.reserve(node_capacity);
    }

    /// 文字列 s を 1 個追加する。計算量 O(|s|)。
    void insert(const std::string& s) {
        int v = 0;
        _nodes[v].pass_count++;
        for (char c : s) {
            int x = char_id(c);
            if (_nodes[v].next[x] == -1) {
                _nodes[v].next[x] = int(_nodes.size());
                _nodes.emplace_back();
            }
            v = _nodes[v].next[x];
            _nodes[v].pass_count++;
        }
        if (_nodes[v].terminal_count++ == 0) _distinct_size++;
        _size++;
    }

    /// 文字列 s を 1 個削除する。s が存在すること。計算量 O(|s|)。
    void erase(const std::string& s) {
        assert(count(s) > 0);
        int v = 0;
        _nodes[v].pass_count--;
        for (char c : s) {
            v = _nodes[v].next[char_id(c)];
            _nodes[v].pass_count--;
        }
        if (--_nodes[v].terminal_count == 0) _distinct_size--;
        _size--;
    }

    /// 文字列 s の登録個数を返す。計算量 O(|s|)。
    int count(const std::string& s) const {
        int v = find_node(s);
        return v == -1 ? 0 : _nodes[v].terminal_count;
    }

    /// 文字列 s が 1 個以上登録されていれば true。計算量 O(|s|)。
    bool contains(const std::string& s) const { return count(s) > 0; }

    /// prefix から始まる登録文字列の個数を返す。重複を含む。計算量 O(|prefix|)。
    int count_starts_with(const std::string& prefix) const {
        int v = find_node(prefix);
        return v == -1 ? 0 : _nodes[v].pass_count;
    }

    /// s のprefixとして登録されている文字列の個数を返す。重複を含む。計算量 O(|s|)。
    int count_prefixes_of(const std::string& s) const {
        int v = 0;
        int res = _nodes[v].terminal_count;
        for (char c : s) {
            int x = char_id(c);
            v = _nodes[v].next[x];
            if (v == -1) break;
            res += _nodes[v].terminal_count;
        }
        return res;
    }

    /// s のprefixとして登録されている最長文字列の長さを返す。存在しなければ -1。
    int longest_prefix_length(const std::string& s) const {
        int v = 0;
        int res = _nodes[v].terminal_count > 0 ? 0 : -1;
        for (int i = 0; i < int(s.size()); i++) {
            int x = char_id(s[i]);
            v = _nodes[v].next[x];
            if (v == -1) break;
            if (_nodes[v].terminal_count > 0) res = i + 1;
        }
        return res;
    }

    /// 登録されている文字列の総数を返す。重複を含む。
    int size() const { return _size; }

    /// 登録されている異なる文字列の個数を返す。
    int distinct_size() const { return _distinct_size; }

    /// 文字列が 1 個も登録されていなければ true。
    bool empty() const { return _size == 0; }

    /// Trie が確保している頂点数を返す。
    int node_count() const { return int(_nodes.size()); }

  private:
    static int char_id(char c) {
        int x = int(c) - int(margin);
        assert(0 <= x && x < alphabet_size);
        return x;
    }

    int find_node(const std::string& s) const {
        int v = 0;
        for (char c : s) {
            int x = char_id(c);
            v = _nodes[v].next[x];
            if (v == -1) return -1;
        }
        return v;
    }

    std::vector<node> _nodes;
    int _size = 0;
    int _distinct_size = 0;
};

#endif  // LIBRARY_STRING_TRIE_HPP
