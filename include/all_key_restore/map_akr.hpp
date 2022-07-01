#ifndef POPLAR_TRIE_MAP_AKR_HPP
#define POPLAR_TRIE_MAP_AKR_HPP

#include <array>
#include <iostream>

#include "../poplar/bit_tools.hpp"
#include "../poplar/exception.hpp"

#include <deque>
#include <map>
#include <algorithm>

namespace poplar {

// static std::map<int, int> cnt_hash;
// static uint64_t all_cnt = 0;

// This class implements an updatable associative array whose keys are strings.
// The data structure is based on a dynamic path-decomposed trie described in the following paper,
// - "Dynamic Path-Decomposed Tries" available at https://arxiv.org/abs/1906.06015.
template <typename Trie, typename NLM>
class map_akr {
    static_assert(Trie::trie_type_id == NLM::trie_type_id);

  public:
    using this_type = map_akr<Trie, NLM>;
    using trie_type = Trie;
    using value_type = typename NLM::value_type;

    static constexpr auto trie_type_id = Trie::trie_type_id;
    static constexpr uint32_t min_capa_bits = Trie::min_capa_bits;

  public:
    // Generic constructor.
    map_akr() = default;

    // Class constructor. Initially allocates the hash table of length
    // 2**capa_bits.
    explicit map_akr(uint32_t capa_bits, uint64_t lambda = 32) {
        POPLAR_THROW_IF(!is_power2(lambda), "lambda must be a power of 2.");

        is_ready_ = true;
        lambda_ = lambda;
        hash_trie_ = Trie{capa_bits, 8 + bit_tools::ceil_log2(lambda_)};
        label_store_ = NLM{hash_trie_.capa_bits()};
        codes_.fill(UINT8_MAX);
        restore_codes_.fill(UINT8_MAX);
        restore_codes_[0] = static_cast<uint8_t>(num_codes_);
        codes_[0] = static_cast<uint8_t>(num_codes_++);  // terminator
    }

    // Generic destructor.
    ~map_akr() = default;

    // Searches the given key and returns the value pointer if registered;
    // otherwise returns nullptr.
    const value_type* find(const std::string& key) const {
        return find(make_char_range(key));
    }
    const value_type* find(char_range key) const {
        POPLAR_THROW_IF(key.empty(), "key must be a non-empty string.");
        POPLAR_THROW_IF(*(key.end - 1) != '\0', "The last character of key must be the null terminator.");

        if (!is_ready_ or hash_trie_.size() == 0) {
            return nullptr;
        }

        auto node_id = hash_trie_.get_root();
        all_cnt = 0;
        int cnt = 0;

        while (!key.empty()) {
            auto [vptr, match] = label_store_.compare(node_id, key);
            if (vptr != nullptr) {
                cnt_hash[cnt] += 1;
                return vptr;
            }

            key.begin += match;

            while (lambda_ <= match) {
                node_id = hash_trie_.find_child(node_id, step_symb);
                if (node_id == nil_id) {
                    return nullptr;
                }
                match -= lambda_;
                cnt += 1;
            }

            if (codes_[*key.begin] == UINT8_MAX) {
                // Detecting an useless character
                return nullptr;
            }

            node_id = hash_trie_.find_child(node_id, make_symb_(*key.begin, match));
            if (node_id == nil_id) {
                return nullptr;
            }
            cnt += 1;

            ++key.begin;
        }

        cnt_hash[cnt] += 1;
        return label_store_.compare(node_id, key).first;
    }

    // Inserts the given key and returns the value pointer.
    value_type* update(const std::string& key) {
        return update(make_char_range(key));
    }
    value_type* update(char_range key) {
        POPLAR_THROW_IF(key.empty(), "key must be a non-empty string.");
        POPLAR_THROW_IF(*(key.end - 1) != '\0', "The last character of key must be the null terminator.");

        char_range tmp_key = key;
        if (hash_trie_.size() == 0) {
            if (!is_ready_) {
                *this = this_type{0};
            }
            // The first insertion
            ++size_;
            hash_trie_.add_root();

            if constexpr (trie_type_id == trie_type_ids::FKHASH_TRIE) {
                // assert(hash_trie_.get_root() == label_store_.size());
                return label_store_.append(key);
            }
            if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
                return label_store_.insert(hash_trie_.get_root(), key);
            }
            // should not come
            assert(false);
        }

        auto node_id = hash_trie_.get_root();

        while (!key.empty()) {
            auto [vptr, match] = label_store_.compare(node_id, key);
            if (vptr != nullptr) {
                return const_cast<value_type*>(vptr);
            }

            key.begin += match;

            while (lambda_ <= match) {
                if (hash_trie_.add_child(node_id, step_symb)) {
                    // expand_if_needed_(node_id);
                    if(is_need_expand()) {
                        // std::cout << "key : " << tmp_key.begin << std::endl;
                        return dynamic_replacement(tmp_key);
                    }
#ifdef POPLAR_EXTRA_STATS
                    ++num_steps_;
#endif
                    if constexpr (trie_type_id == trie_type_ids::FKHASH_TRIE) {
                        assert(node_id == label_store_.size());
                        label_store_.append_dummy();
                    }
                }
                match -= lambda_;
            }

            if (codes_[*key.begin] == UINT8_MAX) {
                // Update table
                restore_codes_[num_codes_] = *key.begin;
                codes_[*key.begin] = static_cast<uint8_t>(num_codes_++);
                POPLAR_THROW_IF(UINT8_MAX == num_codes_, "");
            }

            if (hash_trie_.add_child(node_id, make_symb_(*key.begin, match))) {
                // expand_if_needed_(node_id);
                if(is_need_expand()) {
                    // std::cout << "key : " << tmp_key.begin << std::endl;
                    return dynamic_replacement(tmp_key);
                }
                ++key.begin;
                ++size_;

                if constexpr (trie_type_id == trie_type_ids::FKHASH_TRIE) {
                    assert(node_id == label_store_.size());
                    return label_store_.append(key);
                }
                if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
                    return label_store_.insert(node_id, key);
                }
                // should not come
                assert(false);
            }

            ++key.begin;
        }

        auto vptr = label_store_.compare(node_id, key).first;
        return vptr ? const_cast<value_type*>(vptr) : nullptr;
    }

    std::vector<std::string> all_key_restore() {
        std::vector<std::string> restore_keys;                      // 復元したキーを保存
        std::vector<std::string> part_keys(hash_trie_.capa_size()); // 部分的な文字列を保存

        uint64_t first_node = hash_trie_.get_root();
        std::string restore_str = "";

        // 先頭ノードに対する処理
        auto fs = label_store_.return_string_pointer(first_node);
        for(uint64_t i=0;;i++) {
            if(fs[i] == 0x00) break; // 判定処理としてint型で4つ分使用して1を格納しているので、この判定
            restore_str += fs[i];
        }
        restore_keys.emplace_back(restore_str);

        // 先頭ノード以外のノードに対する処理
        for(uint64_t i=first_node+1; i < hash_trie_.capa_size(); i++) {
            // ハッシュテーブルが使用されている → 終端となる文字列が存在(ダミーノードを除いて)
            if(!hash_trie_.is_use_table(i)) continue;

            uint64_t node_id = i;
            restore_str.clear();
            std::deque<uint64_t> save_route;
            auto str = label_store_.return_string_pointer(node_id); // ダミーノードで遷移した際には、中身が何もないのでcontinue
            if(str == nullptr) continue;
            for(uint64_t j=0;;j++) {
                if(str[j] == 0x00) break;   // 最後にint型の1が格納されていることを利用する
                restore_str += str[j];
            }
            while(node_id != first_node) {
                auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id); // 親ノードとsymbを取得
                auto [c, match] = restore_symb_(symb); // symbから、遷移に失敗した箇所とlabelを取得する
                restore_str = c + restore_str; // 遷移文字を追加
                
                uint64_t dummy_step = 0; // ダミーノード数をカウント
                while(1) { // 親がダミーノードの場合の処理
                    auto str1 = label_store_.return_string_pointer(parent);
                    if(str1 == nullptr) {
                        dummy_step++;
                        auto [parent_tmp, tmp2] = hash_trie_.get_parent_and_symb(parent);
                        parent = parent_tmp;
                    } else {
                        break;
                    }
                }
                match += dummy_step * lambda_; // スキップした回数分足してあげる

                // 親ノードの文字列からmatchまでを抜き出し、追加する
                if(match != 0) { // matchがゼロの時は、文字列比較の一文字目で遷移しているので
                    auto str1 = label_store_.return_string_pointer(parent);
                    std::string str2 = "";
                    for(uint64_t j=0; j < match; j++) {
                        str2 += str1[j];
                    }
                    restore_str = str2 + restore_str;
                    part_keys[node_id] = str2 + c;
                } else {
                    part_keys[node_id] = c;
                }
                save_route.push_back(node_id);
                node_id = parent;
                if(part_keys[node_id].size() != 0) { // 一度、通ったことのあるノードだった場合
                    restore_str = part_keys[node_id] + restore_str;
                    break;
                }
            }
            // part_keysの格納
            while(!save_route.empty()) {
                uint64_t node = save_route.back();
                save_route.pop_back();
                part_keys[node] = part_keys[node_id] + part_keys[node];
                node_id = node;
            }
            restore_keys.emplace_back(restore_str);
        }

        return restore_keys;
    }

    // 途中の文字列を保存しないバージョン
    std::vector<std::string> all_key_restore_simple() {
        std::vector<std::string> restore_keys;  // 復元したキーを全て保存

        uint64_t first_node = hash_trie_.get_root();
        std::string restore_str = "";

        // 先頭ノードに対する処理
        auto fs = label_store_.return_string_pointer(first_node);
        for(uint64_t i=0;;i++) {
            if(fs[i] == 0x00) break; // 判定処理としてint型で4つ分使用して1を格納しているので、この判定
            restore_str += fs[i];
        }
        restore_keys.emplace_back(restore_str);

        for(uint64_t i=first_node+1; i < hash_trie_.capa_size(); i++) {
            if(!hash_trie_.is_use_table(i)) continue;
            std::string restore_str = "";
            uint64_t node_id = i;
            auto str = label_store_.return_string_pointer(node_id);
            if(str == nullptr) continue; // ダミーノード判定
            for(uint64_t j=0;;j++) {
                if(str[j] == 0x00) break;
                restore_str += str[j];
            }

            while(node_id != first_node) {
                auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id); // 親ノードとsymbを取得
                auto [c, match] = restore_symb_(symb); // symbから、遷移に失敗した箇所とlabelを取得する
                restore_str = c + restore_str; // 遷移文字を追加

                uint64_t dummy_step = 0; // ダミーノード数をカウント
                while(1) { // 親がダミーノードの場合の処理
                    auto str1 = label_store_.return_string_pointer(parent);
                    if(str1 == nullptr) {
                        dummy_step++;
                        auto [parent_tmp, tmp2] = hash_trie_.get_parent_and_symb(parent);
                        parent = parent_tmp;
                    } else {
                        break;
                    }
                }
                match += dummy_step * lambda_; // スキップした回数分足してあげる

                // 親ノードの文字列からmatchまでを抜き出し、追加する
                if(match != 0) { // matchがゼロの時は、文字列比較の一文字目で遷移しているので
                    auto str1 = label_store_.return_string_pointer(parent);
                    std::string str2 = "";
                    for(uint64_t j=0; j < match; j++) {
                        str2 += str1[j];
                    }
                    restore_str = str2 + restore_str;
                }

                node_id = parent;
            }
            restore_keys.emplace_back(restore_str);
        }

        return restore_keys;
    }

    // CPD順にソートし、新しい辞書に追加
    template <class Map, class It>
    void insert_by_centroid_path_order(Map& new_map, It begin, It end, uint64_t depth) {
        assert(end-begin > 0);
        if (end-begin == 1) {
            int* ptr = new_map.update(*begin);
            *ptr = 1;
            return;
        }
        std::vector<std::tuple<It,It,char>> ranges;
        auto from = begin;
        auto to = begin;
        if (from->length() == depth) {
            ranges.emplace_back(from, ++to, '\0');
            from = to;
        }
        while (from != end) {
            assert(from->length() > depth);
            char c = (*from)[depth];
            while (to != end and (*to)[depth] == c) {
                ++to;
            }
            ranges.emplace_back(from, to, c);
            from = to;
        }
        std::sort(ranges.begin(), ranges.end(), [](auto l, auto r) {
            auto [fl,tl,cl] = l;
            auto [fr,tr,cr] = r;
            return tl-fl > tr-fr;
        });
        for (auto [f,t,c] : ranges) {
            insert_by_centroid_path_order(new_map, f, t, depth+1);
        }
    }

    value_type* dynamic_replacement(char_range& key) {
    // void dynamic_replacement() {
        // std::vector<std::string> all_keys = all_key_restore();  // 文字列の復元(部分文字列保存あり)
        std::vector<std::string> all_keys = all_key_restore_simple(); // 文字列の復元

        // map_akr new_map(hash_trie_.capa_bits());                // 新しい辞書の作成
        map_akr new_map(hash_trie_.capa_bits()+1);
        std::sort(all_keys.begin(), all_keys.end());
        insert_by_centroid_path_order(new_map, all_keys.begin(), all_keys.end(), 0);
        std::swap(*this, new_map);  // 辞書を入れ替える

        return update(key);
    }

    // Gets the number of registered keys.
    uint64_t size() const {
        return size_;
    }
    // Gets the capacity of the hash table.
    uint64_t capa_size() const {
        return hash_trie_.capa_size();
    }
#ifdef POPLAR_EXTRA_STATS
    double rate_steps() const {
        return double(num_steps_) / size_;
    }
    uint64_t num_resize() const {
        return hash_trie_.num_resize();
    }
#endif
    uint64_t alloc_bytes() const {
        uint64_t bytes = 0;
        bytes += hash_trie_.alloc_bytes();
        bytes += label_store_.alloc_bytes();
        bytes += codes_.size();
        return bytes;
    }

    void reset_cnt_hash() {
        std::cout << "--- reset_cnt_hash ---" << std::endl;
        cnt_hash.clear();
        all_cnt = 0;
    }

    void show_stats(std::ostream& os, int n = 0) const {
        auto indent = get_indent(n);
        show_stat(os, indent, "name", "map_akr");
        show_stat(os, indent, "lambda", lambda_);
        show_stat(os, indent, "size", size());
        show_stat(os, indent, "alloc_bytes", alloc_bytes());
#ifdef POPLAR_EXTRA_STATS
        show_stat(os, indent, "rate_steps", rate_steps());
#endif
        show_member(os, indent, "hash_trie_");
        hash_trie_.show_stats(os, n + 1);
        show_member(os, indent, "label_store_");
        label_store_.show_stats(os, n + 1);
    }

    void show_cnt_hash() {
        std::cout << "--- show_cnt_hash ---" << std::endl;
        std::cout << "all_cnt : " << all_cnt << std::endl;
        uint64_t all = 0;
        uint64_t sum = 0;
        for(auto p : cnt_hash) {
            // std::cout << p.first << " : " << p.second << std::endl;
            // std::cout << p.second << std::endl;
            std::cout << p.first << std::endl;
            all += p.second;
        }
        std::cout << "----" << std::endl;
        for(auto p : cnt_hash) {
            std::cout << p.second << std::endl;
            sum += p.first * p.second;
        }

        std::cout << "all : " << all << std::endl;
        std::cout << "ave_hash : " << (long double)(sum) / (long double)(11414967) << std::endl; // enwiki
        // std::cout << "ave_hash : " << (long double)(sum) / (long double)(52616588) << std::endl; // DS5
        // std::cout << "ave_hash : " << (long double)(sum) / (long double)(10154743) << std::endl; // AOL
        // std::cout << "ave_hash : " << (long double)(sum) / (long double)(7609320) << std::endl; // GeoNames
        // std::cout << "ave_hash : " << (long double)(sum) / (long double)(1382908) << std::endl; // in-2004
        // std::cout << "ave_hash : " << (long double)(sum) / (long double)(39459925) << std::endl; // uk-2005
        // std::cout << "ave_hash : " << (long double)(sum) / (long double)(118142155) << std::endl; // webbase-2001
    }

    map_akr(const map_akr&) = delete;
    map_akr& operator=(const map_akr&) = delete;

    map_akr(map_akr&&) noexcept = default;
    map_akr& operator=(map_akr&&) noexcept = default;

  private:
    static constexpr uint64_t nil_id = Trie::nil_id;
    static constexpr uint64_t step_symb = UINT8_MAX;  // (UINT8_MAX, 0)

    bool is_ready_ = false;
    uint64_t lambda_ = 32;

    Trie hash_trie_;
    NLM label_store_;
    std::array<uint8_t, 256> codes_ = {};
    std::array<uint8_t, 256> restore_codes_ = {};
    uint32_t num_codes_ = 0;
    uint64_t size_ = 0;
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_steps_ = 0;
#endif

    uint64_t make_symb_(uint8_t c, uint64_t match) const {
        assert(codes_[c] != UINT8_MAX);
        return static_cast<uint64_t>(codes_[c]) | (match << 8);
    }

    // labelからcとmatchを復元する
    std::pair<char, uint64_t> restore_symb_(uint64_t label) const  {
        return std::pair{char(restore_codes_[label % 256]), label/256};
    }

    void expand_if_needed_(uint64_t& node_id) {
        if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
            if (!hash_trie_.needs_to_expand()) {
                return;
            }
            auto node_map = hash_trie_.expand();
            node_id = node_map[node_id];
            label_store_.expand(node_map);
        }
    }

    // ハッシュテーブルの拡張が必要かを調べるための関数
    bool is_need_expand() {
        if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
            if(hash_trie_.needs_to_expand()) return true;
        }
        return false;
    }
};

}  // namespace poplar

#endif  // POPLAR_TRIE_MAP_AKR_HPP
