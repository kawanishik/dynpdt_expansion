#ifndef POPLAR_TRIE_MAP_TABLE_HPP
#define POPLAR_TRIE_MAP_TABLE_HPP

#include <array>
#include <iostream>

#include "../poplar/bit_tools.hpp"
#include "../poplar/exception.hpp"

static int cnt_expand = 0;

namespace poplar {

// This class implements an updatable associative array whose keys are strings.
// The data structure is based on a dynamic path-decomposed trie described in the following paper,
// - "Dynamic Path-Decomposed Tries" available at https://arxiv.org/abs/1906.06015.
template <typename Trie, typename NLM>
class map_table {
    static_assert(Trie::trie_type_id == NLM::trie_type_id);

  public:
    using this_type = map_table<Trie, NLM>;
    using trie_type = Trie;
    using value_type = typename NLM::value_type;

    static constexpr auto trie_type_id = Trie::trie_type_id;
    static constexpr uint32_t min_capa_bits = Trie::min_capa_bits;

  public:
    // Generic constructor.
    map_table() = default;

    // Class constructor. Initially allocates the hash table of length
    // 2**capa_bits.
    explicit map_table(uint32_t capa_bits, uint64_t lambda = 32) {
        POPLAR_THROW_IF(!is_power2(lambda), "lambda must be a power of 2.");

        is_ready_ = true;
        lambda_ = lambda;
        hash_trie_ = Trie{capa_bits, 8 + bit_tools::ceil_log2(lambda_)};
        label_store_ = NLM{hash_trie_.capa_bits()};
        codes_.fill(UINT8_MAX);
        codes_[0] = static_cast<uint8_t>(num_codes_++);  // terminator
    }

    // Generic destructor.
    ~map_table() = default;

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
                if(hash_trie_.is_node_use_table(node_id)) node_id = hash_trie_.transiton_not_use_hash_table(node_id, step_symb);
                else node_id = hash_trie_.find_child(node_id, step_symb);
                
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
            if(hash_trie_.is_node_use_table(node_id)) node_id = hash_trie_.transiton_not_use_hash_table(node_id, make_symb_(*key.begin, match));
            else node_id = hash_trie_.find_child(node_id, make_symb_(*key.begin, match));

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
                    expand_if_needed_(node_id);
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
                codes_[*key.begin] = static_cast<uint8_t>(num_codes_++);
                POPLAR_THROW_IF(UINT8_MAX == num_codes_, "");
            }
            
            if (hash_trie_.add_child(node_id, make_symb_(*key.begin, match))) {
                expand_if_needed_(node_id);
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

    // void write_file(std::vector<std::string>& tmp_keys) {
    void write_file(std::map<uint64_t, uint64_t>& mp) {
        std::ofstream of;
        // std::string filename = "../../../dataset/in-2004_restore.txt";
        // of.open(filename, std::ios::out);
        // for(auto b : tmp_keys) {
        //     of << b << std::endl;
        // }

        std::string filename = "../../../memo/";
        filename += std::to_string(cnt_expand);
        filename += ".txt";
        std::cout << filename << std::endl;

        of.open(filename, std::ios::out);
        for(auto m : mp) {
            of << m.first << ":" << m.second << std::endl;
        }

        of.close();
    }

    // それぞれのノードに対して、子ノードの数をカウントする
    // 葉ノードは子ノードが存在しないので、個数は0となる
    // 
    void compute_number_of_children_per_node(std::vector<uint64_t>& cnt_number_of_children_per_node) {
        uint64_t table_size = hash_trie_.capa_size();
        cnt_number_of_children_per_node.resize(table_size, 0);
        // カウントする際には、ダミーノードもひとつのノードとしてカウントする
        for(uint64_t i=hash_trie_.get_root()+1; i < table_size; i++) {
            if(hash_trie_.is_use_table(i)) {
                // if(label_store_.return_string_pointer(i) == nullptr) continue;
                auto [p, label] = hash_trie_.get_parent_and_symb(i);
                // if(label == 255) continue;  // ダミーノード
                // while(label_store_.return_string_pointer(p) == nullptr) {
                //     auto [tmp_parent, tmp_label] = hash_trie_.get_parent_and_symb(p);
                //     p = tmp_parent;
                // }
                cnt_number_of_children_per_node[p] += 1;
            }
        }
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
        show_stat(os, indent, "name", "map_table");
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

    void check_not_use_hash() {
        hash_trie_.check_not_use_hash();
    }

    map_table(const map_table&) = delete;
    map_table& operator=(const map_table&) = delete;

    map_table(map_table&&) noexcept = default;
    map_table& operator=(map_table&&) noexcept = default;

  private:
    static constexpr uint64_t nil_id = Trie::nil_id;
    static constexpr uint64_t step_symb = UINT8_MAX;  // (UINT8_MAX, 0)

    bool is_ready_ = false;
    uint64_t lambda_ = 32;

    Trie hash_trie_;
    NLM label_store_;
    std::array<uint8_t, 256> codes_ = {};
    uint32_t num_codes_ = 0;
    uint64_t size_ = 0;
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_steps_ = 0;
#endif

    uint64_t make_symb_(uint8_t c, uint64_t match) const {
        assert(codes_[c] != UINT8_MAX);
        return static_cast<uint64_t>(codes_[c]) | (match << 8);
    }

    void expand_if_needed_(uint64_t& node_id) {
        if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
            if (!hash_trie_.needs_to_expand()) {
                return;
            }

            // 子ノードの数を計算する
            cnt_expand++;
            std::vector<uint64_t> cnt_number_of_children_per_node;
            compute_number_of_children_per_node(cnt_number_of_children_per_node);
            std::map<uint64_t, uint64_t> number_of_children_per_node_map;
            uint64_t cnt = 0;   // 子を持っているノード数の合計
            for(uint64_t i=0; i < cnt_number_of_children_per_node.size(); i++) {
                if(cnt_number_of_children_per_node[i] == 0) continue;
                number_of_children_per_node_map[cnt_number_of_children_per_node[i]] += 1;
                cnt++;
            }

            // 子の数の最低ラインを数える
            uint64_t limit = cnt / 1000;
            uint64_t min_number_of_children = nil_id;
            cnt = 0;
            std::map<uint64_t, uint64_t>::reverse_iterator it = number_of_children_per_node_map.rbegin();
            for(; it != number_of_children_per_node_map.rend(); it++) {
                cnt += it->second;
                if(cnt > limit) break;
                min_number_of_children = it->first;
            }

            // write_file(mp);
            // std::cout << "min_number_of_children : " << min_number_of_children << std::endl;
            // for(int i=0; i < cnt_number_of_children_per_node.size(); i++) {
            //     if(cnt_number_of_children_per_node[i] != 0) std::cout << i << " : " << cnt_number_of_children_per_node[i] << std::endl;
            // }

            // 通常のハッシュテーブルの拡張
            auto node_map = hash_trie_.expand(cnt_number_of_children_per_node, min_number_of_children, lambda_);
            node_id = node_map[node_id];
            label_store_.expand(node_map);
        }
    }
};

}  // namespace poplar

#endif  // POPLAR_TRIE_MAP_TABLE_HPP
