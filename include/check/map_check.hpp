#ifndef POPLAR_TRIE_MAP_CHECK_HPP
#define POPLAR_TRIE_MAP_CHECK_HPP

#include <array>
#include <iostream>
#include <queue>
#include <algorithm>

#include "../poplar/bit_tools.hpp"
#include "../poplar/exception.hpp"

#include "CP_info.hpp"

namespace poplar {

// This class implements an updatable associative array whose keys are strings.
// The data structure is based on a dynamic path-decomposed trie described in the following paper,
// - "Dynamic Path-Decomposed Tries" available at https://arxiv.org/abs/1906.06015.
template <typename Trie, typename NLM>
class map_check {
    static_assert(Trie::trie_type_id == NLM::trie_type_id);

  public:
    using this_type = map_check<Trie, NLM>;
    using trie_type = Trie;
    using value_type = typename NLM::value_type;

    static constexpr auto trie_type_id = Trie::trie_type_id;
    static constexpr uint32_t min_capa_bits = Trie::min_capa_bits;

  public:
    // Generic constructor.
    map_check() = default;

    // Class constructor. Initially allocates the hash table of length
    // 2**capa_bits.
    explicit map_check(uint32_t capa_bits, uint64_t lambda = 32) {
        POPLAR_THROW_IF(!is_power2(lambda), "lambda must be a power of 2.");

        is_ready_ = true;
        lambda_ = lambda;
        hash_trie_ = Trie{capa_bits, 8 + bit_tools::ceil_log2(lambda_)};
        label_store_ = NLM{hash_trie_.capa_bits()};
        codes_.fill(UINT8_MAX);
        codes_[0] = static_cast<uint8_t>(num_codes_++);  // terminator
    }

    // Generic destructor.
    ~map_check() = default;

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
                restore_codes_[num_codes_] = *key.begin;
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

    // 新しい辞書に挿入する際に使用
    // どのノードの部分から比較を開始するのか(node_id)
    // 最初はノードに格納されている文字列の比較から始まるので、何文字目からの比較をするのかを表す(first_match)
    value_type* update_new(char_range key, uint64_t node_id, uint64_t first_match) {
        // std::cout << "--- update_new ---" << std::endl;
        if(!hash_trie_.checK_first_insert()) {
            hash_trie_.set_first_insert(true);
            if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
                //std::cout << "length : " << key.length() << ", " << key.empty() << std::endl;
                // pre_node = node_id;
                return label_store_.insert_new_table(hash_trie_.get_root(), key); // ルートから出ている
            }
            // should not come
            assert(false);
        }
        
        // 2つ目以降のキー追加の際の処理
        // auto [vptr, match] = label_store_.compare_new_ptrs(node_id, key, first_match);
        // if(vptr == nullptr) return const_cast<value_type*>(vptr);
        bool first_check = true;
        std::vector<uint64_t> shelter;
        while(!key.empty()) {
            shelter.push_back(node_id);
            uint64_t match;
            if(first_check) {
                // auto [vptr, match_tmp] = label_store_.compare_new_ptrs(node_id, key, first_match);
                auto [vptr, match_tmp] = label_store_.compare_new_ptrs(node_id, key, 0);
                if(vptr != nullptr) return const_cast<value_type*>(vptr);
                match = match_tmp;
                first_check = false;
            } else {
                auto [vptr, match_tmp] = label_store_.compare_new_ptrs(node_id, key, 0);
                if(vptr != nullptr) return const_cast<value_type*>(vptr);
                match = match_tmp;
            }
            
            key.begin += match;
            while(lambda_ <= match) {
                if (hash_trie_.add_child_new_table(node_id, step_symb)) { // step_symbはuint8_tの最大値，つまり255(つまり，ダミーノード)
                    // expand_if_needed_(node_id); // 再配置の際に使用するので、必要なし
#ifdef POPLAR_EXTRA_STATS
                    ++num_steps_;
#endif
                    // if constexpr (trie_type_id == trie_type_ids::FKHASH_TRIE) {
                    //     assert(node_id == label_store_.size());
                    //     label_store_.append_dummy();
                    // }
                    shelter.push_back(node_id);
                }
                match -= lambda_;
            }

            if (codes_[*key.begin] == UINT8_MAX) {
                // Update table
                restore_codes_[num_codes_] = *key.begin;
                codes_[*key.begin] = static_cast<uint8_t>(num_codes_++);
                POPLAR_THROW_IF(UINT8_MAX == num_codes_, "");
            }

            if (hash_trie_.add_child_new_table(node_id, make_symb_(*key.begin, match))) { // make_symb_は match << 8 | *key.begin
                // 条件を満たすのは，繊維が失敗したとき
                // expand_if_needed_(node_id);
                ++key.begin;
                // ++size_;

                // std::cout << "node_id : " << node_id << std::endl;
                if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
                    return label_store_.insert_new_table(node_id, key);
                }
                // should not come
                assert(false);
            }

            ++key.begin;
        }


        auto vptr = label_store_.compare(node_id, key).first;
        return vptr ? const_cast<value_type*>(vptr) : nullptr;
    }

    // 新しい辞書に登録する
    // 条件によって、辞書に対する登録方法を変更する
    void insert_new_dic(uint64_t node_id, bool flag, std::string& compare_str) {
        // std::cout << "--- insert_new_dic ---" << std::endl;
        if(flag) {
            std::string insert_key = restore_insert_string(node_id);
            compare_str = insert_key;
            // std::cout << "insert : " << insert_key << std::endl;
            int* ptr = update_new(make_char_range(insert_key), hash_trie_.get_root(), 0); // 先頭から比較
            *ptr = 1;
            // std::cout << "insert_key : " << insert_key << std::endl;
        } else {
            auto fs = label_store_.return_string(node_id);
            if(fs == nullptr) return;
            std::string insert_key = restore_insert_string(node_id);
            // std::cout << "insert : " << insert_key << std::endl;
            int* ptr = update_new(make_char_range(insert_key), hash_trie_.get_root(), 1); // 特定の箇所から比較
            *ptr = 1;
            // std::cout << "insert_key : " << insert_key << std::endl;
        }
        // std::cout << "end" << std::endl;
    }

    // centroid_path_を求めて、新しい辞書に格納する
    void require_centroid_path_order_and_insert_dictionaly(std::vector<std::vector<info_fp>>& fp,
                                                            uint64_t node_id,
                                                            const std::vector<uint64_t>& bn,
                                                            const std::vector<bool>& check_bottom,
                                                            std::string& compare_str) {
        if(fp[node_id].size() == 0) { // 一番まで、たどり着いた時の処理
            return;
        }
        bool check_zero = false;
        if(fp[node_id][0].match == 0) check_zero = true;
        // 0分岐とそれ以外の個数を比べる
        if(check_zero) { // 0分岐がある時
            std::sort(fp[node_id].begin()+1, fp[node_id].end(), [] (auto l, auto r) {
                return l.cnt > r.cnt;
            });
            // [0]は0分岐を示している
            // bn[]には、そのノードよりも後ろにあるノードの数を保存している(0分岐を除いて)
            // bnを使用している理由は、累積和の時間を省略するため
            if(fp[node_id][0].cnt > bn[node_id]) {
                // それぞれの子を持ってくる、その中から、cntが多い順に処理する
                // 下のshelterに格納する部分を時間短縮できそう(here)
                std::sort(fp[node_id][0].children.begin(), fp[node_id][0].children.end(), [] (auto l, auto r) {
                    return l.second > r.second;
                });
                for(auto s : fp[node_id][0].children) {
                    require_centroid_path_order_and_insert_dictionaly(fp, s.first, bn, check_bottom, compare_str);
                    insert_new_dic(s.first, check_bottom[s.first], compare_str);
                }

                // 1分岐以上の処理
                for(uint64_t i=1; i < fp[node_id].size(); i++) {
                    std::sort(fp[node_id][i].children.begin(), fp[node_id][i].children.end(), [] (auto l, auto r) {
                        return l.second > r.second;
                    });
                    for(auto s : fp[node_id][i].children) {
                        require_centroid_path_order_and_insert_dictionaly(fp, s.first, bn, check_bottom, compare_str);
                        insert_new_dic(s.first, check_bottom[s.first], compare_str);
                    }
                }
            } else { // 最後に0分岐を処理する
                // 1分岐以上の処理
                // std::vector<std::pair<uint64_t, uint64_t>> shelter;
                for(uint64_t i=1; i < fp[node_id].size(); i++) {
                    std::sort(fp[node_id][i].children.begin(), fp[node_id][i].children.end(), [] (auto l, auto r) {
                        return l.second > r.second;
                    });
                    for(auto s : fp[node_id][i].children) {
                        require_centroid_path_order_and_insert_dictionaly(fp, s.first, bn, check_bottom, compare_str);
                        insert_new_dic(s.first, check_bottom[s.first], compare_str);
                    }
                }

                // 0分岐の処理
                std::sort(fp[node_id][0].children.begin(), fp[node_id][0].children.end(), [] (auto l, auto r) {
                    return l.second > r.second;
                });
                for(auto s : fp[node_id][0].children) {
                    require_centroid_path_order_and_insert_dictionaly(fp, s.first, bn, check_bottom, compare_str);
                    insert_new_dic(s.first, check_bottom[s.first], compare_str);
                }
            }
        } else { // 0分岐がないとき
            std::sort(fp[node_id].begin(), fp[node_id].end(), [] (auto l, auto r) {
                return l.cnt > r.cnt;
            });
            for(uint64_t i=0; i < fp[node_id].size(); i++) {
                std::sort(fp[node_id][i].children.begin(), fp[node_id][i].children.end(), [] (auto l, auto r) {
                    return l.second > r.second;
                });
                for(auto s : fp[node_id][i].children) {
                    require_centroid_path_order_and_insert_dictionaly(fp, s.first, bn, check_bottom, compare_str);
                    insert_new_dic(s.first, check_bottom[s.first], compare_str);
                }
            }
        }
        if(node_id == hash_trie_.get_root()) insert_new_dic(node_id, check_bottom[node_id], compare_str);
    }

    // 特定のノードから、get_root()までの文字列を復元する
    std::string restore_insert_string(uint64_t node_id) {
        std::string insert_string = ""; // ここに文字列を格納して、新しい辞書に挿入する

        // とりあえず、文字列を復元する
        auto fs = label_store_.return_string(node_id);
        if(fs == nullptr) return insert_string;
        for(uint64_t i=0;; i++) {
            if(fs[i] == 0x00) break;
            insert_string += fs[i];
        }

        // get_root()まで、文字列を復元する
        while(node_id != hash_trie_.get_root()) {
            auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id); // 親ノードとsymbを取得
            auto [c, match] = restore_symb_(symb); // symbから、遷移に失敗した箇所とlabelを取得する

            insert_string = c + insert_string;

            uint64_t dummy_step = 0; // ダミーノードの数を数える
            while(1) {
                fs = label_store_.return_string(parent);
                if(fs == nullptr) {
                    dummy_step++;
                    auto [tmp1, tmp2] = hash_trie_.get_parent_and_symb(parent);
                    parent = tmp1;
                } else {
                    break;
                }
            }

            match += dummy_step * lambda_; // スキップした回数分足してあげる

            if(match != 0) {
                fs = label_store_.return_string(parent);
                std::string tmp_str = "";
                for(uint64_t j=0; j < match; j++) {
                    tmp_str += fs[j];
                }
                insert_string = tmp_str + insert_string;
            }

            node_id = parent;
        }
        return insert_string;
    }

    // データをリセットするための関数
    void data_reset() {
        hash_trie_.reset_data_();
        label_store_.reset_data_();
    }

    // トポロジカルソートを呼び出す
    void call_topo() {
        std::cout << "--- call_topo ---" << std::endl;
        
        auto [fork_info, blanch_num, check_bottom] = hash_trie_.return_partial_CP_info(restore_codes_); // 
        std::cout << "size : " << fork_info.size() << std::endl;

        // fork_infoの情報を元に、CP順を求め、新しい辞書に格納していく
        data_reset();
        hash_trie_.expand_tmp_table();
        label_store_.expand_tmp_ptrs();
        std::string compare_str = "";
        // std::cout << "aaa" << std::endl;
        require_centroid_path_order_and_insert_dictionaly(fork_info, hash_trie_.get_root(), blanch_num, check_bottom, compare_str);
        // std::cout << "bbb" << std::endl;
        hash_trie_.move_table();    // 辞書の移動
        label_store_.move_ptrs();
        hash_trie_.set_first_insert(false);

        // 特定の深さからの線形探索回数を調べる
        // hash_trie_.reset_cnt_compare();
        // std::queue<std::pair<uint64_t, uint64_t>> que;
        // que.push(std::pair{hash_trie_.get_root(), 0});
        // while(!que.empty()) {
        //     auto [node_id, deep] = que.front();
        //     que.pop();
        //     if(deep <= 1) {
        //         for(auto node : fork_info[node_id]) {
        //             // std::cout << node.match << ", " << node.cnt << ", " << node.children.size() << std::endl;
        //             for(auto [next_id, leaf_cnt] : node.children) {
        //                 que.push(std::pair{next_id, deep+1});
        //             }
        //         }
        //     } else {
        //         auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id);
        //         auto next_id = hash_trie_.find_child(parent, symb);
        //     }
        // }
        // hash_trie_.show_cnt_compare();
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
        show_stat(os, indent, "name", "map_check");
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
    }

    map_check(const map_check&) = delete;
    map_check& operator=(const map_check&) = delete;

    map_check(map_check&&) noexcept = default;
    map_check& operator=(map_check&&) noexcept = default;

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
};

}  // namespace poplar

#endif  // POPLAR_TRIE_MAP_CHECK_HPP
