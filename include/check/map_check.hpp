#ifndef POPLAR_TRIE_MAP_CHECK_HPP
#define POPLAR_TRIE_MAP_CHECK_HPP

#include <array>
#include <iostream>
#include <queue>
#include <deque>
#include <stack>
#include <algorithm>
#include <map>

#include "../poplar/bit_tools.hpp"
#include "../poplar/exception.hpp"

#include "CP_info.hpp"

// static std::vector<std::string> CPD_keys;

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
        restore_codes_.fill(UINT8_MAX);
        restore_codes_[0] = static_cast<uint8_t>(num_codes_);
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
                        std::cout << "key : " << tmp_key.begin << std::endl;
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
                    std::cout << "key : " << tmp_key.begin << std::endl;
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
            auto fs = label_store_.return_string_pointer(node_id);
            if(fs == nullptr) return;
            std::string insert_key = restore_insert_string(node_id);
            // std::cout << "insert : " << insert_key << std::endl;
            int* ptr = update_new(make_char_range(insert_key), hash_trie_.get_root(), 1); // 特定の箇所から比較
            *ptr = 1;
            // std::cout << "insert_key : " << insert_key << std::endl;
        }
        // std::cout << "end" << std::endl;
    }

    template<class Map>
    void insert_new_dic_CPD(Map& new_map, uint64_t node_id, uint64_t common_prefix_length, std::string& store_string) {
    // void insert_new_dic_CPD(Map& new_map, uint64_t node_id, uint64_t common_prefix_length, std::string& store_string, char c) {
        // std::cout << "prefix_str : " << prefix_str << std::endl;
        if(node_id == 117464) { // 32246
            std::cout << "^^^^^^^^^^^^" << std::endl;
            std::cout << "particular_part" << std::endl;
            std::cout << "common_prefix_length : " << common_prefix_length << std::endl;
            std::cout << "store_string : " << store_string << std::endl;
            std::cout << "restore_node_string : " << restore_node_string(node_id) << std::endl;
            std::cout << "^^^^^^^^^^^^" << std::endl;
        }
        if(node_id == 176082) { // 84361
            std::cout << "~~~~~~~~~~~~~" << std::endl;
            std::cout << "particular_part" << std::endl;
            std::cout << "common_prefix_length : " << common_prefix_length << std::endl;
            std::cout << "store_string : " << store_string << std::endl;
            std::cout << "restore_node_string : " << restore_node_string(node_id) << std::endl;
            std::cout << "~~~~~~~~~~~~~" << std::endl;
        }
        std::string restore_key = "";
        if(common_prefix_length == 0) {
            restore_key = restore_insert_string(node_id);
            if(restore_key.size() == 0) return;
            // if(restore_key[0] == 'I') {
            //     std::cout << "node_id : " << node_id << std::endl;
            //     std::cout << "restore_key : " << restore_key << std::endl;
            // }
            store_string = restore_key;
            // CPD_keys.emplace_back(restore_key);
            int* ptr = new_map.update(restore_key);
            *ptr = 1;
            return;
        } else {
            // restore_key = store_string.substr(0, common_prefix_length-1) + c + restore_node_string(node_id);
            // restore_key = store_string.substr(0, common_prefix_length) + restore_node_string(node_id);

            // 遷移文字を取得し、文字列を復元する
            auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id); // 親ノードとsymbを取得
            auto [c, match] = restore_symb_(symb);                         // symbから、遷移に失敗した箇所とlabelを取得する

            if(common_prefix_length == 1) restore_key = c + restore_node_string(node_id);
            else restore_key = store_string.substr(0, common_prefix_length-1) + c + restore_node_string(node_id);

            // if(restore_key[0] == 'I') {
            //     std::cout << "node_id : " << node_id << std::endl;
            //     std::cout << "restore_key : " << restore_key << std::endl;
            // }
            // if(node_id == 32246) {
            //     std::cout << "store_string : " << store_string << std::endl;
            //     std::cout << "restore_key  : " << restore_key << std::endl;
            //     std::cout << "node_string  : " << restore_node_string(node_id) << std::endl;
            //     std::cout << "c : " << c << std::endl;
            // }

            if(restore_key.substr(0, 16) == "Cornu_Lunparkeri") {
                std::cout << "node_id : " << node_id << std::endl;
                std::cout << "parent  : " << parent << std::endl;
                std::cout << "c       : " << c << std::endl;
                std::cout << "match   : " << match << std::endl;
                std::cout << restore_key << std::endl;
            }
        }
        // std::string restore_key = restore_insert_string(node_id);
        if(restore_key.size() == 0) return;
        // CPD_keys.emplace_back(restore_key);

        int* ptr = new_map.update(restore_key);
        *ptr = 1;
        // std::cout << "restore_key : " << restore_key << std::endl;
    }

    // 特定のノードに紐づいている文字列を返す
    std::string restore_node_string(uint64_t node_id) {
        std::string node_string = "";
        auto fs = label_store_.return_string_pointer(node_id);
        if(fs == nullptr) return node_string;
        for(uint64_t i=0;; i++) {
            if(fs[i] == 0x00) break;
            node_string += fs[i];
        }
        return node_string;
    }

    void insert_new_dic_using_middle_string(uint64_t node_id, std::vector<std::string>& part_keys) {
        // std::cout << "--- insert_new_dic ---" << std::endl;
        std::string insert_key = restore_insert_string_using_middle_string(node_id, part_keys);
        int* ptr = update_new(make_char_range(insert_key), hash_trie_.get_root(), 0); // 先頭から比較
        *ptr = 1;
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
        bool exist_zero_blanch   = (fp[node_id][0].match == 0 ? true : false);
        bool do_zero_blanch_last = false; // zero分岐を最後に処理する場合は、true
        if(exist_zero_blanch) {
            if(fp[node_id][0].cnt <= bn[node_id]) do_zero_blanch_last = true;
        }

        if(exist_zero_blanch) {
            std::sort(fp[node_id].begin()+1, fp[node_id].end(), [] (auto l, auto r) {
                return l.cnt > r.cnt;
            });
        } else {
            std::sort(fp[node_id].begin(), fp[node_id].end(), [] (auto l, auto r) {
                return l.cnt > r.cnt;
            });
        }

        for(uint64_t i=0; i < fp[node_id].size(); i++) {
            uint64_t pos = (do_zero_blanch_last ? (i + 1) % fp[node_id].size() : i);
            std::sort(fp[node_id][pos].children.begin(), fp[node_id][pos].children.end(), [] (auto l, auto r) {
                return l.second > r.second;
            });
            for(auto s : fp[node_id][pos].children) {
                require_centroid_path_order_and_insert_dictionaly(fp, s.first, bn, check_bottom, compare_str);
                insert_new_dic(s.first, check_bottom[s.first], compare_str);
            }
        }

        if(node_id == hash_trie_.get_root()) insert_new_dic(node_id, check_bottom[node_id], compare_str);
    }

    // centroid_path_を求めて、新しい辞書に格納する
    void require_centroid_path_order_and_insert_dictionaly_not_fork_info(std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                                            uint64_t node_id,
                                                            const std::vector<uint64_t>& blanch_num_except_zero,
                                                            const std::vector<uint64_t>& cnt_leaf,
                                                            const std::vector<bool>& check_bottom,
                                                            // std::vector<std::string>& part_keys) {
                                                            std::string& compare_str) {
        // 一番下まで、たどり着いた時の処理
        if(children[node_id].size() == 0) return;

        std::vector<std::pair<uint64_t, uint64_t>> match_per_leaf_num;  // first:match, second:cnt_leaf[]
        std::map<uint64_t, uint64_t> start_pos;                         // children[node_id]内のmatchのスタート位置を保存
        uint64_t pre_match = UINT64_MAX;                                // 以前、比較したmatch(文字列比較失敗位置)を保存
        uint64_t zero_blanch_num = 0;                                   // ノード(node_id)のzero分岐の総数
        uint64_t children_size = children[node_id].size();              // 子ノードの数

        // matchごとに、分岐後の葉の数をカウント
        // children[]内は、matchの値ごとにソートされている
        // children[](first:match, second:next_id)
        // ex. children[i] : ({0, 4}, {0, 5}, {3, 54}, {13, 23})
        for(uint64_t i=0; i < children_size; i++) {
            std::pair<uint64_t, uint64_t> child = children[node_id][i];
            if(child.first != pre_match) {
                match_per_leaf_num.push_back(std::pair{child.first, cnt_leaf[child.second]});
                pre_match = child.first;
                start_pos[pre_match] = i + 1;
                if(child.first == 0) zero_blanch_num += cnt_leaf[child.second];
            } else {
                match_per_leaf_num[match_per_leaf_num.size()-1].second += cnt_leaf[child.second];
                if(child.first == 0) zero_blanch_num += cnt_leaf[child.second];
            }
        }

        bool exist_zero_blanch = (zero_blanch_num == 0 ? false : true);
        // zero分岐が存在するときは、その部分以外をソートする
        if(exist_zero_blanch) {
            if(match_per_leaf_num[0].first != 0) std::cout << match_per_leaf_num[0].first << std::endl;
            std::sort(match_per_leaf_num.begin()+1, match_per_leaf_num.end(), [] (auto l, auto r) {
                return l.second > r.second;
            });
        } else {
            std::sort(match_per_leaf_num.begin(), match_per_leaf_num.end(), [] (auto l, auto r) {
                return l.second > r.second;
            });
        }

        // zero分岐を最後に処理するのかを判定
        bool do_zero_blanch_last = false;
        if(exist_zero_blanch) {
            if(zero_blanch_num <= blanch_num_except_zero[node_id]) do_zero_blanch_last = true;
        }

        // 分岐後の葉の数が多い順に処理していく
        // 初めに、対象とするmatchの子ノードとそれ以降の葉の数を保存
        uint64_t match_per_leaf_num_size = match_per_leaf_num.size();
        for(uint64_t i=0; i < match_per_leaf_num_size; i++) {
            uint64_t pos = (do_zero_blanch_last ? (i+1)%match_per_leaf_num_size : i);   // zero分岐を最後に処理するのかによって、posの値が変化する
            std::vector<std::pair<uint64_t, uint64_t>> children_shelter;                // first:cnt_leaf[], second:next_id
            pre_match = match_per_leaf_num[pos].first;
            for(uint64_t j=start_pos[pre_match]-1; j < children_size; j++) {
                if(pre_match != children[node_id][j].first) break;
                uint64_t next_id = children[node_id][j].second;
                children_shelter.push_back({cnt_leaf[next_id], next_id});
            }
            std::sort(children_shelter.begin(), children_shelter.end(), [] (auto l, auto r) {
                return l.first > r.first;
            });
            for(auto& s : children_shelter) {
                require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, s.second, blanch_num_except_zero, cnt_leaf, check_bottom, compare_str);
                // insert_new_dic(s.second, check_bottom[s.first], compare_str); // 新しい辞書に登録する
                // require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, s.second, blanch_num_except_zero, cnt_leaf, check_bottom, part_keys);
                // insert_new_dic_using_middle_string(s.second, part_keys);
            }
        }

        // if(node_id == hash_trie_.get_root()) insert_new_dic(node_id, check_bottom[node_id], compare_str);
        // if(node_id == hash_trie_.get_root()) insert_new_dic_using_middle_string(node_id, part_keys);
    }

    // 特定のノードから、get_root()までの文字列を復元する
    std::string restore_insert_string(uint64_t node_id) {
        std::string insert_string = ""; // ここに文字列を格納して、新しい辞書に挿入する

        // とりあえず、文字列を復元する
        auto fs = label_store_.return_string_pointer(node_id);
        if(fs == nullptr) return insert_string;
        for(uint64_t i=0;; i++) {
            if(fs[i] == 0x00) break;
            insert_string += fs[i];
        }

        // get_root()まで、文字列を復元する
        while(node_id != hash_trie_.get_root()) {
            auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id); // 親ノードとsymbを取得
            auto [c, match] = restore_symb_(symb);                         // symbから、遷移に失敗した箇所とlabelを取得する

            insert_string = c + insert_string;

            uint64_t dummy_step = 0; // ダミーノードの数を数える
            while(1) {
                fs = label_store_.return_string_pointer(parent);
                if(fs == nullptr) {
                    dummy_step++;
                    auto [tmp1, tmp2] = hash_trie_.get_parent_and_symb(parent);
                    parent = tmp1;
                } else {
                    break;
                }
            }

            match += dummy_step * lambda_; // スキップした回数分だけ足す

            if(match != 0) {
                fs = label_store_.return_string_pointer(parent);
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

    // 特定のノードに対応する文字列を復元するのは同じ
    // 途中の部分を保存しておく
    std::string restore_insert_string_using_middle_string(uint64_t node_id, std::vector<std::string>& part_keys) {
        std::string insert_string = "";
        std::deque<uint64_t> save_route;
        auto str_pointer = label_store_.return_string_pointer(node_id);
        for(uint64_t i=0;; i++) {
            if(str_pointer[i] == 0x00) break;
            insert_string += str_pointer[i];
        }
        while(node_id != hash_trie_.get_root()) {
            auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id);
            auto [c, match] = restore_symb_(symb);
            insert_string = c + insert_string;

            uint64_t dummy_step = 0;
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
            match += dummy_step * lambda_;
            
            if(match != 0) { // matchがゼロの時は、文字列比較の一文字目で遷移しているので
                auto str1 = label_store_.return_string_pointer(parent);
                std::string str2 = "";
                for(uint64_t j=0; j < match; j++) {
                    str2 += str1[j];
                }
                insert_string = str2 + insert_string;
                part_keys[node_id] = str2 + c;
            } else {
                part_keys[node_id] = c;
            }
            save_route.push_back(node_id);
            node_id = parent;
             if(part_keys[node_id].size() != 0) { // 一度、通ったことのあるノードだった場合
                insert_string = part_keys[node_id] + insert_string;
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
        return insert_string;
    }

    // データをリセットするための関数
    void data_reset() {
        hash_trie_.reset_data_();
        label_store_.reset_data_();
    }

    // trie.hppに書いてあった処理をこちらで行う
    // 二分探索の実装(lower_bound)
    std::pair<bool, uint64_t> BinarySearch(const std::vector<info_fp>& data, uint64_t match) {
        uint64_t size = data.size();
        if(size == 0) return {false, 0};
        int64_t left = 0;
        int64_t right = size - 1;
        int64_t mid;
        while(left <= right) {
            mid = (left + right) / 2;
            if(data[mid].match == match) {
                return {true, mid};
            } else if(data[mid].match > match) {
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        return {false, uint64_t(left)};
    }

    // 全探索の実装
    std::pair<bool, uint64_t> FullSearch(const std::vector<info_fp>& data, uint64_t match) {
        uint64_t size = data.size();
        for(uint64_t i=0; i < size; i++) {
            if(data[i].match == match) {
                return {true, i};
            }
        }
        return {false, size};
    }

    // CP順を求めるために部分的な情報を求める
    std::tuple<std::vector<std::vector<info_fp>>, std::vector<uint64_t>, std::vector<bool>> return_partial_CP_info() {
        // std::cout << "--- return_partial_CP_info ---" << std::endl;
        uint64_t table_size = hash_trie_.capa_size();
        std::vector<parent_info> parent(table_size);            // 親ノード、いくつめの分岐、どの文字かを保存
        std::vector<uint64_t> partial_num(table_size, 0);       // 子の数を格納するための配列(自身の数も含む)
        
        std::vector<std::vector<info_fp>> fork_pos(table_size); // std::vector<std::vector<std::pair<uint64_t, uint64_t>>> fork_pos(table_size); // それぞれの分岐位置で個数を求めるためのもの(分岐位置)
        std::vector<uint64_t> cnt_leaf(table_size, 0);          // それぞれのノードから繋がっている葉ノードの数をカウント
        std::vector<uint64_t> all_branch(table_size, 0);        // 特定のノード以下に何個のノードが存在するのか

        // O(n)で、親の位置、子の数(ノード番号も)を数える
        // 現在はダミーノードの数も計測してしまっている → 完了
        for(uint64_t i=0; i < table_size; i++) {
            if(hash_trie_.is_use_table(i) != 0) {
                auto [p, label] = hash_trie_.get_parent_and_symb(i);                            // 親と遷移情報の取得
                if(label == 255) continue;                                                      // ダミーノードをスキップ
                auto [c, match] = std::pair{uint8_t(restore_codes_[label % 256]), label/256};   // 遷移文字と分岐位置を取得
                while(1) { // 親がダミーノードの時の処理
                    // label_store_からポインタを参照
                    auto string_pointer = label_store_.return_string_pointer(p);
                    if(string_pointer != nullptr) break;
                    // ダミーノードの時、上の処理を繰り返す（ダミーではない親を探す）
                    auto [tmp_parent, tmp_label] = hash_trie_.get_parent_and_symb(p);
                    match += lambda_;
                    p = tmp_parent;
                }
                partial_num[i] += 1; // 自身の数をカウントする
                partial_num[p] += 1; // 子から親の数をカウントする
                parent[i].parent = p;
                parent[i].match = match;
                parent[i].c = c;
            }
        }

        // queueを使用して、一番下のものから処理していく
        // 初めに、一番下のノードのみを取得 (O(n))
        std::queue<uint64_t> que;
        std::vector<bool> check_bottom(table_size, false);
        // 対象のデータを集めてくる
        for(uint64_t i=0; i < table_size; i++) {
            if(partial_num[i] == 1) {
                que.push(i);
                check_bottom[i] = true;
            }
        }

        // queueに追加した順に処理（DynPDT上の一番下のノードから）
        check_bottom[hash_trie_.get_root()] = true; // 先頭文字列も底とした方が計算しやすいため
        while(!que.empty()) {                       // 追加された順に処理
            uint64_t q = que.front();
            que.pop();
            auto [node_id, match, c] = parent[q];
            // if(c != 0) cnt_leaf[q] += 1; // ダミーノード以外のとき、加算(上で処理済み)
            cnt_leaf[q] += 1;
            uint64_t p = node_id;
            cnt_leaf[p] += cnt_leaf[q];
            if(match != 0) {    // match=0は、親ノードからの分岐位置が同じことを示している
                all_branch[p] += cnt_leaf[q];
            } else {            // もう一つ上の親に足す
                if(p != hash_trie_.get_root()) {
                    auto [pp, m , cc] = parent[p];
                    all_branch[pp] += cnt_leaf[q];
                }
            }
            // fork_posに対して、分岐の位置に対して、葉がいくつあるのかをカウント
            // ソートした状態を保つために、このような処理をしている(ここに時間がかかっている)
            auto [flag, pos] = BinarySearch(fork_pos[p], match);    // matchの失敗位置が以前に存在していたのかなど
            if(flag) {
                fork_pos[p][pos].cnt += cnt_leaf[q];
                fork_pos[p][pos].children.push_back({q, cnt_leaf[q]});
            } else {
                fork_pos[p].insert(fork_pos[p].begin()+pos, info_fp{match, cnt_leaf[q]});
                fork_pos[p][pos].children.push_back({q, cnt_leaf[q]});
            }

            partial_num[p]--;
            if(partial_num[p] == 1) { // 子供の処理がすべて終了すると追加
                que.push(p);
            }
        }

        // std::map<uint64_t, uint64_t> mp;
        // for(uint64_t i=0; i < table_size; i++) {
        //     if(fork_pos[i].size() == 0) continue;
        //     for(uint64_t j=0; j < fork_pos[i].size(); j++) {
        //         mp[fork_pos[i][j].children.size()] += 1;
        //     }
        // }

        // for(auto m : mp) std::cout << m.first << " : " << m.second << std::endl;

        return std::tuple{fork_pos, all_branch, check_bottom};
    }

    // mapを使用して、CPの情報を求める
    std::tuple<std::vector<std::multimap<uint64_t, std::pair<uint64_t, uint64_t>>>, std::vector<uint64_t>, std::vector<bool>> return_partial_CP_info_using_map() {
        // std::cout << "--- return_partial_CP_info_using_map ---" << std::endl;
        uint64_t table_size = hash_trie_.capa_size();
        std::vector<parent_info> parent(table_size);            // 親ノード、いくつめの分岐、どの文字かを保存
        std::vector<uint64_t> partial_num(table_size, 0);       // 子の数を格納するための配列(自身の数も含む)
        
        // std::vector<std::vector<info_fp>> fork_pos(table_size); // std::vector<std::vector<std::pair<uint64_t, uint64_t>>> fork_pos(table_size); // それぞれの分岐位置で個数を求めるためのもの(分岐位置)
        // std::vector<std::map<uint64_t, std::vector<std::pair<uint64_t, uint64_t>>>> fork_pos(table_size);
        std::vector<std::multimap<uint64_t, std::pair<uint64_t, uint64_t>>> fork_pos(table_size);
        std::vector<std::map<uint64_t, uint64_t>> every_match_leaf(table_size);
        std::vector<uint64_t> cnt_leaf(table_size, 0);          // それぞれのノードから繋がっている葉ノードの数をカウント
        std::vector<uint64_t> all_branch(table_size, 0);        // 特定のノード以下に何個のノードが存在するのか

        // O(n)で、親の位置、子の数(ノード番号も)を数える
        // 現在はダミーノードの数も計測してしまっている → 完了
        for(uint64_t i=0; i < table_size; i++) {
            if(hash_trie_.is_use_table(i) != 0) {
                auto [p, label] = hash_trie_.get_parent_and_symb(i);                            // 親と遷移情報の取得
                if(label == 255) continue;                                                      // ダミーノードをスキップ
                auto [c, match] = std::pair{uint8_t(restore_codes_[label % 256]), label/256};   // 遷移文字と分岐位置を取得
                while(1) { // 親がダミーノードの時の処理
                    // label_store_からポインタを参照
                    auto string_pointer = label_store_.return_string_pointer(p);
                    if(string_pointer != nullptr) break;
                    // ダミーノードの時、上の処理を繰り返す（ダミーではない親を探す）
                    auto [tmp_parent, tmp_label] = hash_trie_.get_parent_and_symb(p);
                    match += lambda_;
                    p = tmp_parent;
                }
                partial_num[i] += 1; // 自身の数をカウントする
                partial_num[p] += 1; // 子から親の数をカウントする
                parent[i].parent = p;
                parent[i].match = match;
                parent[i].c = c;
            }
        }

        // queueを使用して、一番下のものから処理していく
        // 初めに、一番下のノードのみを取得 (O(n))
        std::queue<uint64_t> que;
        std::vector<bool> check_bottom(table_size, false);
        // 対象のデータを集めてくる
        for(uint64_t i=0; i < table_size; i++) {
            if(partial_num[i] == 1) {
                que.push(i);
                check_bottom[i] = true;
            }
        }

        // queueに追加した順に処理（DynPDT上の一番下のノードから）
        check_bottom[hash_trie_.get_root()] = true; // 先頭文字列も底とした方が計算しやすいため
        while(!que.empty()) {                       // 追加された順に処理
            uint64_t q = que.front();
            que.pop();
            auto [node_id, match, c] = parent[q];
            // if(c != 0) cnt_leaf[q] += 1; // ダミーノード以外のとき、加算(上で処理済み)
            cnt_leaf[q] += 1;
            uint64_t p = node_id;
            cnt_leaf[p] += cnt_leaf[q];
            if(match != 0) all_branch[p] += cnt_leaf[q];    // match=0は、親ノードからの分岐位置が同じことを示している
            // fork_posに対して、分岐の位置に対して、葉がいくつあるのかをカウント
            every_match_leaf[p][match] += cnt_leaf[q];
            // fork_pos[p][match].push_back({q, cnt_leaf[q]});
            fork_pos[p].insert({p, {q, cnt_leaf[q]}});

            partial_num[p]--;
            if(partial_num[p] == 1) { // 子供の処理がすべて終了すると追加
                que.push(p);
            }
        }

        return std::tuple{fork_pos, all_branch, check_bottom};
    }

    // CP順を求めるために部分的な情報を求める
    void return_partial_CP_info_not_fork_info(std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                              std::vector<uint64_t>& blanch_num_except_zero,
                                              std::vector<uint64_t>& cnt_leaf_per_node,
                                              std::vector<bool>& is_node_bottom) {
        // std::cout << "--- return_partial_CP_info_not_fork_info ---" << std::endl;
        uint64_t table_size = hash_trie_.capa_size();
        std::vector<parent_info> parent(table_size);        // 親ノード、文字列比較での失敗位置、遷移文字を保存
        std::vector<uint64_t> partial_num(table_size, 0);   // 子の数を格納するための配列(自身の数も含む)
        
        cnt_leaf_per_node.resize(table_size, 0);            // それぞれのノードから繋がっている葉ノードの数をカウント
        blanch_num_except_zero.resize(table_size, 0);       // 特定のノード以下に何個のノードが存在するのか
        children.resize(table_size);                        // 子集合を保存

        // O(n)で、親の位置、子の数(ノード番号も)を数える
        // 現在はダミーノードの数も計測してしまっている → 完了
        partial_num[hash_trie_.get_root()] += 1;
        for(uint64_t i=hash_trie_.get_root()+1; i < table_size; i++) {
            if(hash_trie_.is_use_table(i)) {                                                    // テーブルが使用されているとき
                auto [p, label] = hash_trie_.get_parent_and_symb(i);                            // 親と遷移情報の取得
                if(label == 255) continue;                                                      // ダミーノードをスキップ
                auto [c, match] = std::pair{uint8_t(restore_codes_[label % 256]), label/256};   // 遷移文字と分岐位置を取得
                while(label_store_.return_string_pointer(p) == nullptr) {                       // 親がダミーノードの時の処理(ダミーノードに文字列は保存されていないため)
                    auto [tmp_parent, tmp_label] = hash_trie_.get_parent_and_symb(p);
                    match += lambda_;
                    p = tmp_parent;
                }
                partial_num[i] += 1; // 自身の数をカウントする
                partial_num[p] += 1; // 子から親の数をカウントする
                children[p].push_back({match, i});
                parent[i].parent = p;
                parent[i].match = match;
                parent[i].c = c;
            }
        }

        // childrenをmatch毎にソートする
        for(uint64_t i=0; i < table_size; i++) {
            std::sort(children[i].begin(), children[i].end(), [] (auto l, auto r) {
                return l.first < r.first;
            });
        }

        // queueを使用して、一番下のものから処理していく
        // 初めに、一番下のノードのみを取得 (O(n))
        std::queue<uint64_t> que;
        is_node_bottom.resize(table_size, false);
        // 対象のデータを集めてくる
        for(uint64_t i=0; i < table_size; i++) {
            if(partial_num[i] == 1) {
                que.push(i);
                is_node_bottom[i] = true;
            }
        }

        // queueに追加した順に処理（DynPDT上の一番下のノードから）
        is_node_bottom[hash_trie_.get_root()] = true;   // 先頭文字列も底とした方が計算しやすいため
        while(!que.empty()) {                           // 追加された順に処理
            uint64_t q = que.front();
            que.pop();
            auto [node_id, match, c] = parent[q];
            cnt_leaf_per_node[q] += 1;
            uint64_t p = node_id;
            cnt_leaf_per_node[p] += cnt_leaf_per_node[q];
            if(match != 0) blanch_num_except_zero[p] += cnt_leaf_per_node[q];    // match=0は、親ノードからの分岐位置が同じことを示している
            else { // もう一つ上の親に足す(zero分岐のトライ上で考えるとそのようになるから)
                if(p != hash_trie_.get_root()) {
                    auto [pp, m , cc] = parent[p];
                    blanch_num_except_zero[pp] += cnt_leaf_per_node[q];
                }
            }
            partial_num[p]--;
            if(partial_num[p] == 1) { // 子供の処理がすべて終了すると追加
                que.push(p);
            }
        }
    }

    // トポロジカルソートを呼び出す
    void call_topo() {
        std::cout << "--- call_topo ---" << std::endl;
        
        // auto [fork_info, blanch_num, check_bottom] = hash_trie_.return_partial_CP_info(restore_codes_); // trie.hppで計算
        // auto [fork_info, blanch_num, check_bottom] = return_partial_CP_info(); // vectorですべて計算
        // auto [fork_info, blanch_num, check_bottom] = return_partial_CP_info_using_map(); // map.hppで計算
        // std::cout << "size : " << fork_info.size() << std::endl;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children;
        std::vector<uint64_t> blanch_num_except_zero;
        std::vector<uint64_t> cnt_leaf_per_node;
        std::vector<bool> is_node_bottom;
        return_partial_CP_info_not_fork_info(children, blanch_num_except_zero, cnt_leaf_per_node, is_node_bottom); // fork_infoを使用しない方法
        std::cout << "size : " << children.size() << std::endl;

        // fork_infoの情報を元に、CP順を求め、新しい辞書に格納していく
        // data_reset();
        // hash_trie_.expand_tmp_table();
        // label_store_.expand_tmp_ptrs();
        std::string compare_str = "";
        // require_centroid_path_order_and_insert_dictionaly(fork_info, hash_trie_.get_root(), blanch_num, check_bottom, compare_str);
        // std::vector<std::string> part_keys(hash_trie_.capa_size());
        // require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, hash_trie_.get_root(), blanch_num_except_zero, cnt_leaf_per_node, is_node_bottom, part_keys);
        require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, hash_trie_.get_root(), blanch_num_except_zero, cnt_leaf_per_node, is_node_bottom, compare_str);
        // hash_trie_.move_table();    // 辞書の移動
        // label_store_.move_ptrs();
        // hash_trie_.set_first_insert(false);

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

    void restore_string(std::vector<std::string>& keys) {
        std::vector<std::string> part_keys(hash_trie_.capa_size()); // 部分的な文字列を保存しておく
        
        uint64_t first_node = hash_trie_.get_root();
        std::string restore_str;

        // 先頭ノードに対する処理
        auto fs = label_store_.return_string_pointer(first_node);
        for(uint64_t i=0;;i++) {
            if(fs[i] == 0x00) break; // 判定処理としてint型で4つ分使用して1を格納しているので、この判定
            restore_str += fs[i];
        }
        keys.push_back(restore_str);

        // 先頭ノード以外のノードに対する処理
        for(int i=first_node+1; i < hash_trie_.capa_size(); i++) {
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
            keys.push_back(restore_str);
        }
    }

    void reset() {
        size_ = 0;
        for(int i=0; i < 256; i++) {
            codes_[i] = static_cast<uint8_t>(UINT8_MAX);
            restore_codes_[i] = static_cast<uint8_t>(UINT8_MAX);
        }
        num_codes_ = 0;
        codes_[0] = static_cast<uint8_t>(num_codes_++);
    }

    // 再配置するための関数
    template <class It>
    void insert_by_centroid_path_order(It begin, It end, int depth) {
        assert(end-begin > 0);
        if (end-begin == 1) {
            int* ptr = update(*begin);
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
            insert_by_centroid_path_order(f, t, depth+1);
        }
    }

    void call_restore_string_CP() {
        std::cout << "--- call_restore_string_CP ---" << std::endl;
        std::vector<std::string> restore_keys;
        restore_string(restore_keys);
        std::cout << hash_trie_.capa_size() << ", " << label_store_.num_ptrs() << std::endl;
        sort(restore_keys.begin(), restore_keys.end()); // 取得した文字列をソートする
        hash_trie_.expand_restore_string();
        label_store_.expand_restore_string();
        reset();
        insert_by_centroid_path_order(restore_keys.begin(), restore_keys.end(), 0); // 新しい辞書にCP順で格納する
    }

    // compute_...関数内での使用する一時変数を記憶するための構造体
    struct compute_temporary_set_of_variable {
        compute_temporary_set_of_variable(uint64_t size) {
            parent.resize(size);
            partial_num.resize(size, 0);
            children_size.resize(size, 0);
        }

        std::vector<parent_info> parent;    // 親ノード、文字列比較での失敗位置、遷移文字を保存
        std::vector<uint64_t> partial_num;  // 子の数を格納するための配列(自身の数も含む)
        std::vector<uint64_t> children_size;// reserveで要素を確保するために使用
    };

    // ノードのつながりと分岐数(葉の数)を調べるための関数
    void compute_node_connect_and_blanch_num(std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                            //  std::vector<uint64_t>& blanch_num_except_zero,
                                             std::vector<uint64_t>& cnt_leaf_per_node) {
        uint64_t table_size = hash_trie_.capa_size();
        cnt_leaf_per_node.resize(table_size, 0);            // それぞれのノードから繋がっている葉ノードの数をカウント
        // blanch_num_except_zero.resize(table_size, 0);    // 特定のノード以下に何個のノードが存在するのか
        children.resize(table_size);                        // 子集合を保存
        compute_temporary_set_of_variable variable(table_size);

        // O(n)で、親の位置、子の数(ノード番号も)を数える
        // 現在はダミーノードの数も計測してしまっている → 完了
        variable.partial_num[hash_trie_.get_root()] += 1;
        for(uint64_t i=hash_trie_.get_root()+1; i < table_size; i++) {
            if(hash_trie_.is_use_table(i)) {                                                    // テーブルが使用されているとき
                auto [p, label] = hash_trie_.get_parent_and_symb(i);                            // 親と遷移情報の取得
                if(label == 255) continue;                                                      // ダミーノードをスキップ
                auto [c, match] = std::pair{uint8_t(restore_codes_[label % 256]), label/256};   // 遷移文字と分岐位置を取得
                while(label_store_.return_string_pointer(p) == nullptr) {                       // 親がダミーノードの時の処理(ダミーノードに文字列は保存されていないため)
                    auto [tmp_parent, tmp_label] = hash_trie_.get_parent_and_symb(p);
                    match += lambda_;
                    p = tmp_parent;
                }
                assert(p != nil_id);
                variable.partial_num[p] += 1; // 子から親の数をカウントする
                variable.children_size[p]++;
                variable.partial_num[i] += 1; // 自身の数をカウントする
                variable.parent[i].parent = p;
                variable.parent[i].match = match;
                variable.parent[i].c = c;
            }
        }
        for(uint64_t i=hash_trie_.get_root()+1; i < table_size; i++) {
          if (!hash_trie_.is_use_table(i)) continue;
          children[i].reserve(variable.children_size[i]);
        }
        for(uint64_t i=hash_trie_.get_root()+1; i < table_size; i++) {
          if (!hash_trie_.is_use_table(i)) continue;
          auto p = variable.parent[i].parent;
          auto match = variable.parent[i].match;
          children[p].emplace_back(match, i);
        }

        // queueを使用して、一番下のものから処理していく
        // 初めに、一番下のノードのみを取得 (O(n))
        std::vector<uint64_t> que;
        que.reserve(hash_trie_.size());
        // 対象のデータを集めてくる
        for(uint64_t i=0; i < table_size; i++) {
            if(variable.partial_num[i] == 1) {
                que.emplace_back(i);
            }
        }

        // queueに追加した順に処理（DynPDT上の一番下のノードから）
        for(uint64_t k = 0; k < que.size(); ++k) {                           // 追加された順に処理
            uint64_t q = que[k];
            auto [node_id, match, c] = variable.parent[q];
            cnt_leaf_per_node[q] += 1;
            uint64_t p = node_id;
            assert(p != nil_id);
            cnt_leaf_per_node[p] += cnt_leaf_per_node[q];
            // if(match != 0) blanch_num_except_zero[p] += cnt_leaf_per_node[q];    // match=0は、親ノードからの分岐位置が同じことを示している
            // else { // もう一つ上の親に足す(zero分岐のトライ上で考えるとそのようになるから)
            //     if(p != hash_trie_.get_root()) {
            //         auto [pp, m , cc] = parent[p];
            //         blanch_num_except_zero[pp] += cnt_leaf_per_node[q];
            //     }
            // }
            variable.partial_num[p]--;
            if(variable.partial_num[p] == 1) { // 子供の処理がすべて終了すると追加
                que.emplace_back(p);
            }
        }
    }

    struct temporary_set_of_variable_require_backward {
        temporary_set_of_variable_require_backward(uint64_t pm, uint64_t zbn, uint64_t cs, uint64_t p, uint64_t cus) :
                        pre_match{pm}, zero_blanch_num{zbn}, children_size{cs}, pos{p}, cumulative_sum{cus} {}

        uint64_t pre_match;                 // 以前、比較したmatch(文字列比較失敗位置)を保存
        uint64_t zero_blanch_num;           // ノード(node_id)のzero分岐の総数
        uint64_t children_size;             // 子ノードの数
        uint64_t pos;                       // 位置を示す際に使用
        uint64_t cumulative_sum;            // 累積和(0分岐を除く) 
        uint64_t match_per_leaf_num_size;   // vector match_per_leaf_numのサイズ
        std::vector<std::pair<uint64_t, uint64_t>> match_per_leaf_num;  // first:match, second:cnt_leaf_per_node[]
        std::vector<uint64_t> backward_cumulative_num;
    };

    // 後方累積和を使用して、CPD順を求めたもの
    template<class Map>
    bool require_centroid_path_order_and_insert_dictionary(
                                                            Map& new_map,
                                                            std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                                            uint64_t node_id,
                                                            const std::vector<uint64_t>& cnt_leaf_per_node,
                                                            uint64_t common_prefix_length,
                                                            std::string& store_string) {
        // std::cout << "common_prefix_length : " << common_prefix_length << std::endl;
        // 一番下まで、たどり着いた時の処理
        if(children[node_id].size() == 0) return true;

        std::map<uint64_t, uint64_t> start_pos;                         // children[node_id]内のmatchのスタート位置を保存
        temporary_set_of_variable_require_backward variable(UINT64_MAX, 0, children[node_id].size(), 0, 0);

        // 現在のchildrent[node_id]内は以下のようになっている
        // ex. ({3, 54}, {0, 4}, {13, 23}, {0, 5})
        // ソートして、下のように修正
        // ex. ({0, 4}, {0, 5}, {3, 54}, {13, 23})
        std::sort(children[node_id].begin(), children[node_id].end(), [] (auto l, auto r) {
            return l.first < r.first;
        });
        
        // TODO: children の match を座標圧縮して処理する
        std::vector<uint64_t> matches(children[node_id].size());
        std::transform(children[node_id].begin(), children[node_id].end(), matches.begin(), [](auto& c) {return c.first;});
        matches.erase(std::unique(matches.begin(), matches.end()), matches.end());

        // メモリを最初に確保し、childrenの先頭から順に2分探索で探す
        // matchごとに、分岐後の葉の数をカウント
        // children           : ({0, 4}, {0, 5}, {3, 54}, {13, 23})
        // match_per_leaf_num : ({0, 8}, {3, 7}, {13, 15})
        // start_pos          : ({0, 1}, {3, 3}, {13, 4})
        variable.match_per_leaf_num.resize(matches.size());
        variable.pre_match = children[node_id][0].first;
        start_pos[variable.pre_match] = 1;
        for(uint64_t i=0; i < variable.children_size; i++) {
            auto [match, next_id] = children[node_id][i];
            if(variable.pre_match != match) {
                variable.pre_match = match;
                start_pos[match] = i + 1;
                variable.pos++;
            }
            // uint64_t pos = std::lower_bound(matches.begin(), matches.end(), match) - matches.begin();
            variable.match_per_leaf_num[variable.pos].first = match;
            variable.match_per_leaf_num[variable.pos].second += cnt_leaf_per_node[next_id];
            if(match == 0) variable.zero_blanch_num += cnt_leaf_per_node[next_id];  // 0分岐の和
            else variable.cumulative_sum += cnt_leaf_per_node[next_id];             // 0分岐以外の累積和
        }

        bool exist_zero_blanch = variable.zero_blanch_num != 0;
        // // 後方累積和を計算する
        uint64_t size;
        if(exist_zero_blanch) { // zero分岐が存在する場合、zero分岐以外の和を取る必要があるため
            size = matches.size() - 1;
            variable.backward_cumulative_num.resize(size, 0);
            for(int64_t i=size-1; i >= 0; i--) {
                if(i == (size-1)) variable.backward_cumulative_num[i] = variable.match_per_leaf_num[i+1].second;
                else variable.backward_cumulative_num[i] = variable.backward_cumulative_num[i+1] + variable.match_per_leaf_num[i+1].second;
            }
        } else {
            size = matches.size();
            variable.backward_cumulative_num.resize(size, 0);
            for(int64_t i=size-1; i >= 0; i--) {
                if(i == (size-1)) variable.backward_cumulative_num[i] = variable.match_per_leaf_num[i].second;
                else variable.backward_cumulative_num[i] = variable.backward_cumulative_num[i+1] + variable.match_per_leaf_num[i].second;
            }
        }

        // 処理する順番を決める
        std::queue<uint64_t> front;
        std::stack<uint64_t> back;

        if(exist_zero_blanch) { // zero分岐を最初、または最後に処理するのか
            if(variable.zero_blanch_num <= variable.cumulative_sum) back.push(0);
            else front.push(0);
        }

        if(size > 0) {
            for(uint64_t i=0; i < size; i++) {
                uint64_t pos = exist_zero_blanch ? (i+1) : (i);
                if(i == (size-1)) {
                    front.push(pos);
                } else {
                    if(variable.match_per_leaf_num[pos].second > variable.backward_cumulative_num[i+1]) {
                        front.push(pos);
                    } else {
                        back.push(pos);
                    }
                }
            }
        }        
        std::string compared_string = ""; // 帰りがけで取得した文字列を保存しておく文字列
        while(!front.empty()) {
            uint64_t num = front.front();
            front.pop();
            variable.pre_match = variable.match_per_leaf_num[num].first;
            std::vector<std::pair<uint64_t, uint64_t>> children_shelter;    // first:cnt_leaf_per_node[], second:next_id
            for(uint64_t i=start_pos[variable.pre_match]-1; i < variable.children_size; i++) {
                if(variable.pre_match != children[node_id][i].first) break;
                uint64_t next_id = children[node_id][i].second;
                children_shelter.emplace_back(cnt_leaf_per_node[next_id], next_id);
            }

            // 特定の分岐内の内、分岐後の葉の数が多いnext_idを取得
            auto max_pos = std::max_element(children_shelter.begin(), children_shelter.end());
            uint64_t next_common_prefix_length = common_prefix_length+variable.pre_match+1; // 行き掛けで渡す共通のprefixの長さを伝搬する
            bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, max_pos->second, cnt_leaf_per_node, next_common_prefix_length, store_string);
            // ここで、original_stringからcompare_prefix_length分のみを抜き出し、insert_new_dicの関数の引数とする
            if(is_leaf_node) {
                insert_new_dic_CPD(new_map, max_pos->second, 0, store_string);
            } else {
                insert_new_dic_CPD(new_map, max_pos->second, next_common_prefix_length, store_string);
            }
            // insert_new_dic(new_map, max_pos->second, "");
            for(auto& s : children_shelter) {
                if(s.second == max_pos->second) continue;
                bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string);
                if(is_leaf_node) {
                    insert_new_dic_CPD(new_map, s.second, 0, store_string);
                } else {
                    insert_new_dic_CPD(new_map, s.second, next_common_prefix_length, store_string);
                }
            }
        }
        while(!back.empty()) {
            uint64_t num = back.top();
            back.pop();
            variable.pre_match = variable.match_per_leaf_num[num].first;
            std::vector<std::pair<uint64_t, uint64_t>> children_shelter;                // first:cnt_leaf_per_node[], second:next_id
            for(uint64_t i=start_pos[variable.pre_match]-1; i < variable.children_size; i++) {
                if(variable.pre_match != children[node_id][i].first) break;
                uint64_t next_id = children[node_id][i].second;
                children_shelter.emplace_back(cnt_leaf_per_node[next_id], next_id);
            }

            // 特定の分岐内の内、分岐後の葉の数が多いnext_idを取得
            auto max_pos = std::max_element(children_shelter.begin(), children_shelter.end());
            uint64_t next_common_prefix_length = common_prefix_length+variable.pre_match+1;
            bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, max_pos->second, cnt_leaf_per_node, next_common_prefix_length, store_string);
            if(is_leaf_node) {
                insert_new_dic_CPD(new_map, max_pos->second, 0, store_string);
            } else {
                insert_new_dic_CPD(new_map, max_pos->second, next_common_prefix_length, store_string);
            }
            for(auto& s : children_shelter) {
                if(s.second == max_pos->second) continue;
                bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string);
                if(is_leaf_node) {
                    insert_new_dic_CPD(new_map, s.second, 0, store_string);
                } else {
                    insert_new_dic_CPD(new_map, s.second, next_common_prefix_length, store_string);
                }
            }
        }

        if(node_id == hash_trie_.get_root()) insert_new_dic_CPD(new_map, node_id, 0, store_string); // 新しい辞書に登録（未実装）

        return false;
    }

    value_type* dynamic_replacement(char_range& key) {
    // void dynamic_replacement() {
        // std::cout << "--- dynamic_replacement ---" << std::endl;
        // CPD_keys.clear();

        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children;   // 子ノードの集合
        // std::vector<uint64_t> blanch_num_except_zero;                    // 0分岐を除く累積和
        std::vector<uint64_t> cnt_leaf_per_node;
        compute_node_connect_and_blanch_num(children, cnt_leaf_per_node);

        // std::cout << "is_leaf : " << children[108491].size() << std::endl;
        // std::cout << cnt_leaf_per_node[108491] << std::endl;
        // std::cout << children[32246][0].first << " : " << children[32246][0].second << std::endl;

        map_check new_map(hash_trie_.capa_bits()+1);
        // std::cout << "now_map_capa_size : " << capa_size() <<std::endl;
        // std::cout << "new_map_capa_size : " << new_map.capa_size() << std::endl;
        std::string store_string = "";
        require_centroid_path_order_and_insert_dictionary(new_map, children, hash_trie_.get_root(), cnt_leaf_per_node, 0, store_string);
        std::swap(*this, new_map); // 時間がかかるので、注意

        // write_file(CPD_keys);

        return update(key);
    }

    void write_file(std::vector<std::string>& tmp_keys) {
        std::ofstream of;
        std::string filename = "../../../dataset/enwiki_restore.txt";
        of.open(filename, std::ios::out);
        for(auto b : tmp_keys) {
            of << b << std::endl;
        }
        of.close();
    }

    // 辞書に含まれている全ての文字列を復元し、正確に復元されているのかを確かめる
    void restore_and_compare(std::vector<std::string>& keys) {
        std::cout << "--- restore_and_compare ---" << std::endl;
        std::vector<std::string> restore_keys = all_key_restore_simple();
        // uint64_t table_size = hash_trie_.size();
        // int cnt = 0;
        // for(uint64_t i=hash_trie_.get_root(); i < table_size; i++) {
        //     if(!hash_trie_.is_use_table(i)) continue;
        //     std::string tmp_key = restore_insert_string(i);
        //     cnt++;
        //     if(tmp_key.size() != 0) {
        //         restore_keys.emplace_back(tmp_key);
        //     }
        // }

        write_file(restore_keys);

        if(keys.size() != restore_keys.size()) {
            std::cout << "size is different" << std::endl;
            std::cout << "keys_size : " << keys.size() << std::endl;
            std::cout << "restore_keys_size : " << restore_keys.size() << std::endl;
            // std::cout << "cnt : " << cnt << std::endl;
            return;
        }

        std::sort(keys.begin(), keys.end());
        std::sort(restore_keys.begin(), restore_keys.end());

        // 文字列の比較
        uint64_t cnt = 0;
        uint64_t size = keys.size();
        for(uint64_t i=0; i < size; i++) {
            if(keys[i] != restore_keys[i]) {
                cnt++;
                // std::cout << keys[i] << " : " << restore_keys[i] << std::endl;
                std::cout << cnt << std::endl;
                std::cout << keys[i] << std::endl;
                std::cout << restore_keys[i] << std::endl;
                const int* ptr1 = find(keys[i]);
                if(not (ptr1 != nullptr and *ptr1 == 1)) {
                    std::cout << "fails" << std::endl;
                } else {
                    std::cout << "success" << std::endl;
                }
                const int* ptr2 = find(restore_keys[i]);
                if(not (ptr2 != nullptr and *ptr2 == 1)) {
                    std::cout << "fails" << std::endl;
                } else {
                    std::cout << "success" << std::endl;
                }
            }
        }
        std::cout << "cnt : " << cnt << std::endl;

        std::cout << "ok." << std::endl;
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

    double get_label_search_time() {
        return label_store_.get_label_search_time();
    }

    double get_node_transition_search_time() {
        return hash_trie_.get_node_transition_search_time();
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

    // ハッシュテーブルの拡張が必要かを調べるための関数
    bool is_need_expand() {
        if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
            if(hash_trie_.needs_to_expand()) return true;
        }
        return false;
    }
};

}  // namespace poplar

#endif  // POPLAR_TRIE_MAP_CHECK_HPP
