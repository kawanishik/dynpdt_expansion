#ifndef POPLAR_TRIE_MAP_CHECK_HPP
#define POPLAR_TRIE_MAP_CHECK_HPP

#include <array>
#include <iostream>
#include <queue>
#include <deque>
#include <algorithm>
#include <map>

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
                    // insert_new_dic(s.first, check_bottom[s.first], compare_str);
                }

                // 1分岐以上の処理
                for(uint64_t i=1; i < fp[node_id].size(); i++) {
                    std::sort(fp[node_id][i].children.begin(), fp[node_id][i].children.end(), [] (auto l, auto r) {
                        return l.second > r.second;
                    });
                    for(auto s : fp[node_id][i].children) {
                        require_centroid_path_order_and_insert_dictionaly(fp, s.first, bn, check_bottom, compare_str);
                        // insert_new_dic(s.first, check_bottom[s.first], compare_str);
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
                        // insert_new_dic(s.first, check_bottom[s.first], compare_str);
                    }
                }

                // 0分岐の処理
                std::sort(fp[node_id][0].children.begin(), fp[node_id][0].children.end(), [] (auto l, auto r) {
                    return l.second > r.second;
                });
                for(auto s : fp[node_id][0].children) {
                    require_centroid_path_order_and_insert_dictionaly(fp, s.first, bn, check_bottom, compare_str);
                    // insert_new_dic(s.first, check_bottom[s.first], compare_str);
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
                    // insert_new_dic(s.first, check_bottom[s.first], compare_str);
                }
            }
        }
        // if(node_id == hash_trie_.get_root()) insert_new_dic(node_id, check_bottom[node_id], compare_str);
    }

    // centroid_path_を求めて、新しい辞書に格納する
    void require_centroid_path_order_and_insert_dictionaly_not_fork_info(std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                                            uint64_t node_id,
                                                            const std::vector<uint64_t>& bn,
                                                            const std::vector<bool>& check_bottom,
                                                            std::string& compare_str) {
        if(children[node_id].size() == 0) { // 一番下まで、たどり着いた時の処理
            return;
        }
        bool check_zero = false;
        if(children[node_id][0].first == 0) check_zero = true;
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
                    require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, s.first, bn, check_bottom, compare_str);
                    // insert_new_dic(s.first, check_bottom[s.first], compare_str);
                }

                // 1分岐以上の処理
                for(uint64_t i=1; i < fp[node_id].size(); i++) {
                    std::sort(fp[node_id][i].children.begin(), fp[node_id][i].children.end(), [] (auto l, auto r) {
                        return l.second > r.second;
                    });
                    for(auto s : fp[node_id][i].children) {
                        require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, s.first, bn, check_bottom, compare_str);
                        // insert_new_dic(s.first, check_bottom[s.first], compare_str);
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
                        require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, s.first, bn, check_bottom, compare_str);
                        // insert_new_dic(s.first, check_bottom[s.first], compare_str);
                    }
                }

                // 0分岐の処理
                std::sort(fp[node_id][0].children.begin(), fp[node_id][0].children.end(), [] (auto l, auto r) {
                    return l.second > r.second;
                });
                for(auto s : fp[node_id][0].children) {
                    require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, s.first, bn, check_bottom, compare_str);
                    // insert_new_dic(s.first, check_bottom[s.first], compare_str);
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
                    require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, s.first, bn, check_bottom, compare_str);
                    // insert_new_dic(s.first, check_bottom[s.first], compare_str);
                }
            }
        }
        if(node_id == hash_trie_.get_root()) insert_new_dic(node_id, check_bottom[node_id], compare_str);
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
            auto [c, match] = restore_symb_(symb); // symbから、遷移に失敗した箇所とlabelを取得する

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

            match += dummy_step * lambda_; // スキップした回数分足してあげる

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
            if(match != 0) all_branch[p] += cnt_leaf[q];    // match=0は、親ノードからの分岐位置が同じことを示している
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
    std::tuple<std::vector<std::vector<std::pair<uint64_t, uint64_t>>>, std::vector<uint64_t>, std::vector<bool>> return_partial_CP_info_not_fork_info() {
        // std::cout << "--- return_partial_CP_info ---" << std::endl;
        uint64_t table_size = hash_trie_.capa_size();
        std::vector<parent_info> parent(table_size);            // 親ノード、いくつめの分岐、どの文字かを保存
        std::vector<uint64_t> partial_num(table_size, 0);       // 子の数を格納するための配列(自身の数も含む)
        
        // std::vector<std::vector<info_fp>> fork_pos(table_size); // std::vector<std::vector<std::pair<uint64_t, uint64_t>>> fork_pos(table_size); // それぞれの分岐位置で個数を求めるためのもの(分岐位置)
        std::vector<uint64_t> cnt_leaf(table_size, 0);          // それぞれのノードから繋がっている葉ノードの数をカウント
        std::vector<uint64_t> all_branch(table_size, 0);        // 特定のノード以下に何個のノードが存在するのか
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children(table_size);// 子集合を保存

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

        // for(uint64_t i=0; i < table_size; i++) {
        //     uint64_t size = children[i].size();
        //     if(size > 5) {
        //         for(uint64_t j=0; j < size; j++) std::cout << children[i][j].first << " : " << children[i][j].second << std::endl;
        //         break;
        //     }
        // }

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
        }

        return std::tuple{children, all_branch, check_bottom};
    }

    // トポロジカルソートを呼び出す
    void call_topo() {
        std::cout << "--- call_topo ---" << std::endl;
        
        // auto [fork_info, blanch_num, check_bottom] = hash_trie_.return_partial_CP_info(restore_codes_); // trie.hppで計算
        // auto [fork_info, blanch_num, check_bottom] = return_partial_CP_info(); // vectorですべて計算
        // auto [fork_info, blanch_num, check_bottom] = return_partial_CP_info_using_map(); // map.hppで計算
        // std::cout << "size : " << fork_info.size() << std::endl;
        auto [children, blanch_num, check_bottom] = return_partial_CP_info_not_fork_info(); // fork_infoを使用しない方法
        std::cout << "size : " << children.size() << std::endl;

        // fork_infoの情報を元に、CP順を求め、新しい辞書に格納していく
        // data_reset();
        // hash_trie_.expand_tmp_table();
        // label_store_.expand_tmp_ptrs();
        std::string compare_str = "";
        // // std::cout << "aaa" << std::endl;
        // require_centroid_path_order_and_insert_dictionaly(fork_info, hash_trie_.get_root(), blanch_num, check_bottom, compare_str);
        require_centroid_path_order_and_insert_dictionaly_not_fork_info(children, hash_trie_.get_root(), blanch_num, check_bottom, compare_str);
        // // std::cout << "bbb" << std::endl;
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
