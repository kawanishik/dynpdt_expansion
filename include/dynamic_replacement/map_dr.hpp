#ifndef POPLAR_TRIE_MAP_HPP_DR
#define POPLAR_TRIE_MAP_HPP_DR

#include <array>
#include <iostream>
#include <queue>
#include <deque>
#include <algorithm>
#include <map>
#include <stack>
#include <string_view>

#include "../poplar/bit_tools.hpp"
#include "../poplar/exception.hpp"

namespace poplar {

static std::map<int, int> cnt_hash;
static uint64_t all_cnt = 0;

// This class implements an updatable associative array whose keys are strings.
// The data structure is based on a dynamic path-decomposed trie described in the following paper,
// - "Dynamic Path-Decomposed Tries" available at https://arxiv.org/abs/1906.06015.
template <typename Trie, typename NLM>
class map_dr {
    static_assert(Trie::trie_type_id == NLM::trie_type_id);

  public:
    using this_type = map_dr<Trie, NLM>;
    using trie_type = Trie;
    using value_type = typename NLM::value_type;

    static constexpr auto trie_type_id = Trie::trie_type_id;
    static constexpr uint32_t min_capa_bits = Trie::min_capa_bits;

  public:
    // Generic constructor.
    map_dr() = default;

    // Class constructor. Initially allocates the hash table of length
    // 2**capa_bits.
    explicit map_dr(uint32_t capa_bits, uint64_t lambda = 32) {
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
    ~map_dr() = default;

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

    // 新しい辞書に登録するための関数
    template<class Map>
    // void insert_new_dic(Map& new_map, uint64_t node_id, std::string prefix_str) {
    void insert_new_dic(Map& new_map, uint64_t node_id, uint64_t common_prefix_length, std::string& store_string) {
        // std::cout << "prefix_str : " << prefix_str << std::endl;
        std::string restore_key = "";
        if(common_prefix_length == 0) {
            restore_key = restore_insert_string(node_id);
            // std::cout << "restore_key : " << restore_key << std::endl;
            store_string = restore_key;
            int* ptr = new_map.update(restore_key);
            *ptr = 1;
            return;
        } else {
            restore_key = store_string.substr(0, common_prefix_length) + restore_node_string(node_id);
        }
        // std::string restore_key = restore_insert_string(node_id);
        if(restore_key.size() == 0) return;
        int* ptr = new_map.update(restore_key);
        *ptr = 1;
        // std::cout << "restore_key : " << restore_key << std::endl;
    }

    template<class Map>
    void test_insert_new_dic(Map& new_map,
                             uint64_t node_id,
                             uint64_t common_prefix_length,
                             std::string& store_string,
                             std::vector<std::string>& restore_keys_order,
                             std::vector<std::pair<std::string, uint64_t>>& not_leaf_node_compare_string) {
        // std::cout << "prefix_str : " << prefix_str << std::endl;
        std::string restore_key = "";
        if(common_prefix_length == 0) {
            restore_key = restore_insert_string(node_id);
            // std::cout << "restore_key : " << restore_key << std::endl;
            store_string = restore_key;
            int* ptr = new_map.update(restore_key);
            *ptr = 1;
            restore_keys_order.emplace_back(restore_key);
            return;
            // return restore_key;
        }
        else restore_key = store_string.substr(0, common_prefix_length) + restore_node_string(node_id);
        // std::string restore_key = restore_insert_string(node_id);
        if(restore_key.size() == 0) return;
        int* ptr = new_map.update(restore_key);
        *ptr = 1;
        restore_keys_order.emplace_back(restore_key);
        not_leaf_node_compare_string.emplace_back(restore_key, common_prefix_length);
        // std::cout << "store_string: " << store_string << std::endl;
        // std::cout << "restore_key : " << restore_key << std::endl;
        // return "";
    }

    // 親情報などを保存するときに使用
    struct parent_info {
        uint64_t parent;    // 親のノード番号
        uint64_t match;     // 文字列比較での失敗位置
        uint8_t c;          // 親からの遷移文字
    };

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

    // 一時変数をまとめたもの
    struct temporary_set_of_variable_require {
        temporary_set_of_variable_require(uint64_t pm, uint64_t zbn, uint64_t cs, uint64_t p, uint64_t cus) :
                        pre_match{pm}, zero_blanch_num{zbn}, children_size{cs}, pos{p}, cumulative_sum{cus} {}

        uint64_t pre_match;                 // 以前、比較したmatch(文字列比較失敗位置)を保存
        uint64_t zero_blanch_num;           // ノード(node_id)のzero分岐の総数
        uint64_t children_size;             // 子ノードの数
        uint64_t pos;                       // 位置を示す際に使用
        uint64_t cumulative_sum;            // 累積和(0分岐を除く) 
        uint64_t match_per_leaf_num_size;   // vector match_per_leaf_numのサイズ
        std::vector<std::pair<uint64_t, uint64_t>> match_per_leaf_num;  // first:match, second:cnt_leaf_per_node[]
    };


    // require_centroid_path_order_and_insert_dictionaryで使用する一時変数をまとめたもの
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
                insert_new_dic(new_map, max_pos->second, 0, store_string);
            } else {
                insert_new_dic(new_map, max_pos->second, next_common_prefix_length, store_string);
            }
            // insert_new_dic(new_map, max_pos->second, "");
            for(auto& s : children_shelter) {
                if(s.second == max_pos->second) continue;
                bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string);
                if(is_leaf_node) {
                    insert_new_dic(new_map, s.second, 0, store_string);
                } else {
                    insert_new_dic(new_map, s.second, next_common_prefix_length, store_string);
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
                insert_new_dic(new_map, max_pos->second, 0, store_string);
            } else {
                insert_new_dic(new_map, max_pos->second, next_common_prefix_length, store_string);
            }
            for(auto& s : children_shelter) {
                if(s.second == max_pos->second) continue;
                bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string);
                if(is_leaf_node) {
                    insert_new_dic(new_map, s.second, 0, store_string);
                } else {
                    insert_new_dic(new_map, s.second, next_common_prefix_length, store_string);
                }
            }
        }

        if(node_id == hash_trie_.get_root()) insert_new_dic(new_map, node_id, 0, store_string); // 新しい辞書に登録（未実装）

        return false;
    }

    template<class Map>
    bool test_CPD_order_and_insert_new_dic(
                                            Map& new_map,
                                            std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                            uint64_t node_id,
                                            const std::vector<uint64_t>& cnt_leaf_per_node,
                                            uint64_t common_prefix_length,
                                            std::string& store_string,
                                            std::vector<std::string>& restore_keys_order,
                                            std::vector<std::pair<std::string, uint64_t>>& not_leaf_node_compare_string) {
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
            bool is_leaf_node = test_CPD_order_and_insert_new_dic(new_map, children, max_pos->second, cnt_leaf_per_node, next_common_prefix_length, store_string, restore_keys_order, not_leaf_node_compare_string);
            // ここで、original_stringからcompare_prefix_length分のみを抜き出し、insert_new_dicの関数の引数とする
            if(is_leaf_node) {
                test_insert_new_dic(new_map, max_pos->second, 0, store_string, restore_keys_order, not_leaf_node_compare_string);
            } else {
                test_insert_new_dic(new_map, max_pos->second, next_common_prefix_length, store_string, restore_keys_order, not_leaf_node_compare_string);
            }
            for(auto& s : children_shelter) {
                if(s.second == max_pos->second) continue;
                bool is_leaf_node = test_CPD_order_and_insert_new_dic(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string, restore_keys_order, not_leaf_node_compare_string);
                if(is_leaf_node) {
                    test_insert_new_dic(new_map, s.second, 0, store_string, restore_keys_order, not_leaf_node_compare_string);
                } else {
                    test_insert_new_dic(new_map, s.second, next_common_prefix_length, store_string, restore_keys_order, not_leaf_node_compare_string);
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
            bool is_leaf_node = test_CPD_order_and_insert_new_dic(new_map, children, max_pos->second, cnt_leaf_per_node, next_common_prefix_length, store_string, restore_keys_order, not_leaf_node_compare_string);
            if(is_leaf_node) {
                test_insert_new_dic(new_map, max_pos->second, 0, store_string, restore_keys_order, not_leaf_node_compare_string);
            } else {
                test_insert_new_dic(new_map, max_pos->second, next_common_prefix_length, store_string, restore_keys_order, not_leaf_node_compare_string);
            }
            for(auto& s : children_shelter) {
                if(s.second == max_pos->second) continue;
                bool is_leaf_node = test_CPD_order_and_insert_new_dic(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string, restore_keys_order, not_leaf_node_compare_string);
                if(is_leaf_node) {
                    test_insert_new_dic(new_map, s.second, 0, store_string, restore_keys_order, not_leaf_node_compare_string);
                } else {
                    test_insert_new_dic(new_map, s.second, next_common_prefix_length, store_string, restore_keys_order, not_leaf_node_compare_string);
                }
            }
        }

        if(node_id == hash_trie_.get_root()) test_insert_new_dic(new_map, node_id, 0, store_string, restore_keys_order, not_leaf_node_compare_string); // 新しい辞書に登録（未実装）

        return false;
    }

    // 再帰関数を使用せずにCPD順を求めるための関数
    // cpd順を追加する順序の最後から順に求めている
    // 再帰を使用していた時とは異なり、cpd順を求めながら、新しい辞書に登録することはできない
    // template <class Map>
    std::stack<uint64_t> require_centroid_path_order_not_using_recursion(//Map& new_map,
                                                        std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                                        const std::vector<uint64_t>& cnt_leaf_per_node) {
                                                        // uint64_t common_prefix_length,
                                                        // std::string& store_string,
                                                        // std::vector<std::string>& restore_keys_order,
                                                        // std::vector<std::pair<std::string, uint64_t>>& not_leaf_node_compare_string) {
        uint64_t node_id = hash_trie_.get_root();   
        std::stack<uint64_t> process_ord;   // 処理するノード順
        process_ord.push(node_id);
        std::stack<uint64_t> cpd_ord;       // 求めたCPD順

        while(!process_ord.empty()) {
            node_id = process_ord.top();
            process_ord.pop();
            cpd_ord.push(node_id);      // cpd順を下から順に求めるため
            // std::cout << node_id << " : " << node_id << std::endl;

            uint64_t children_size = children[node_id].size();
            if(children_size == 0) continue;

            std::map<uint64_t, uint64_t> start_pos; // children[node_id]内のmatchのスタート位置を保存

            std::sort(children[node_id].begin(), children[node_id].end(), [] (auto l, auto r) {  // 子ノードの中で、分岐順でソート
                return l.first < r.first;
            });

            // 座標圧縮
            std::vector<uint64_t> matches(children[node_id].size());
            std::transform(children[node_id].begin(), children[node_id].end(), matches.begin(), [](auto& c) {return c.first;});
            matches.erase(std::unique(matches.begin(), matches.end()), matches.end());

            // 
            std::vector<std::pair<uint64_t, uint64_t>> match_per_leaf_num(matches.size());   // 分岐位置ごとの子の数の合計
            uint64_t pre_match = children[node_id][0].first;            // どこに保存するのか
            start_pos[pre_match] = 1;
            uint64_t pos = 0;
            uint64_t zero_blanch_num = 0;
            uint64_t cumulative_sum = 0;
            for(uint64_t i=0; i < children_size; i++) {
                auto [match, next_id] = children[node_id][i];
                if(pre_match != match) {
                    pre_match = match;
                    start_pos[match] = i+1;
                    pos++;
                }
                match_per_leaf_num[pos].first = match;
                match_per_leaf_num[pos].second += cnt_leaf_per_node[next_id];
                if(match == 0) zero_blanch_num += cnt_leaf_per_node[next_id];
                else cumulative_sum += cnt_leaf_per_node[next_id];
            }

            bool exist_zero_blanch = zero_blanch_num != 0;
            // 後方累積和
            uint64_t size;
            std::vector<uint64_t> backward_cumulative_num;
            if(exist_zero_blanch) { // Zero分岐が存在する場合
                size = matches.size() - 1;
                backward_cumulative_num.resize(size, 0);
                for(int64_t i=size-1; i >= 0; i--) {
                    if(i == (size-1)) backward_cumulative_num[i] = match_per_leaf_num[i].second;
                    else backward_cumulative_num[i] = backward_cumulative_num[i+1] + match_per_leaf_num[i].second;
                }
            } else {
                size = matches.size();
                backward_cumulative_num.resize(size, 0);
                for(int64_t i=size-1; i >= 0; i--) {
                    if(i == (size-1)) backward_cumulative_num[i] = match_per_leaf_num[i].second;
                    else backward_cumulative_num[i] = backward_cumulative_num[i+1] + match_per_leaf_num[i].second;
                }
            }

            // 処理する順番を決める
            // 元の処理だと、優先度が高い順に求めていた
            // ここでは、優先度が低い順に求めることによって、処理する順番が遅い順にスタックに詰めることによって、処理を実現する
            // 下のstackとqueueに格納するのは、matchなので、注意
            std::queue<uint64_t> front; // 辞書に追加する際に、対象とするノードが早い順にstackに格納するので、このままでよい
            std::stack<uint64_t> back;

            if(exist_zero_blanch) {
                if(zero_blanch_num <= cumulative_sum) back.push(0);
                else front.push(0);
            }

            if(size > 0) {
                for(uint64_t i=0; i < size; i++) {
                    uint64_t pos = exist_zero_blanch ? (i+1) : (i);
                    if(i == (size-1)) {
                        front.push(pos);
                    } else {
                        if(match_per_leaf_num[pos].second > backward_cumulative_num[i+1]) {
                            front.push(pos);
                        } else {
                            back.push(pos);
                        }
                    }
                }
            }

            while(!front.empty()) {
                uint64_t num = front.front();
                front.pop();
                pre_match = match_per_leaf_num[num].first;
                std::vector<std::pair<uint64_t, uint64_t>> children_shelter;    // first:cnt_per_node[], second:next_id
                for(uint64_t i=start_pos[pre_match]-1; i < children_size; i++) {
                    if(pre_match != children[node_id][i].first) break;
                    uint64_t next_id = children[node_id][i].second;
                    children_shelter.emplace_back(cnt_leaf_per_node[next_id], next_id);
                }

                // 特定の分岐内の内、分岐後の葉の数が多い順にソート
                std::sort(children_shelter.begin(), children_shelter.end(), [] (auto l, auto r) {
                    return l.first > r.first;
                });

                // 求めた順にprocess_ordに格納していく
                for(auto child : children_shelter) process_ord.push(child.second);
            }
            while(!back.empty()) {
                uint64_t num = back.top();
                back.pop();
                pre_match = match_per_leaf_num[num].first;
                std::vector<std::pair<uint64_t, uint64_t>> children_shelter;    // first:cnt_per_node[], second:next_id
                for(uint64_t i=start_pos[pre_match]-1; i < children_size; i++) {
                    if(pre_match != children[node_id][i].first) break;
                    uint64_t next_id = children[node_id][i].second;
                    children_shelter.emplace_back(cnt_leaf_per_node[next_id], next_id);
                }
                // 特定の分岐内の内、分岐後の葉の数が多い順にソート
                std::sort(children_shelter.begin(), children_shelter.end(), [] (auto l, auto r) {
                    return l.first > r.first;
                });
                // 求めた順にprocess_ordに格納していく
                for(auto child : children_shelter) process_ord.push(child.second);
            }
        }
        return cpd_ord;
    }

    template <class Map>
    void restore_keys_and_insert_dictionary(Map& new_map, std::stack<uint64_t>& cpd_ord) {
        while(!cpd_ord.empty()) {
            uint64_t node_id = cpd_ord.top();
            cpd_ord.pop();
            // 新しい辞書に登録
            int* ptr = new_map.update(restore_insert_string(node_id));
            *ptr = 1;
        }
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

    // 動的にいれかえるための関数
    value_type* dynamic_replacement(char_range& key) {
    // void dynamic_replacement() {
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children;   // 子ノードの集合
        // std::vector<uint64_t> blanch_num_except_zero;                    // 0分岐を除く累積和
        std::vector<uint64_t> cnt_leaf_per_node;
        compute_node_connect_and_blanch_num(children, cnt_leaf_per_node);

        // map_dr new_map(hash_trie_.capa_bits());
        map_dr new_map(hash_trie_.capa_bits()+1);
        // std::cout << "now_map_capa_size : " << capa_size() <<std::endl;
        // std::cout << "new_map_capa_size : " << new_map.capa_size() << std::endl;
        
        // 再帰関数を使用する場合
        // std::string store_string = "";
        // require_centroid_path_order_and_insert_dictionary(new_map, children, hash_trie_.get_root(), cnt_leaf_per_node, 0, store_string);
        
        // 再帰関数を使用しない場合
        std::stack<uint64_t> cpd_ord = require_centroid_path_order_not_using_recursion(children, cnt_leaf_per_node);
        restore_keys_and_insert_dictionary(new_map, cpd_ord);

        std::swap(*this, new_map); // 時間がかかるので、注意

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
        show_stat(os, indent, "name", "map_dr");
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

    map_dr(const map_dr&) = delete;
    map_dr& operator=(const map_dr&) = delete;

    map_dr(map_dr&&) noexcept = default;
    map_dr& operator=(map_dr&&) noexcept = default;

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

#endif  // POPLAR_TRIE_MAP_HPP_DR
