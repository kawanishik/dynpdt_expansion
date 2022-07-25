#ifndef POPLAR_TRIE_MAP_SNR_HPP
#define POPLAR_TRIE_MAP_SNR_HPP

#include <array>
#include <iostream>

#include "../poplar/bit_tools.hpp"
#include "../poplar/exception.hpp"

static int expand_cnt = 0;

namespace poplar {

// This class implements an updatable associative array whose keys are strings.
// The data structure is based on a dynamic path-decomposed trie described in the following paper,
// - "Dynamic Path-Decomposed Tries" available at https://arxiv.org/abs/1906.06015.
template <typename Trie, typename NLM>
class map_SNR {
    static_assert(Trie::trie_type_id == NLM::trie_type_id);

  public:
    using this_type = map_SNR<Trie, NLM>;
    using trie_type = Trie;
    using value_type = typename NLM::value_type;

    static constexpr auto trie_type_id = Trie::trie_type_id;
    static constexpr uint32_t min_capa_bits = Trie::min_capa_bits;

  public:
    // Generic constructor.
    map_SNR() = default;

    // Class constructor. Initially allocates the hash table of length
    // 2**capa_bits.
    explicit map_SNR(uint32_t capa_bits, uint64_t lambda = 32) {
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
    ~map_SNR() = default;

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
                if(c == 0x00) restore_str.clear();  // 遷移文字が0x00(empty時に発生する処理)
                else restore_str = c + restore_str; // 遷移文字を追加

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

    // 新しい辞書に登録するための関数
    template<class Map>
    // void insert_new_dic(Map& new_map, uint64_t node_id, std::string prefix_str) {
    void insert_new_dic(Map& new_map, uint64_t node_id, uint64_t common_prefix_length, std::string& store_string) {
        // std::cout << "prefix_str : " << prefix_str << std::endl;
        std::string restore_key = "";
        if(common_prefix_length == 0) {
            restore_key = restore_insert_string(node_id);
            if(restore_key.size() == 0) return;
            // std::cout << "restore_key : " << restore_key << std::endl;
            store_string = restore_key;
            int* ptr = new_map.update(restore_key);
            *ptr = 1;
            return;
        } else {
            // restore_key = store_string.substr(0, common_prefix_length) + restore_node_string(node_id);
            // 遷移文字を取得し、文字列を復元する
            auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id); // 親ノードとsymbを取得
            auto [c, match] = restore_symb_(symb);                         // symbから、遷移に失敗した箇所とlabelを取得する

            if(common_prefix_length == 1) restore_key = c + restore_node_string(node_id);
            else restore_key = store_string.substr(0, common_prefix_length-1) + c + restore_node_string(node_id);
        }
        // std::string restore_key = restore_insert_string(node_id);
        if(restore_key.size() == 0) return;
        int* ptr = new_map.update(restore_key);
        *ptr = 1;
        // std::cout << "restore_key : " << restore_key << std::endl;
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

            if(c == 0x00) insert_string.clear();
            else insert_string = c + insert_string;

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

    // ノードのつながりと分岐数(葉の数)を調べるための関数
    void compute_node_connect_and_blanch_num(std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                            //  std::vector<uint64_t>& blanch_num_except_zero,
                                             std::vector<uint64_t>& cnt_leaf_per_node) {
        uint64_t table_size = hash_trie_.capa_size();
        cnt_leaf_per_node.resize(table_size, 0);            // それぞれのノードから繋がっている葉ノードの数をカウント
        // blanch_num_except_zero.resize(table_size, 0);    // 特定のノード以下に何個のノードが存在するのか
        children.resize(table_size);                        // 子集合を保存
        std::vector<parent_info> parent(table_size);        // 親情報を保存
        std::vector<uint64_t> partial_num(table_size, 0);   // 自身を含めた子の数を保存(子が存在しないとき、1)
        std::vector<uint64_t> children_size(table_size, 0); // 自身を含めない子ノードの数

        // O(n)で、親の位置、子の数(ノード番号も)を数える
        // 現在はダミーノードの数も計測してしまっている → 完了
        partial_num[hash_trie_.get_root()] += 1;
        for(uint64_t i=hash_trie_.get_root()+1; i < table_size; i++) {
            if(hash_trie_.is_use_table(i)) {                                                    // テーブルが使用されているとき
                if(label_store_.return_string_pointer(i) == nullptr) continue;                  // 対象としているノードがダミーノードの時
                auto [p, label] = hash_trie_.get_parent_and_symb(i);                            // 親と遷移情報の取得
                if(label == 255) continue;                                                      // ダミーノードをスキップ
                auto [c, match] = std::pair{uint8_t(restore_codes_[label % 256]), label/256};   // 遷移文字と分岐位置を取得
                while(label_store_.return_string_pointer(p) == nullptr) {                       // 親がダミーノードの時の処理(ダミーノードに文字列は保存されていないため)
                    auto [tmp_parent, tmp_label] = hash_trie_.get_parent_and_symb(p);
                    match += lambda_;
                    p = tmp_parent;
                }
                assert(p != nil_id);
                partial_num[p] += 1; // 子から親の数をカウントする
                children_size[p]++;
                partial_num[i] += 1; // 自身の数をカウントする
                parent[i].parent = p;
                parent[i].match = match;
                parent[i].c = c;
            }
        }
        for(uint64_t i=hash_trie_.get_root()+1; i < table_size; i++) {
          if (!hash_trie_.is_use_table(i)) continue;
          children[i].reserve(children_size[i]);
        }
        for(uint64_t i=hash_trie_.get_root()+1; i < table_size; i++) {
          if (!hash_trie_.is_use_table(i)) continue;
          auto p = parent[i].parent;
          auto match = parent[i].match;
          children[p].emplace_back(match, i);
        }

        // queueを使用して、一番下のものから処理していく
        // 初めに、一番下のノードのみを取得 (O(n))
        std::vector<uint64_t> que;
        que.reserve(hash_trie_.size());
        // 対象のデータを集めてくる
        for(uint64_t i=0; i < table_size; i++) {
            if(partial_num[i] == 1) {
                que.emplace_back(i);
            }
        }

        // queueに追加した順に処理（DynPDT上の一番下のノードから）
        for(uint64_t k = 0; k < que.size(); ++k) {                           // 追加された順に処理
            uint64_t q = que[k];
            auto [node_id, match, c] = parent[q];
            cnt_leaf_per_node[q] += 1;
            uint64_t p = node_id;
            assert(p != nil_id);
            cnt_leaf_per_node[p] += cnt_leaf_per_node[q];
            partial_num[p]--;
            if(partial_num[p] == 1) { // 子供の処理がすべて終了すると追加
                que.emplace_back(p);
            }
        }
    }

    // 後方累積和を使用して、CPD順を求めたもの
    template<class Map>
    bool require_centroid_path_order_and_insert_dictionary(Map& new_map,
                                                           std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children, // [node_id] {first:match, second:next_node_id}
                                                           uint64_t node_id,    // 処理するノード番号
                                                           const std::vector<uint64_t>& cnt_leaf_per_node, // 特定のノード以下に何個のノードがあるのか
                                                           uint64_t common_prefix_length,   // 共通の文字列の長さ
                                                           std::string& store_string,       // 引き継いだ文字列の長さ
                                                           int standard_height,             // 基準となる高さ(ここまでが入れ替えの高さ)
                                                           std::vector<int>& node_height,   // ノードごとの高さ(使用されていない部分は0)
                                                           /*std::vector<std::pair<bool, uint64_t>>& store_node,*/  // 基準値以上の高さを持っておく
                                                           std::vector<std::pair<uint64_t, uint64_t>>& store_leading_node) {    // 基準となる先頭ノードのみを保存する{node_id, common_prefix}
        // std::cout << "common_prefix_length : " << common_prefix_length << std::endl;

        // store_nodeを使用する場合
        // if(node_height[node_id] >= standard_height) {
        //     for(auto child : children[node_id]) {
        //         if(node_height[node_id] == standard_height) store_node.emplace_back(true, child.second);
        //         else store_node.emplace_back(false, child.second);
        //         require_centroid_path_order_and_insert_dictionary(new_map, children, child.second, cnt_leaf_per_node, common_prefix_length, store_string, standard_height, node_height, store_node);
        //     }
        //     if(node_height[node_id] == standard_height) return true;
        //     else return false;
        // }

        // store_leading_nodeを使用する場合
        if(node_height[node_id] >= standard_height) {   // 基準となる高さに到達した場合
            for(auto child : children[node_id]) {
                store_leading_node.emplace_back(child.second, common_prefix_length);
            }
            return true;
        }
        
        // 一番下まで、たどり着いた時の処理(これにあたるのは基準値よりも低い場合)
        if(children[node_id].size() == 0) return true;

        std::map<uint64_t, uint64_t> start_pos;     // children[node_id]内のmatchのスタート位置を保存

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
        std::vector<std::pair<uint64_t, uint64_t>> match_per_leaf_num(matches.size());  // first:match, second:cnt_leaf_per_node[]
        uint64_t pre_match = children[node_id][0].first;    // 以前、比較したmatch(文字列比較失敗位置)を保存
        start_pos[pre_match] = 1;
        uint64_t pos = 0;                                   // 位置を示す際に使用
        uint64_t zero_blanch_num = 0;                       // ノード(node_id)のzero分岐の総数
        uint64_t children_size = children[node_id].size();  // 子ノードの数
        uint64_t cumulative_sum = 0;                        // 累積和(0分岐を除く) 
        for(uint64_t i=0; i < children_size; i++) {
            auto [match, next_id] = children[node_id][i];
            if(pre_match != match) {
                pre_match = match;
                start_pos[match] = i + 1;
                pos++;
            }
            // uint64_t pos = std::lower_bound(matches.begin(), matches.end(), match) - matches.begin();
            match_per_leaf_num[pos].first = match;
            match_per_leaf_num[pos].second += cnt_leaf_per_node[next_id];
            if(match == 0) zero_blanch_num += cnt_leaf_per_node[next_id];  // 0分岐の和
            else cumulative_sum += cnt_leaf_per_node[next_id];             // 0分岐以外の累積和
        }

        bool exist_zero_blanch = zero_blanch_num != 0;
        // // 後方累積和を計算する
        uint64_t size;
        std::vector<uint64_t> backward_cumulative_num;  // 後方累積和の値の保存場所
        if(exist_zero_blanch) { // zero分岐が存在する場合、zero分岐以外の和を取る必要があるため
            size = matches.size() - 1;
            backward_cumulative_num.resize(size, 0);
            for(int64_t i=size-1; i >= 0; i--) {
                if(i == (size-1)) backward_cumulative_num[i] = match_per_leaf_num[i+1].second;
                else backward_cumulative_num[i] = backward_cumulative_num[i+1] + match_per_leaf_num[i+1].second;
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
        std::queue<uint64_t> front;
        std::stack<uint64_t> back;

        if(exist_zero_blanch) { // zero分岐を最初、または最後に処理するのか
            if(zero_blanch_num <= cumulative_sum) back.push(0);
            else front.push(0);
        }

        if(size > 0) {
            for(uint64_t i=0; i < size; i++) {
                pos = exist_zero_blanch ? (i+1) : (i);
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
        std::string compared_string = ""; // 帰りがけで取得した文字列を保存しておく文字列
        while(!front.empty()) {
            uint64_t num = front.front();
            front.pop();
            pre_match = match_per_leaf_num[num].first;
            std::vector<std::pair<uint64_t, uint64_t>> children_shelter;    // first:cnt_leaf_per_node[], second:next_id
            for(uint64_t i=start_pos[pre_match]-1; i < children_size; i++) {
                if(pre_match != children[node_id][i].first) break;
                uint64_t next_id = children[node_id][i].second;
                children_shelter.emplace_back(cnt_leaf_per_node[next_id], next_id);
            }

            // 特定の分岐内の内、分岐後の葉の数が多いnext_idを取得
            auto max_pos = std::max_element(children_shelter.begin(), children_shelter.end());
            uint64_t next_common_prefix_length = common_prefix_length+pre_match+1; // 行き掛けで渡す共通のprefixの長さを伝搬する
            // bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, max_pos->second, cnt_leaf_per_node, next_common_prefix_length, store_string, standard_height, node_height, store_node);
            bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, max_pos->second, cnt_leaf_per_node, next_common_prefix_length, store_string, standard_height, node_height, store_leading_node);
            // ここで、original_stringからcompare_prefix_length分のみを抜き出し、insert_new_dicの関数の引数とする
            if(is_leaf_node) {
                insert_new_dic(new_map, max_pos->second, 0, store_string);
            } else {
                insert_new_dic(new_map, max_pos->second, next_common_prefix_length, store_string);
            }
            // insert_new_dic(new_map, max_pos->second, "");
            for(auto& s : children_shelter) {
                if(s.second == max_pos->second) continue;
                // bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string, standard_height, node_height, store_node);
                bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string, standard_height, node_height, store_leading_node);
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
            pre_match = match_per_leaf_num[num].first;
            std::vector<std::pair<uint64_t, uint64_t>> children_shelter;                // first:cnt_leaf_per_node[], second:next_id
            for(uint64_t i=start_pos[pre_match]-1; i < children_size; i++) {
                if(pre_match != children[node_id][i].first) break;
                uint64_t next_id = children[node_id][i].second;
                children_shelter.emplace_back(cnt_leaf_per_node[next_id], next_id);
            }

            // 特定の分岐内の内、分岐後の葉の数が多いnext_idを取得
            auto max_pos = std::max_element(children_shelter.begin(), children_shelter.end());
            uint64_t next_common_prefix_length = common_prefix_length+pre_match+1;
            // bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, max_pos->second, cnt_leaf_per_node, next_common_prefix_length, store_string, standard_height, node_height, store_node);
            bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, max_pos->second, cnt_leaf_per_node, next_common_prefix_length, store_string, standard_height, node_height, store_leading_node);
            if(is_leaf_node) {
                insert_new_dic(new_map, max_pos->second, 0, store_string);
            } else {
                insert_new_dic(new_map, max_pos->second, next_common_prefix_length, store_string);
            }
            for(auto& s : children_shelter) {
                if(s.second == max_pos->second) continue;
                // bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string, standard_height, node_height, store_node);
                bool is_leaf_node = require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, cnt_leaf_per_node, next_common_prefix_length, store_string, standard_height, node_height, store_leading_node);
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

    // store_nodeに保存されているノードを処理するための関数
    // 現在は、一つずつ文字列を復元することによって、代用している
    template<class Map>
    void process_store_node(Map& new_map, std::vector<std::pair<bool, uint64_t>>& store_node) {
        std::vector<std::string> restore_keys;
        for(auto node : store_node)  {
            std::string restore_key = restore_insert_string(node.second);
            if(restore_key.size() == 0) continue;
            restore_keys.emplace_back(restore_key);
            int* ptr = new_map.update(restore_key);
            *ptr = 1;
        }
        write_file(restore_keys);
    }

    // store_leading_nodeを処理するための関数
    template<class Map>
    void process_store_leading_node(Map& new_map,
                                    std::vector<std::pair<uint64_t, uint64_t>>& store_leading_node, // {node_id, common_prefix_length}
                                    std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children) {    // [node_id] {match, node_id}
        std::string restore_key = "";
        for(auto [node_id, common_prefix_length] : store_leading_node) {
            subtree_add_new_map(new_map, node_id, restore_key, common_prefix_length, children);
        }
    }

    // 与えられたノードに対して、再帰関数を使用して、新しいmapに追加する
    template<class Map>
    void subtree_add_new_map(Map& new_map,
                             uint64_t node_id,
                             std::string restore_key,
                             uint64_t common_prefix_length,
                             std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children) {
        insert_new_dic(new_map, node_id, 0, restore_key);
        for(auto child : children[node_id]) {
            subtree_add_new_map(new_map, child.second, restore_key, common_prefix_length+child.first, children);
        }
    }

    int calc_height(std::vector<int>& node_height, std::vector<uint64_t>& each_node_height_sum) {
        uint64_t table_size = hash_trie_.capa_size();
        node_height.resize(table_size, 0);
        node_height[hash_trie_.get_root()] = 1;
        each_node_height_sum.resize(10, 0);
        each_node_height_sum[1] = 1;
        // std::vector<bool> done_flags(table_size, false);    // 処理が完了したノードかの判定
        bit_vector done_flags(capa_size());
        done_flags.set(hash_trie_.get_root());
        int max_height = 1;

        std::vector<std::pair<uint64_t, uint64_t>> path;
        path.reserve(256);

        for(uint64_t i=2; i < table_size; ++i) {
            if(done_flags[i] || hash_trie_.is_use_table(i) == false) {   // 処理が完了している
                continue;
            }

            path.clear();
            uint64_t node_id = i;

            do {
                auto [parent, label] = hash_trie_.get_parent_and_symb(node_id);
                assert(parent != nil_id);
                path.emplace_back(std::make_pair(node_id, label));
                node_id = parent;
            } while(!done_flags[node_id]);

            int now_height = node_height[node_id];

            for(auto rit = std::rbegin(path); rit != std::rend(path); ++rit) {
                now_height++;
                node_id = rit->first;
                node_height[node_id] = now_height;
                if(now_height >= int(each_node_height_sum.size())) each_node_height_sum.resize(now_height+1, 0);
                each_node_height_sum[now_height] += 1;
                done_flags.set(node_id);
            }

            max_height = std::max(max_height, now_height);
        }

        // std::cout << "max_height : " << max_height << std::endl;
        return max_height;
    }

    void calc_height_sum(const std::vector<int>& height) {
        expand_cnt += 1;
        std::vector<int> height_sum;
        for(auto h : height) {
            if(h >= int(height_sum.size())) height_sum.resize(h+1, 0);
            height_sum[h] += 1;
        }
        write_file(height_sum);
    }

    void check_store_node(std::vector<std::pair<bool, uint64_t>>& store_node,
                          std::vector<int>& node_height,
                          int standard_height,
                          std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children) {
        std::cout << "--- check_store_node ---" << std::endl;
        std::cout << "size : " << size() << std::endl;  // キーの数

        // 想定している数がstore_nodeに格納されているのか
        // int sum = 0;
        // // for(int i=0; i < int(node_height.size()); i++) {
        // //     if(node_height[i] > standard_height) sum++;
        // // }

        // std::map<uint64_t, int> mp;
        // for(auto node : store_node) mp[node.second] = 1;
        // for(int i=0; i < int(node_height.size()); i++) {
        //     if(node_height[i] > standard_height) {
        //         if(mp[i] == 0) {
        //             sum++;
        //             std::cout << i << std::endl;
        //             std::cout << "child_num1 : " << children[i].size() << std::endl;
        //             auto [parent, label] = hash_trie_.get_parent_and_symb(i);
        //             // std::cout << "parent : " << parent << std::endl;
        //             // std::cout << "label : " << label << std::endl;
        //             // std::cout << "child_num2 : " << children[parent].size() << std::endl;
        //             std::cout << parent << ", " << label << ", " << children[parent].size() << std::endl;
        //             std::string restore_key = restore_insert_string(i);
        //             if(restore_key.size() == 0) std::cout << "dummy" << std::endl;
        //             else std::cout << "not dummy" << std::endl;
        //         }
        //     }
        // }
        // std::cout << store_node.size() << " , " << sum << std::endl;

        // 全てのノードが基準の高さよりも高いことの確認 → 特に問題なし
        // int failed_node_num = 0;
        // for(auto node : store_node) {
        //     if(node_height[node.second] <= standard_height) {
        //         failed_node_num++;
        //     }
        // }
        // std::cout << "failed_node_num : " << failed_node_num << std::endl;

        // boolがtrueの部分から次のtrueの部分まで、一つの繋がりになっているのかどうか(分岐があるので注意が必要)

    }

    void check_store_leading_node(std::vector<std::pair<uint64_t, uint64_t>>& store_leading_node,   // {node_id, common_prefix}
                                  std::vector<int>& node_height,
                                  int standard_height,
                                  std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children) {  // [node_id] {first:match, second:node_id}
        // store_leading_nodeとstandard_height+1の数が同じであることの確認
        uint64_t table_size = hash_trie_.capa_size();
        int size = store_leading_node.size();
        int sum = 0;
        for(uint64_t i=1; i < table_size; i++) {
            if(node_height[i] == (standard_height+1)) sum++;
        }
        if(size == sum) {
            std::cout << "same" << std::endl;
        } else {    // children内でダミーノードはノードとして、認められていないから
            std::cout << "different" << std::endl;
        }
        std::cout << size << ", " << sum << std::endl;

        int failed_node = 0;
        for(auto node : store_leading_node) {
            if(node_height[node.first] != (standard_height+1)) failed_node++;
        }
        std::cout << failed_node << std::endl;
        
        // 格納されているノードの高さがすべてstandard_height+1であることの確認→ダミーノードの関係上、そうとは限らない
    }

    void write_file(std::vector<int>& height_sum) {
        std::ofstream of;
        std::string filename = "../../../memo/enwiki/";
        // std::string filename = "../../../memo/AOL/";
        // std::string filename = "../../../memo/GeoNames/";
        // std::string filename = "../../../memo/in-2004/";
        filename += std::to_string(expand_cnt);
        filename += ".txt";
        of.open(filename, std::ios::out);
        // for(auto b : tmp_keys) {
        //     of << b << std::endl;
        // }
        for(int i=0; i < int(height_sum.size()); i++) {
            of << i << ":" << height_sum[i] << std::endl;
        }
        of.close();
    }

    void write_file(std::vector<std::string>& restore_keys) {
        std::ofstream of;
        std::string filename = "../../../memo/enwiki/restore.txt";
        std::sort(restore_keys.begin(), restore_keys.end());
        of.open(filename, std::ios::out);
        for(int i=0; i < int(restore_keys.size()); i++) {
            of << restore_keys[i] << std::endl;
        }
    } 

    void dynamic_replacement() {
        std::vector<int> node_height;
        std::vector<uint64_t> each_node_height_sum;
        int max_height = calc_height(node_height, each_node_height_sum);  // それぞれのノードの高さと最大のノードの高さを取得
        // calc_height_sum(node_height);
        uint64_t standard_height = 0;           // 基準となる木の高さ(ここまでは入れ替えの対象)
        uint64_t node_sum = hash_trie_.size();  // ノードの合計数
        double ratio_sum = 0.0;                 // 比率の合計

        for(int i=1; i < int(each_node_height_sum.size()); i++) {
            // std::cout << i << " : " << each_node_height_sum[i] << std::endl;
            double ratio = (long double)(each_node_height_sum[i]) / (long double)(node_sum);
            if((ratio_sum+ratio) > 0.5) {
                standard_height = i-1;
                break;
            }
            ratio_sum += ratio;
        }

        std::cout << "max_height : " << max_height << std::endl;
        std::cout << "standard_height : " << standard_height << std::endl;
        // std::cout << "ratio_sum : " << ratio_sum << std::endl;

        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children;   // 子ノードの集合
        std::vector<uint64_t> cnt_leaf_per_node;
        compute_node_connect_and_blanch_num(children, cnt_leaf_per_node);

        map_SNR new_map(hash_trie_.capa_bits());

        // std::vector<std::pair<bool, uint64_t>> store_node;  // standard_height以上の高さの木を処理する順に保存しておく
        std::vector<std::pair<uint64_t, uint64_t>> store_leading_node;  // 基準のノードの先頭を保持する
        std::string store_string = "";
        require_centroid_path_order_and_insert_dictionary(new_map, children, hash_trie_.get_root(), cnt_leaf_per_node, 0, store_string, standard_height, node_height, /*store_node,*/ store_leading_node);

        // store_nodeの確認(想定通りの動きをしているのか)
        // check_store_node(store_node, node_height, standard_height, children);

        // store_leading_nodeの確認
        check_store_leading_node(store_leading_node, node_height, standard_height, children);

        // store_nodeを辞書に登録する
        // process_store_node(new_map, store_node);

        // store_leading_nodeを新しい辞書に登録する
        process_store_leading_node(new_map, store_leading_node, children);

        std::swap(*this, new_map);
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
        show_stat(os, indent, "name", "map_SNR");
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

    map_SNR(const map_SNR&) = delete;
    map_SNR& operator=(const map_SNR&) = delete;

    map_SNR(map_SNR&&) noexcept = default;
    map_SNR& operator=(map_SNR&&) noexcept = default;

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
            // std::cout << "--- expand ---" << std::endl;
            // ここですべての文字列を復元する
            // std::vector<std::string> all_keys = all_key_restore_simple(); // 文字列の復元

            // ノードがどの高さにいるのかを調べる
            // std::vector<int> node_height;
            // std::vector<uint64_t> each_node_height_sum;
            // int max_height = calc_height(node_height, each_node_height_sum);  // それぞれのノードの高さと最大のノードの高さを取得
            // // calc_height_sum(node_height);
            // uint64_t standard_height = 0;           // 基準となる木の高さ
            // uint64_t node_sum = hash_trie_.size();  // ノードの合計数
            // double ratio_sum = 0.0;                 // 比率の合計

            // for(int i=1; i < int(each_node_height_sum.size()); i++) {
            //     // std::cout << i << " : " << each_node_height_sum[i] << std::endl;
            //     double ratio = (long double)(each_node_height_sum[i]) / (long double)(node_sum);
            //     if((ratio_sum+ratio) > 0.5) {
            //         standard_height = i-1;
            //         break;
            //     }
            //     ratio_sum += ratio;
            // }

            // std::cout << "max_height : " << max_height << std::endl;
            // std::cout << "standard_height : " << standard_height << std::endl;
            // // std::cout << "ratio_sum : " << ratio_sum << std::endl;

            // std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children;   // 子ノードの集合
            // std::vector<uint64_t> cnt_leaf_per_node;
            // compute_node_connect_and_blanch_num(children, cnt_leaf_per_node);

            auto node_map = hash_trie_.expand();
            node_id = node_map[node_id];
            label_store_.expand(node_map);
        }
    }
};

}  // namespace poplar

#endif  // POPLAR_TRIE_MAP_SNR_HPP
