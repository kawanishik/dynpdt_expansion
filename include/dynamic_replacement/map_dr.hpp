#ifndef POPLAR_TRIE_MAP_HPP_DR
#define POPLAR_TRIE_MAP_HPP_DR

#include <array>
#include <iostream>
#include <queue>
#include <deque>
#include <algorithm>
#include <map>

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

    // 新しい辞書に登録するための関数
    template<class Map>
    void insert_new_dic(Map& new_map, uint64_t node_id) {
        std::string restore_key = restore_insert_string(node_id);
        int* ptr = new_map.update(restore_key);
        *ptr = 1;
    }

    struct parent_info {
        uint64_t parent;    // 親のノード番号
        uint64_t match;     // 文字列比較での失敗位置
        uint8_t c;          // 親からの遷移文字
    };

    // ノードのつながりと分岐数(葉の数)を調べるための関数
    void compute_node_connect_and_blanch_num(std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                             std::vector<uint64_t>& blanch_num_except_zero,
                                             std::vector<uint64_t>& cnt_leaf_per_node) {
        uint64_t table_size = hash_trie_.capa_size();
        std::vector<parent_info> parent(table_size);        // 親ノード、文字列比較での失敗位置、遷移文字を保存
        std::vector<uint64_t> partial_num(table_size, 0);   // 子の数を格納するための配列(自身の数も含む)
        cnt_leaf_per_node.resize(table_size, 0);            // それぞれのノードから繋がっている葉ノードの数をカウント
        blanch_num_except_zero.resize(table_size, 0);       // 特定のノード以下に何個のノードが存在するのか
        children.resize(table_size);                        // 子集合を保存
        std::vector<uint64_t> children_size(table_size);    // reserveで要素を確保するために使用

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
//        std::queue<uint64_t> que;
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
            if(match != 0) blanch_num_except_zero[p] += cnt_leaf_per_node[q];    // match=0は、親ノードからの分岐位置が同じことを示している
            else { // もう一つ上の親に足す(zero分岐のトライ上で考えるとそのようになるから)
                if(p != hash_trie_.get_root()) {
                    auto [pp, m , cc] = parent[p];
                    blanch_num_except_zero[pp] += cnt_leaf_per_node[q];
                }
            }
            partial_num[p]--;
            if(partial_num[p] == 1) { // 子供の処理がすべて終了すると追加
                que.emplace_back(p);
            }
        }
    }

    // 求めたノードの関係からCPD順を求める
    // 最終的には、ここで新しい辞書に対して、キーを追加する(まだ、未実装)
    template<class Map>
    void require_centroid_path_order_and_insert_dictionary(Map& new_map,
                                                           std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& children,
                                                           uint64_t node_id,
                                                           const std::vector<uint64_t>& blanch_num_except_zero,
                                                           const std::vector<uint64_t>& cnt_leaf) {
        // 一番下まで、たどり着いた時の処理
        if(children[node_id].size() == 0) return;

        std::vector<std::pair<uint64_t, uint64_t>> match_per_leaf_num;  // first:match, second:cnt_leaf[]
        std::map<uint64_t, uint64_t> start_pos;                         // children[node_id]内のmatchのスタート位置を保存
        uint64_t pre_match = UINT64_MAX;                                // 以前、比較したmatch(文字列比較失敗位置)を保存
        uint64_t zero_blanch_num = 0;                                   // ノード(node_id)のzero分岐の総数
        uint64_t children_size = children[node_id].size();              // 子ノードの数

        // 現在のchildrent[node_id]内は以下のようになっている
        // ex. ({3, 54}, {0, 4}, {13, 23}, {0, 5})
        // ソートして、下のように修正
        // ex. ({0, 4}, {0, 5}, {3, 54}, {13, 23})
        std::sort(children[node_id].begin(), children[node_id].end(), [] (auto l, auto r) {
            return l.first < r.first;
        });
        /* TODO: children の match を座標圧縮して処理する
        std::vector<uint64_t> matches(children[node_id].size());
        std::transform(children[node_id].begin(), children[node_id].end(), matches.begin(), [](auto& c) {return c.first;});
        matches.erase(std::unique(matches.begin(), matches.end()), matches.end());
         */

        // matchごとに、分岐後の葉の数をカウント
        for(uint64_t i=0; i < children_size; i++) {
            std::pair<uint64_t, uint64_t> child = children[node_id][i];
            if(child.first != pre_match) {
                match_per_leaf_num.emplace_back(child.first, cnt_leaf[child.second]);
                pre_match = child.first;
                start_pos[pre_match] = i + 1;
            } else {
                match_per_leaf_num.back().second += cnt_leaf[child.second];
            }
            if(child.first == 0) zero_blanch_num += cnt_leaf[child.second];
        }

        bool exist_zero_blanch = zero_blanch_num != 0;
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
                children_shelter.emplace_back(cnt_leaf[next_id], next_id);
            }
            std::sort(children_shelter.begin(), children_shelter.end(), [] (auto l, auto r) {
                return l.first > r.first;
            });
            for(auto& s : children_shelter) {
                require_centroid_path_order_and_insert_dictionary(new_map, children, s.second, blanch_num_except_zero, cnt_leaf);
                insert_new_dic(new_map, s.second); // 新しい辞書に登録（未実装）
            }
        }

        if(node_id == hash_trie_.get_root()) insert_new_dic(new_map, node_id); // 新しい辞書に登録（未実装）
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

    // 動的にいれかえるための関数
    void dynamic_replacement() {
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children;
        std::vector<uint64_t> blanch_num_except_zero;
        std::vector<uint64_t> cnt_leaf_per_node;
        compute_node_connect_and_blanch_num(children, blanch_num_except_zero, cnt_leaf_per_node);

        map_dr new_map(hash_trie_.capa_bits()+1);
        // std::cout << "now_map_capa_size : " << capa_size() <<std::endl;
        // std::cout << "new_map_capa_size : " << new_map.capa_size() << std::endl;
        require_centroid_path_order_and_insert_dictionary(new_map, children, hash_trie_.get_root(), blanch_num_except_zero, cnt_leaf_per_node);
        std::swap(*this, new_map); // 時間がかかるので、注意
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
};

}  // namespace poplar

#endif  // POPLAR_TRIE_MAP_HPP_DR
