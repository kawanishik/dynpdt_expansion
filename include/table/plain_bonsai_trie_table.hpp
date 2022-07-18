#ifndef POPLAR_TRIE_PLAIN_BONSAI_TRIE_TABLE_HPP
#define POPLAR_TRIE_PLAIN_BONSAI_TRIE_TABLE_HPP

#include "../poplar/bit_tools.hpp"
#include "../poplar/bit_vector.hpp"
#include "../poplar/compact_vector.hpp"
#include "../poplar/hash.hpp"

namespace poplar {

template <uint32_t MaxFactor = 90, typename Hasher = hash::vigna_hasher>
class plain_bonsai_trie_table {
  private:
    static_assert(0 < MaxFactor and MaxFactor < 100);

    struct table_info {
        uint64_t node_id;
        std::vector<uint64_t> table;
    };

  public:
    static constexpr uint64_t nil_id = UINT64_MAX;
    static constexpr uint32_t min_capa_bits = 16; // 5

    static constexpr auto trie_type_id = trie_type_ids::BONSAI_TRIE;

  public:
    plain_bonsai_trie_table() = default;

    plain_bonsai_trie_table(uint32_t capa_bits, uint32_t symb_bits) {
        capa_size_ = size_p2{std::max(min_capa_bits, capa_bits)};
        symb_size_ = size_p2{symb_bits};
        max_size_ = static_cast<uint64_t>(capa_size_.size() * MaxFactor / 100.0);
        table_ = compact_vector{capa_size_.size(), capa_size_.bits() + symb_size_.bits()};
        node_use_table.resize(capa_size_.size(), false);
    }

    ~plain_bonsai_trie_table() = default;

    uint64_t get_root() const {
        assert(size_ != 0);
        return 1;
    }

    void add_root() {
        assert(size_ == 0);
        size_ = 1;
    }

    uint64_t find_child(uint64_t node_id, uint64_t symb) const {
        assert(node_id < capa_size_.size());
        assert(symb < symb_size_.size());

        if (size_ == 0) {
            return nil_id;
        }

        uint64_t key = make_key_(node_id, symb);
        assert(key != 0);

        // int cnt = 0;
        for (uint64_t i = Hasher::hash(key) & capa_size_.mask();; i = right_(i)) {
            // cnt += 1;
            if (i == 0) {
                // table_[0] is always empty so that table_[i] = 0 indicates to be empty.
                continue;
            }
            if (i == get_root()) {
                continue;
            }
            if (table_[i] == 0) {
                // encounter an empty slot
                return nil_id;
            }
            if (table_[i] == key) {
                // if(symb == 255) { // ダミーノードのみを対象としてするときに
                    // cnt_linear_proving[cnt] += 1;
                    // cnt_linear_proving_all += 1;
                // }
                return i;
            }
        }
    }

    bool add_child(uint64_t& node_id, uint64_t symb) {
        assert(node_id < capa_size_.size());
        assert(symb < symb_size_.size());

        uint64_t key = make_key_(node_id, symb);
        assert(key != 0);

        for (uint64_t i = Hasher::hash(key) & capa_size_.mask();; i = right_(i)) {
            if (i == 0) {
                // table_[0] is always empty so that any table_[i] = 0 indicates to be empty.
                continue;
            }

            if (i == get_root()) {
                continue;
            }

            if (table_[i] == 0) {
                // this slot is empty
                if (size_ == max_size_) {
                    return false;  // needs to expand
                }

                table_.set(i, key);

                ++size_;

                if(is_node_use_table(node_id)) set_not_use_hash_table(node_id, symb, i, 0);
                node_id = i;

                return true;
            }

            if (table_[i] == key) {
                node_id = i;
                return false;  // already stored
            }
        }
    }

    std::pair<uint64_t, uint64_t> get_parent_and_symb(uint64_t node_id) const {
        assert(node_id < capa_size_.size());

        uint64_t key = table_[node_id];
        if (key == 0) {
            // root or not exist
            return {nil_id, 0};
        }
        // Returns pair (parent, label)
        return std::make_pair(key >> symb_size_.bits(), key & symb_size_.mask());
    };

    class node_map {
      public:
        node_map() = default;

        node_map(compact_vector&& map, bit_vector&& done_flags)
            : map_{std::move(map)}, done_flags_{std::move(done_flags)} {}

        ~node_map() = default;

        uint64_t operator[](uint64_t i) const {
            return done_flags_[i] ? map_[i] : UINT64_MAX;
        }

        uint64_t size() const {
            return map_.size();
        }

        node_map(const node_map&) = delete;
        node_map& operator=(const node_map&) = delete;

        node_map(node_map&& rhs) noexcept = default;
        node_map& operator=(node_map&& rhs) noexcept = default;

      private:
        compact_vector map_;
        bit_vector done_flags_;
    };

    bool needs_to_expand() const {
        return max_size() <= size();
    }

    node_map expand(const std::vector<uint64_t>& cnt_number_of_children_per_node, uint64_t min_number_of_children, uint64_t lambda) {
        plain_bonsai_trie_table new_ht{capa_bits() + 1, symb_size_.bits()};
        new_ht.add_root();

#ifdef POPLAR_EXTRA_STATS
        new_ht.num_resize_ = num_resize_ + 1;
#endif

        bit_vector done_flags(capa_size());
        done_flags.set(get_root());

        table_.set(get_root(), new_ht.get_root());
        // if(cnt_number_of_children_per_node[get_root()] >= min_number_of_children) new_ht.set_is_use_table_elements(get_root());

        std::vector<std::pair<uint64_t, uint64_t>> path;
        path.reserve(256);

        // 0 is empty, 1 is root
        for (uint64_t i = 2; i < table_.size(); ++i) {
            if (done_flags[i] || table_[i] == 0) {
                // skip already processed or empty elements
                continue;
            }

            path.clear();
            uint64_t node_id = i;

            do {
                auto [parent, label] = get_parent_and_symb(node_id);
                assert(parent != nil_id);
                path.emplace_back(std::make_pair(node_id, label));
                node_id = parent;
            } while (!done_flags[node_id]);

            uint64_t new_node_id = table_[node_id];

            // before_parent : 元の配列上の親番号
            // after_parent  : new_ht上での親番号
            uint64_t before_parent = node_id;
            for (auto rit = std::rbegin(path); rit != std::rend(path); ++rit) {
                uint64_t after_parent = new_node_id;
                new_ht.add_child(new_node_id, rit->second);
                table_.set(rit->first, new_node_id);
                if(cnt_number_of_children_per_node[before_parent] >= min_number_of_children) {      // 対象とするノードの子ノードの数が条件を満たしているとき
                    new_ht.set_node_use_table(after_parent);                                 // 
                    new_ht.set_not_use_hash_table(after_parent, rit->second, new_node_id, lambda);  // 
                }
                done_flags.set(rit->first);
                before_parent = rit->first;
            }
        }

        node_map node_map{std::move(table_), std::move(done_flags)};
        std::swap(*this, new_ht);

        return node_map;
    }

    // not_use_hash_table内を2分探索するための関数
    std::pair<bool, uint64_t> binary_search_not_use_hash_table(uint64_t node_id) const {
        if(not_use_hash_table.size() == 0) return {false, 0};
        uint64_t left=0, right=not_use_hash_table.size()-1, mid;
        while(left <= right) {
            mid = (left + right) / 2;
            if(not_use_hash_table[mid].node_id == node_id) return {true, mid};
            else if(not_use_hash_table[mid].node_id < node_id) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        return {false, left};
    }

    // not_use_hash_tableに対し、値を更新する
    void set_not_use_hash_table(uint64_t node_id, uint64_t label, uint64_t next_node_id, uint64_t lambda) {
                                // 遷移前の番号,   遷移文字,        遷移後の番号,          matchの最大値
        // 手順
        // node_idでnot_use_hash_tableのnode_idを2分探索
        // 無ければ、insertし、サイズを拡張する(256*lambda_)
        // not_use_hash_table[pos].table[label] = next_node_idとし、完了
        auto [flag_hit, pos] = binary_search_not_use_hash_table(node_id);
        if(!flag_hit) { // 対象とする部分が存在しないとき
            not_use_hash_table.insert(not_use_hash_table.begin()+pos, table_info{node_id});
            not_use_hash_table[pos].table.resize(256*lambda, nil_id);
        }
        not_use_hash_table[pos].table[label] = next_node_id;
    }

    // is_use_table_elementsの値をtrueにする
    // node_id(ノード)がテーブルを必要とした際にtrueにする
    void set_node_use_table(uint64_t node_id) {
        node_use_table[node_id] = true;
    }

    // 使用しているのか(true : 使用, false : 未使用)
    bool is_node_use_table(uint64_t node_id) const {
        return node_use_table[node_id];
    }

    // テーブル内の特定の要素が使用されているのか
    bool is_use_table(uint64_t i) {
        if(table_[i] == 0) return false;
        return true;
    }

    // not_use_hash_tableを使用して、ノード遷移する際に使用する関数
    uint64_t transiton_not_use_hash_table(uint64_t node_id, uint64_t label) const {
        auto [flag_hit, pos] = binary_search_not_use_hash_table(node_id);
        if(!flag_hit) return nil_id; // nil_idを返値としないのは、使用していない部分の値が0だから
        return not_use_hash_table[pos].table[label];
    }

    void reset_cnt_compare() {
        std::cout << "--- reset_cnt_linear_proving ---" << std::endl;
        cnt_linear_proving.clear();
        cnt_linear_proving_all = 0;
    }

    void show_cnt_compare() {
        std::cout << "--- cnt_linear_proving ---" << std::endl;
        std::cout << "all : " << cnt_linear_proving_all << std::endl;
        uint64_t all = 0;
        for(auto p : cnt_linear_proving) {
            // std::cout << p.first << " : " << (long double)(p.second) / (long double)(cnt_linear_proving_all) << std::endl;
            // std::cout << p.first << " : " << p.second << std::endl;
            std::cout << p.first << std::endl;
            all += p.second;
        }
        std::cout << "---" << std::endl;
        for(auto p : cnt_linear_proving) {
            std::cout << p.second << std::endl;
        }
        std::cout << "all : " << all << std::endl;
    }

    // テーブルの値を確認する
    void check_not_use_hash() {
        std::cout << "--- check_not_use_hash ---" << std::endl;
        int cnt = 0;
        int size = int(node_use_table.size());
        for(int i=0; i < size; i++) {
            if(node_use_table[i]) cnt++;
        }
        std::cout << "cnt : " << cnt << std::endl;
        // for(int i=0; i < not_use_hash_table.size(); i++) {
        //     std::cout << "node_id : " << not_use_hash_table[i].node_id << std::endl;
        //     std::cout << "table_size : " << not_use_hash_table[i].table.size() << std::endl;
        //     cnt = 0;
        //     for(int j=0; j < int(not_use_hash_table[i].table.size()); j++) {
        //         if(not_use_hash_table[i].table[j] != 0) cnt++;
        //     }
        //     std::cout << "use_table_num : " << cnt << std::endl;
        // }
    }

    // # of registerd nodes
    uint64_t size() const {
        return size_;
    }
    uint64_t max_size() const {
        return max_size_;
    }
    uint64_t capa_size() const {
        return capa_size_.size();
    }
    uint32_t capa_bits() const {
        return capa_size_.bits();
    }
    uint64_t symb_size() const {
        return symb_size_.size();
    }
    uint32_t symb_bits() const {
        return symb_size_.bits();
    }
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_resize() const {
        return num_resize_;
    }
#endif
    uint64_t alloc_bytes() const {
        return table_.alloc_bytes();
    }

    void show_stats(std::ostream& os, int n = 0) const {
        auto indent = get_indent(n);
        show_stat(os, indent, "name", "plain_bonsai_trie_table");
        show_stat(os, indent, "factor", double(size()) / capa_size() * 100);
        show_stat(os, indent, "max_factor", MaxFactor);
        show_stat(os, indent, "size", size());
        show_stat(os, indent, "alloc_bytes", alloc_bytes());
        show_stat(os, indent, "capa_bits", capa_bits());
        show_stat(os, indent, "symb_bits", symb_bits());
#ifdef POPLAR_EXTRA_STATS
        show_stat(os, indent, "num_resize", num_resize_);
#endif
    }

    plain_bonsai_trie_table(const plain_bonsai_trie_table&) = delete;
    plain_bonsai_trie_table& operator=(const plain_bonsai_trie_table&) = delete;

    plain_bonsai_trie_table(plain_bonsai_trie_table&&) noexcept = default;
    plain_bonsai_trie_table& operator=(plain_bonsai_trie_table&&) noexcept = default;

  private:
    compact_vector table_;
    uint64_t size_ = 0;  // # of registered nodes
    uint64_t max_size_ = 0;  // MaxFactor% of the capacity
    size_p2 capa_size_;
    size_p2 symb_size_;
    std::vector<bool> node_use_table;    // テーブルに吐き出しているのか
    std::vector<table_info> not_use_hash_table; // ハッシュテーブルを使用せずにノード遷移するための配列
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_resize_ = 0;
#endif

    uint64_t make_key_(uint64_t node_id, uint64_t symb) const {
        return (node_id << symb_size_.bits()) | symb;
    }
    uint64_t right_(uint64_t slot_id) const {
        return (slot_id + 1) & capa_size_.mask();
    }
};

}  // namespace poplar

#endif  // POPLAR_TRIE_PLAIN_BONSAI_TRIE_TABLE_HPP
