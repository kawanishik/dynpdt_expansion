#ifndef POPLAR_TRIE_PLAIN_BONSAI_TRIE_SD_HPP
#define POPLAR_TRIE_PLAIN_BONSAI_TRIE_SD_HPP

#include "../poplar/bit_tools.hpp"
#include "../poplar/bit_vector.hpp"
#include "../poplar/compact_vector.hpp"
#include "../poplar/hash.hpp"

namespace poplar {

template <uint32_t MaxFactor = 90, typename Hasher = hash::vigna_hasher>
class plain_bonsai_trie_SD {
  private:
    static_assert(0 < MaxFactor and MaxFactor < 100);

    struct skip_dummy_info {
        uint64_t node_id;
        std::vector<uint64_t> skip;
    };

  public:
    static constexpr uint64_t nil_id = UINT64_MAX;
    static constexpr uint32_t min_capa_bits = 16;

    static constexpr auto trie_type_id = trie_type_ids::BONSAI_TRIE;

  public:
    plain_bonsai_trie_SD() = default;

    plain_bonsai_trie_SD(uint32_t capa_bits, uint32_t symb_bits) {
        capa_size_ = size_p2{std::max(min_capa_bits, capa_bits)};
        symb_size_ = size_p2{symb_bits};
        max_size_ = static_cast<uint64_t>(capa_size_.size() * MaxFactor / 100.0);
        table_ = compact_vector{capa_size_.size(), capa_size_.bits() + symb_size_.bits()};
    }

    ~plain_bonsai_trie_SD() = default;

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
                node_id = i;

                return true;
            }

            if (table_[i] == key) {
                node_id = i;
                return false;  // already stored
            }
        }
    }

    // skip_dummy内を2分探索するための関数
    std::pair<bool, uint64_t> binary_search_skip_dummy(uint64_t node_id) const {
        int left = 0;
        int right = skip_dummy.size()-1;
        int mid;
        while(left <= right) {
            mid = (left + right) / 2;
            if(skip_dummy[mid].node_id == node_id) {
                return {true, mid};
            } else if(skip_dummy[mid].node_id < node_id) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        return {false, left};
    }

    // スキップするダミーノードに対し、追加する必要がある場合は、配列を拡張し、追加する
    // 必要じゃない場合は更新して終了する
    void add_new_skip_dummy(uint64_t pre_node_id, uint64_t node_id, uint64_t index) {
        if(skip_dummy.size() == 0) {    // 初期の動作
            skip_dummy.resize(1);
            skip_dummy[0].node_id = pre_node_id;
            skip_dummy[0].skip.resize(index+1);
            skip_dummy[0].skip[index] = node_id;
            return;
        }
        // pre_node_idで2分探索
        auto [flag_hit, pos] = binary_search_skip_dummy(pre_node_id);
        if(flag_hit) { // 該当するnode_idが見つかったので、indexが存在するのかの確認
            if(index >= skip_dummy[pos].skip.size()) skip_dummy[pos].skip.resize(index+1);
            skip_dummy[pos].skip[index] = node_id;
        } else { // 該当するnode_idが存在しないということなので、insertが必要
            skip_dummy.insert(skip_dummy.begin()+pos, skip_dummy_info{pre_node_id});
            skip_dummy[pos].skip.resize(index+1);
            skip_dummy[pos].skip[index] = node_id;
        }
    }

    // findする際に、スキップするdummyを見つける
    uint64_t find_skip_dummy(uint64_t node_id, uint64_t index) const {
        auto [flag_hit, pos] = binary_search_skip_dummy(node_id);
        if(!flag_hit) return nil_id;    // 対象とするノードが見つからないとき

        if(index >= skip_dummy[pos].skip.size()) return nil_id;
        return skip_dummy[pos].skip[index];
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
    }

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

    node_map expand() {
        plain_bonsai_trie_SD new_ht{capa_bits() + 1, symb_size_.bits()};
        new_ht.add_root();

#ifdef POPLAR_EXTRA_STATS
        new_ht.num_resize_ = num_resize_ + 1;
#endif

        bit_vector done_flags(capa_size());
        done_flags.set(get_root());

        table_.set(get_root(), new_ht.get_root());

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

            for (auto rit = std::rbegin(path); rit != std::rend(path); ++rit) {
                new_ht.add_child(new_node_id, rit->second);
                table_.set(rit->first, new_node_id);
                done_flags.set(rit->first);
            }
        }

        // ここでnew_htのskip_dummyに元の情報からtable_を使用して更新する
        new_ht.skip_dummy_resize(skip_dummy.size());
        for(int i=0; i < int(skip_dummy.size()); i++) {
            new_ht.skip_dummy_set_node_id(i, table_[skip_dummy[i].node_id]);
            new_ht.skip_dummy_skip_resize(i, skip_dummy[i].skip.size());
            std::vector<uint64_t> skip;
            for(int j=0; j < int(skip_dummy[i].skip.size()); j++) {
                skip.emplace_back(table_[skip_dummy[i].skip[j]]);
            }
            new_ht.skip_dummy_set_skip(i, skip);
        }

        new_ht.skip_dummy_sort();

        node_map node_map{std::move(table_), std::move(done_flags)};
        std::swap(*this, new_ht);

        return node_map;
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

    // skip_dummy...から始まる関数はskip_dummyに対する処理を記述している
    void skip_dummy_resize(uint64_t size) {
        skip_dummy.resize(size);
    }

    void skip_dummy_skip_resize(uint64_t pos, uint64_t size) {
        skip_dummy[pos].skip.resize(size);
    }

    void skip_dummy_set_node_id(uint64_t pos, uint64_t node_id) {
        skip_dummy[pos].node_id = node_id;
    }

    void skip_dummy_set_skip(uint64_t pos, const std::vector<uint64_t>& skip) {
        for(int i=0; i < skip.size(); i++) skip_dummy[pos].skip[i] = skip[i];
    }

    void skip_dummy_sort() {
        std::sort(skip_dummy.begin(), skip_dummy.end(), [] (auto l, auto r) {
            return l.node_id < r.node_id;
        });
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
        show_stat(os, indent, "name", "plain_bonsai_trie_SD");
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

    plain_bonsai_trie_SD(const plain_bonsai_trie_SD&) = delete;
    plain_bonsai_trie_SD& operator=(const plain_bonsai_trie_SD&) = delete;

    plain_bonsai_trie_SD(plain_bonsai_trie_SD&&) noexcept = default;
    plain_bonsai_trie_SD& operator=(plain_bonsai_trie_SD&&) noexcept = default;

  private:
    compact_vector table_;
    uint64_t size_ = 0;  // # of registered nodes
    uint64_t max_size_ = 0;  // MaxFactor% of the capacity
    size_p2 capa_size_;
    size_p2 symb_size_;
    std::vector<skip_dummy_info> skip_dummy;
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

#endif  // POPLAR_TRIE_PLAIN_BONSAI_TRIE_SD_HPP
