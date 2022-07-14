/**
 * MIT License
 *
 * Copyright (c) 2018–2019 Shunsuke Kanda
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef POPLAR_TRIE_POPLAR_HPP
#define POPLAR_TRIE_POPLAR_HPP

#include "poplar/compact_bonsai_trie.hpp"
#include "poplar/compact_fkhash_trie.hpp"
#include "poplar/plain_bonsai_trie.hpp"
#include "poplar/plain_fkhash_trie.hpp"

#include "poplar/compact_bonsai_nlm.hpp"
#include "poplar/compact_fkhash_nlm.hpp"
#include "poplar/plain_bonsai_nlm.hpp"
#include "poplar/plain_fkhash_nlm.hpp"

#include "poplar/map.hpp"

#include "check/map_check.hpp"
#include "check/plain_bonsai_nlm_check.hpp"
#include "check/plain_bonsai_tire_check.hpp"

#include "dynamic_replacement/map_dr.hpp"
#include "dynamic_replacement/plain_bonsa_trie_dr.hpp"
#include "dynamic_replacement/plain_bonsai_nlm_dr.hpp"

#include "all_key_restore/map_akr.hpp"
#include "all_key_restore/plain_bonsai_trie_akr.hpp"
#include "all_key_restore/plain_bonsai_nlm_akr.hpp"

#include "table/map_table.hpp"
#include "table/plain_bonsai_trie_table.hpp"
#include "table/plain_bonsai_nlm_table.hpp"

#include "skip_dummy/map_SD.hpp"
#include "skip_dummy/plain_bonsai_trie_SD.hpp"
#include "skip_dummy/plain_bonsai_nlm_SD.hpp"

namespace poplar {

template <typename Value>
using plain_bonsai_map = map<plain_bonsai_trie<>, plain_bonsai_nlm<Value>>;

template <typename Value, uint64_t ChunkSize = 16>
using semi_compact_bonsai_map = map<plain_bonsai_trie<>, compact_bonsai_nlm<Value, ChunkSize>>;

template <typename Value, uint64_t ChunkSize = 16>
using compact_bonsai_map = map<compact_bonsai_trie<>, compact_bonsai_nlm<Value, ChunkSize>>;

template <typename Value>
using plain_fkhash_map = map<plain_fkhash_trie<>, plain_fkhash_nlm<Value>>;

template <typename Value, uint64_t ChunkSize = 16>
using semi_compact_fkhash_map = map<plain_fkhash_trie<>, compact_fkhash_nlm<Value, ChunkSize>>;

template <typename Value, uint64_t ChunkSize = 16>
using compact_fkhash_map = map<compact_fkhash_trie<>, compact_fkhash_nlm<Value, ChunkSize>>;

// 何かの実装をする際に使用する
template <typename Value>
using plain_bonsai_map_check = map_check<plain_bonsai_trie_check<>, plain_bonsai_nlm_check<Value>>;

// 動的に組み替える際に使用する
template <typename Value>
using plain_bonsai_map_dr = map_dr<plain_bonsai_trie_dr<>, plain_bonsai_nlm_dr<Value>>;

// 全てのキーを復元して、新しい辞書に登録する際に使用
template <typename Value>
using plain_bonsai_map_akr = map_akr<plain_bonsai_trie_akr<>, plain_bonsai_nlm_akr<Value>>;

// テーブルを使用
template <typename Value>
using plain_bonsai_map_table = map_table<plain_bonsai_trie_table<>, plain_bonsai_nlm_table<Value>>;

// 検索する際のダミーノードをスキップ
template <typename Value>
using plain_bonsai_map_SD = map_SD<plain_bonsai_trie_SD<>, plain_bonsai_nlm_SD<Value>>;

}  // namespace poplar

#endif  // POPLAR_TRIE_POPLAR_HPP
