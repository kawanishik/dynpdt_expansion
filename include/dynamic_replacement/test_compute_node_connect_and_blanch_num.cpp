#include "map_dr.hpp"
#include "plain_bonsa_trie_dr.hpp"
#include "plain_bonsai_nlm_dr.hpp"
#include <cstdio>
#include <vector>
#include <string>

using namespace poplar;

int main() {
  std::vector<std::string> keys = {
      "aaaaaaaaa",
      "aaaaaabbb",
      "aaaaaaccc",
      "aaabbbbbb",
      "aaabbbccc",
      "bbbbbbbbb",
      "bbbcccccc",
      "bbbccccccddd"
  };
  map_dr<plain_bonsai_trie_dr<>, plain_bonsai_nlm_dr<int>> mp(8);
  for (auto it = keys.rbegin(); it != keys.rend(); ++it) {
    mp.update(*it);
  }
  std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children;
  // std::vector<uint64_t> blanch_num_except_zero;
  std::vector<uint64_t> cnt_leaf_per_node;
  // mp.compute_node_connect_and_blanch_num(children, blanch_num_except_zero, cnt_leaf_per_node);
  mp.compute_node_connect_and_blanch_num(children, cnt_leaf_per_node);

  auto show_internal = [&]() {
    std::cerr << "children = {" << std::endl;
    for (size_t i = 0; i < children.size(); ++i) {
      auto& cl = children[i];
      if (cl.empty()) continue;
      std::cerr << "\t["<<i<<"] ""{" << std::endl << "\t\t";
      for (auto c : cl) {
        std::cerr << "{ " << c.first << ", " << c.second << " }," ;
      }
      std::cerr << std::endl << "\t}," << std::endl;
    }
    std::cerr << "}" << std::endl;
    // std::cerr << "blanch_num_except_zero = {" << std::endl;
    // for (size_t i = 0; i < blanch_num_except_zero.size(); ++i) {
    //   auto v = blanch_num_except_zero[i];
    //   if (v == 0) continue;
    //   std::cerr << "\t[" << i << "] " << v << ", " << std::endl;
    // }
    std::cerr << "}" << std::endl;
    std::cerr << "cnt_leaf_per_node = {" << std::endl;
    for (size_t i = 0; i < cnt_leaf_per_node.size(); ++i) {
      auto v = cnt_leaf_per_node[i];
      if (v == 0) continue;
      std::cerr << "\t[" << i << "] " << v << ", " << std::endl;
    }
    std::cerr << "}" << std::endl;
  };

  std::vector<std::pair<size_t, std::vector<std::pair<uint64_t,uint64_t>>>> children_exp {
      {1, {{ 9, 22042 },{ 0, 38070 },{ 3, 38084 }}},
        {1920, {
                { 2, 16923 },{ 2, 42257 },
        }},
        {38070, {
                { 2, 1920 },{ 5, 51392 },
        }},
  };
  // std::vector<std::pair<size_t,uint64_t>> blanch_num_except_zero_exp {
  //   {0, 8},
  //   {1, 2},
  //   {1920, 2},
  //   {38070, 4},
  // };
  std::vector<std::pair<size_t,uint64_t>> cnt_leaf_per_node_exp {
    {0, 8},
    {1, 8},
    {1920, 3},
    {16923, 1},
    {22042, 1},
    {38070, 5},
    {38084, 1},
    {42257, 1},
    {51392, 1},
  };

  for (auto [i,cle] : children_exp) {
    if (cle != children[i]) {
      show_internal();
      return EXIT_FAILURE;
    }
  }
  // for (auto [i,v] : blanch_num_except_zero_exp) {
  //   if (v != blanch_num_except_zero[i]) {
  //     show_internal();
  //     return EXIT_FAILURE;
  //   }
  // }
  for (auto [i,v] : cnt_leaf_per_node_exp) {
    if (v != cnt_leaf_per_node[i]) {
      show_internal();
      return EXIT_FAILURE;
    }
  }
  std::cout << "OK" << std::endl;
}