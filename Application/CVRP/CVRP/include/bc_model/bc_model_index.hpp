//
// Created by You, Zhengzhong on 2/18/24.
//

#ifndef INCLUDE_NEWMODELINDEX_HPP_
#define INCLUDE_NEWMODELINDEX_HPP_

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include "macro.hpp"

struct TripleHasher {
  size_t operator()(const std::tuple<int, int, int> &V) const {
	return std::get<0>(V) * MAX_NUM_CUSTOMERS * MAX_NUM_CUSTOMERS + std::get<1>(V) * MAX_NUM_CUSTOMERS + std::get<2>(V);
  }
};

class VariableIndexMap {
 public:
  std::unordered_map<char, std::unordered_map<int, int>> single_varMap;
  std::unordered_map<int, std::pair<char, int>> single_index2ArcMap;

  std::unordered_map<char, std::unordered_map<std::pair<int, int>, int, PairHasher>> varMap;
  std::unordered_map<int, std::tuple<char, int, int>> index2ArcMap;

  std::unordered_map<char, std::unordered_map<std::tuple<int, int, int>, int, TripleHasher>>
	  varMap3;
  std::unordered_map<int, std::tuple<char, int, int, int>> index2ArcMap3;

  void insert(int key, char varType, int value) {
	single_varMap[varType].emplace(key, value);
	single_index2ArcMap.emplace(value, std::make_pair(varType, key));
  }

  int getIndex(int key, char varType) {
	return single_varMap[varType].at(key);
  }

  void insert(const std::pair<int, int> &keyPair, char varType, int value) {
	varMap[varType].emplace(keyPair, value);
	index2ArcMap.emplace(value, std::make_tuple(varType, keyPair.first, keyPair.second));
  }

  int getIndex(const std::pair<int, int> &keyPair, char varType) {
	return varMap[varType].at(keyPair);
  }

  void insert3(const std::tuple<int, int, int> &keyPair, char varType, int value) {
	varMap3[varType].emplace(keyPair, value);
	index2ArcMap3.emplace(value,
						  std::make_tuple(varType, std::get<0>(keyPair), std::get<1>(keyPair), std::get<2>(keyPair)));
  }

  int getIndex(const std::tuple<int, int, int> &keyPair, char varType) {
	return varMap3[varType].at(keyPair);
  }

  std::tuple<char, int, int, int> getVar(int index) {
	if (index2ArcMap.find(index) != index2ArcMap.end()) {
	  auto &tmp = index2ArcMap.at(index);
	  return std::make_tuple(std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), -1);
	} else if (single_index2ArcMap.find(index) != single_index2ArcMap.end()) {
	  auto &tmp = single_index2ArcMap.at(index);
	  return std::make_tuple(tmp.first, tmp.second, -1, -1);
	} else return index2ArcMap3.at(index);
  }

  void clear() {
	single_varMap.clear();
	single_index2ArcMap.clear();
	varMap.clear();
	index2ArcMap.clear();
	varMap3.clear();
	index2ArcMap3.clear();
  }

  bool testIfExist(int key, char varType) {
	return single_varMap[varType].find(key) != single_varMap[varType].end();
  }

  bool testIfExist(const std::tuple<int, int, int> &keyPair, char varType) {
	return varMap3[varType].find(keyPair) != varMap3[varType].end();
  }

  bool testIfExist(const std::pair<int, int> &keyPair, char varType) {
	return varMap[varType].find(keyPair) != varMap[varType].end();
  }
};

#endif //INCLUDE_NEWMODELINDEX_HPP_
