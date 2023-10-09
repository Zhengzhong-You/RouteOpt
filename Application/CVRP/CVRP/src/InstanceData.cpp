//
// Created by Zhengzhong You on 12/10/22.
//

#include "InstanceData.hpp"
#include "MACRO.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;
void takeDataFromFile(std::string &file_name, InstanceData &data) {
  std::ifstream file(file_name);
  if (!file.is_open()) {
    std::cout << "File not found!" << std::endl;
    safe_Hyperparameter(EXIT_FAILURE)
  }

  bool if_pop = false;
  string tmp_name = {};
  for (auto i = file_name.end() - 2; i >= file_name.begin(); --i) {
    if (*(i + 1) == '.') {
      if_pop = true;
    }
    if (*i == '/' || *i == '\\')break;
    if (if_pop) {
      data.name.push_back(*i);
    }
    tmp_name.push_back(*i);
  }
  if (!if_pop) data.name = tmp_name;

  std::reverse(data.name.begin(), data.name.end());

  std::string tmp_string, tmp_string_2;
#ifdef SOLVER_VRPTW
  while (true) {
    if (tmp_string.find("NUMBER") != std::string::npos) {
      std::getline(file, tmp_string);
      std::stringstream temp_ss(tmp_string);
      temp_ss >> data.K;
      temp_ss >> data.Cap;
      break;
    } else std::getline(file, tmp_string);
  }
  while (true) {
    std::getline(file, tmp_string);
    if (tmp_string.find("CUST NO.") != std::string::npos) {
      std::getline(file, tmp_string);
      break;
    }
  }
  int num = MaxNum_Customers, tmp_j = 0;
  data.InfoVertex.resize(num, std::vector<double>(7));
  while (std::getline(file, tmp_string) && tmp_j < num) {
    if (tmp_string.empty()) break;
    std::stringstream temp_ss(tmp_string);
    temp_ss >> data.InfoVertex[tmp_j][0];
    temp_ss >> data.InfoVertex[tmp_j][1];
    temp_ss >> data.InfoVertex[tmp_j][2];
    temp_ss >> data.InfoVertex[tmp_j][3];
    temp_ss >> data.InfoVertex[tmp_j][4];
    temp_ss >> data.InfoVertex[tmp_j][5];
    temp_ss >> data.InfoVertex[tmp_j][6];
    ++tmp_j;
  }
  data.Dim = tmp_j;
  data.InfoVertex.resize(data.Dim);
#else
  //read n and k
//  std::getline(file, tmp_string);
//  auto index_k = tmp_string.find('k'), index_n = tmp_string.find('n');
//  data.Dim = std::stoi(tmp_string.substr(index_n + 1, index_k - index_n - 2));
//  data.K = std::stoi(tmp_string.substr(index_k + 1, tmp_string.size() - index_k - 1));
//  cout << "Dim: " << data.Dim << " K: " << data.K << endl;
  string line;
  string name;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string key;
    iss >> key;
    if (key == "NAME") {
      iss >> key; // Discard the ":"
      iss >> name;
      if (name.find('k') != std::string::npos) {
        size_t k_pos = name.find_last_of('k');
        data.K = std::stoi(name.substr(k_pos + 1));
      } else {
        data.K = -1;
      }
    } else if (key == "DIMENSION") {
      iss >> key; // Discard the ":"
      iss >> data.Dim;
      break;
    }
  }

  if (data.K == -1) {
    data.K = data.Dim - 1;
  }

  std::getline(file, tmp_string);
  while (true) {
    if (tmp_string.find("CAPACITY") != std::string::npos) {
      std::stringstream temp_ss(tmp_string);
      temp_ss >> tmp_string_2;
      temp_ss >> tmp_string_2;
      temp_ss >> data.Cap;
      break;
    } else std::getline(file, tmp_string);
  }
  int tmp_j = 0;
  data.InfoVertex.resize(data.Dim, vector<double>(4));
  std::getline(file, tmp_string);
  while (true) {
    std::getline(file, tmp_string);
    if (tmp_string.find("DEMAND_SECTION") != string::npos || tmp_j == data.Dim) {
      break;
    } else {
      stringstream temp_ss(tmp_string);
      temp_ss >> data.InfoVertex[tmp_j][0];
      temp_ss >> data.InfoVertex[tmp_j][1];
      temp_ss >> data.InfoVertex[tmp_j][2];
      ++tmp_j;
    }
  }
  while (true) {
    if (tmp_string.find("DEMAND_SECTION") != string::npos) {
      break;
    }
    std::getline(file, tmp_string);
  }
  int tmp_int;
  tmp_j = 0;
  while (true) {
    std::getline(file, tmp_string);
    if (tmp_string.find("DEPOT_SECTION") != string::npos || tmp_j == data.Dim) {
      break;
    } else {
      stringstream temp_ss(tmp_string);
      temp_ss >> tmp_int;
      temp_ss >> data.InfoVertex[tmp_j][3];
      ++tmp_j;
    }
  }
//  data.Dim = int(tmp_j / 3);
//  data.Cap *= 1.5;
//  data.Cap = int(data.Cap);
//  cout << "Capacity changed!!!!!!! " << data.Cap << endl;
  data.InfoVertex.resize(data.Dim);
#endif
  file.close();
}