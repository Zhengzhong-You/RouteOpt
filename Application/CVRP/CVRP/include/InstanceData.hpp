//
// Created by Zhengzhong You on 12/10/22.
//

#ifndef INSTANCEDATA_HPP_
#define INSTANCEDATA_HPP_

#include <string>
#include <vector>

struct InstanceData {
  std::string name{};
  int K{};
  double Cap{};
  int Dim{};
  std::vector<std::vector<double>> InfoVertex;
};

void takeDataFromFile(std::string &file_name, InstanceData &data);

#endif //INSTANCEDATA_HPP_
