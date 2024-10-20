
#ifndef INSTANCEDATA_HPP_
#define INSTANCEDATA_HPP_

#include <string>
#include <vector>

struct InstanceData {
  std::string name{};
  int k{};
  double cap{};
  int dim{};
  std::vector<std::vector<double>> info_vertex;
};

void takeDataFromFile(std::string &file_name, InstanceData &data);

#endif //INSTANCEDATA_HPP_
