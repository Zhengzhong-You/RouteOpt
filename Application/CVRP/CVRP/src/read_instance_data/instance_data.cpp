
#include "instance_data.hpp"
#include "macro.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "config.hpp"
using namespace std;

void extractName(const std::string &file_name, InstanceData &data);

void parseVRPTWOrRVRPSTW(std::ifstream &file, InstanceData &data);

void parseVRPTWSecondType(std::ifstream &file, InstanceData &data);

void parseCVRPSolver(std::ifstream &file, InstanceData &data);

void takeDataFromFile(std::string &file_name, InstanceData &data) {
  std::ifstream file(file_name);
  if (!file.is_open()) {
	std::cout << "File not found!" << std::endl;
	safe_Hyperparameter(EXIT_FAILURE);
  }

  extractName(file_name, data);

#if defined(SOLVER_VRPTW) || defined(SOLVER_RVRPSTW)
  parseVRPTWSecondType(file, data);
  if (data.dim == 0) parseVRPTWOrRVRPSTW(file, data);
#else
  parseCVRPSolver(file, data);
#endif

  file.close();
}

void extractName(const std::string &file_name, InstanceData &data) {
  bool if_pop = false;
  std::string tmp_name;
  for (auto i = file_name.end() - 2; i >= file_name.begin(); --i) {
	if (*(i + 1) == '.') {
	  if_pop = true;
	}
	if (*i == '/' || *i == '\\') break;
	if (if_pop) {
	  data.name.push_back(*i);
	}
	tmp_name.push_back(*i);
  }
  if (!if_pop) data.name = tmp_name;
  std::reverse(data.name.begin(), data.name.end());
}

void parseVRPTWOrRVRPSTW(std::ifstream &file, InstanceData &data) {
  std::string tmp_string;
  while (true) {
	if (tmp_string.find("NUMBER") != std::string::npos) {
	  std::getline(file, tmp_string);
	  std::stringstream temp_ss(tmp_string);
	  temp_ss >> data.k;
	  temp_ss >> data.cap;
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

  int num = MAX_NUM_CUSTOMERS, tmp_j = 0;
  data.info_vertex.resize(num, std::vector<double>(7));
  while (std::getline(file, tmp_string) && tmp_j < num) {
	if (tmp_string.empty()) break;
	std::stringstream temp_ss(tmp_string);
	temp_ss >> data.info_vertex[tmp_j][0]; // index
	temp_ss >> data.info_vertex[tmp_j][1]; // x
	temp_ss >> data.info_vertex[tmp_j][2]; // y
	temp_ss >> data.info_vertex[tmp_j][3]; // demand
	temp_ss >> data.info_vertex[tmp_j][4]; // ready time
	temp_ss >> data.info_vertex[tmp_j][5]; // due date
	temp_ss >> data.info_vertex[tmp_j][6]; // service time
	++tmp_j;
  }
  data.dim = tmp_j;
  data.info_vertex.resize(data.dim);
}

bool checkForTravelTimeMatrix(std::ifstream &file) {
  std::string line;
  while (std::getline(file, line)) {
	if (line.find("Travel Time Matrix") != std::string::npos) {
	  return true;
	}
  }
  return false;
}

void parseVRPTWSecondType(std::ifstream &file, InstanceData &data) {
  file.clear();
  file.seekg(0, std::ios::beg);

  if (!checkForTravelTimeMatrix(file)) {
	file.clear();
	file.seekg(0, std::ios::beg);
	return;
  }

  file.clear();
  file.seekg(0, std::ios::beg);

  std::string tmp_string;

  std::string line;
  while (std::getline(file, line)) {
	std::stringstream ss(line);
	if (line.find("Dimension:") != std::string::npos) {
	  ss.ignore(11); // Ignore the label part
	  ss >> data.dim;
	} else if (line.find("Capacity:") != std::string::npos) {
	  ss.ignore(10); // Ignore the label part
	  ss >> data.cap;
	} else if (line.find("UB:") != std::string::npos) {
	  ss.ignore(4); // Ignore the label part
	  ss >> Config::ub;
	} else if (line.find("Vertex Information:") != std::string::npos) {
	  break; // Stop after reading all required headers
	}
  }

  std::getline(file, line); // Skip the header line with ID, X, Y, etc.

  data.info_vertex.resize(data.dim, std::vector<double>(7));
  int index = 0;
  while (std::getline(file, line) && index < data.dim) {
	if (line.empty()) continue; // Skip empty lines
	std::stringstream ss(line);
	for (int j = 0; j < 7; ++j) {
	  ss >> data.info_vertex[index][j];
	}
	++index;
  }
}

void parseCVRPSolver(std::ifstream &file, InstanceData &data) {
  std::string line, name, tmp_string, tmp_string_2;
  while (std::getline(file, line)) {
	std::istringstream iss(line);
	std::string key;
	iss >> key;
	if (key == "NAME") {
	  iss >> key; // Discard the ":"
	  iss >> name;
	  if (name.find('k') != std::string::npos) {
		size_t k_pos = name.find_last_of('k');
		data.k = std::stoi(name.substr(k_pos + 1));
	  } else {
		data.k = -1;
	  }
	} else if (key == "DIMENSION") {
	  iss >> key; // Discard the ":"
	  iss >> data.dim;
	  break;
	}
  }

  if (data.k == -1) {
	data.k = data.dim - 1;
  }

  std::getline(file, tmp_string);
  while (true) {
	if (tmp_string.find("CAPACITY") != std::string::npos) {
	  std::stringstream temp_ss(tmp_string);
	  temp_ss >> tmp_string_2;
	  temp_ss >> tmp_string_2;
	  temp_ss >> data.cap;
	  break;
	} else std::getline(file, tmp_string);
  }

  int tmp_j = 0;
  data.info_vertex.resize(data.dim, std::vector<double>(4));
  std::getline(file, tmp_string);
  while (true) {
	std::getline(file, tmp_string);
	if (tmp_string.find("DEMAND_SECTION") != std::string::npos || tmp_j == data.dim) {
	  break;
	} else {
	  std::stringstream temp_ss(tmp_string);
	  temp_ss >> data.info_vertex[tmp_j][0];
	  temp_ss >> data.info_vertex[tmp_j][1];
	  temp_ss >> data.info_vertex[tmp_j][2];
	  ++tmp_j;
	}
  }

  while (true) {
	if (tmp_string.find("DEMAND_SECTION") != std::string::npos) {
	  break;
	}
	std::getline(file, tmp_string);
  }

  int tmp_int;
  tmp_j = 0;
  while (true) {
	std::getline(file, tmp_string);
	if (tmp_string.find("DEPOT_SECTION") != std::string::npos || tmp_j == data.dim) {
	  break;
	} else {
	  std::stringstream temp_ss(tmp_string);
	  temp_ss >> tmp_int;
	  temp_ss >> data.info_vertex[tmp_j][3];
	  ++tmp_j;
	}
  }
  data.info_vertex.resize(data.dim);
}
