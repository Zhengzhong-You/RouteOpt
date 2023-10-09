//
// Created by Zhengzhong You on 3/7/23.
//

#include "CVRP.hpp"
#include "experimental/filesystem"

using namespace std;

#ifdef if_draw_BBT_graph

void checkIfFolderExist(const string &folderName) {
  if (!std::filesystem::is_directory(folderName)) {
    std::filesystem::create_directory(folderName);
    std::cout << "Folder created: " << folderName << std::endl;
  }
}

void CVRP::drawBBTgraph() {
  ofstream file;
  string dir_dot = "BBTgraphDot/";
  string dir_pdf = "BBTgraphPdf/";
  auto tmp_dot_name = dir_dot + "BBTgraph_" + FileName + ".dot";
  auto tmp_pdf_name = dir_pdf + "BBTgraph_" + FileName + ".pdf";

  checkIfFolderExist(dir_dot);
  checkIfFolderExist(dir_pdf);
  // draw
  file.open(tmp_dot_name);
  file << "digraph G {" << endl;
  stringstream ss;

#ifdef AccuracyTest
  ss << "\"" << "Idx | LB | UB | M1 | M2 | Brc" << "\"" << endl;
#else
  ss << "\"" << "Idx | LB | UB | Brc" << "\"" << endl;
#endif

  for (auto &node : BBT_nodes_name) {
    if (node.first == 0) {
      ss << "\"" << node.second.second << "\"" << endl;
    }
    if (node.second.first.first != -1 && !BBT_nodes_name[node.second.first.first].second.empty()) {
      ss << "\"" << node.second.second << "\"" << " -> " << "\"" << BBT_nodes_name[node.second.first.first].second
         << "\"" << endl;
    }
    if (node.second.first.second != -1 && !BBT_nodes_name[node.second.first.second].second.empty()) {
      ss << "\"" << node.second.second << "\"" << " -> " << "\"" << BBT_nodes_name[node.second.first.second].second
         << "\"" << endl;
    }
  }

  file << ss.str();
  file << "}" << endl;

  file.close();

  system(("dot -Tpdf " + tmp_dot_name + " -o " + tmp_pdf_name).c_str());
}

void CVRP::constructBBNodeName(BBNODE *node, const std::string &terminateReason) {
  string name;
  string pattern = " | ";
  name += to_string(NumExploredNodes - 1) + pattern;
  name += doubleToString(ceil_transformed_number_related(node->Val)) + pattern;
  name += doubleToString(UB) + pattern;
  if (terminateReason.empty()) {
#ifdef AccuracyTest
    name += to_string(rank_phase1) + pattern;
    name += to_string(rank_phase2) + pattern;
#endif
    auto &info = node->BrCs.back().Edge;
    name += "(" + to_string(info.first) + "," + to_string(info.second) + ")";
    BBT_nodes_name[node->Idx] = {{IdxNode + 1, IdxNode + 2}, name};
  } else {
#ifdef AccuracyTest
    name += "-1" + pattern;
    name += "-1" + pattern;
#endif
    name += terminateReason;
    BBT_nodes_name[node->Idx] = {{-1, -1}, name};
  }
}
#endif

std::string doubleToString(double d) {
  std::string s = std::to_string(d);
  s.erase(s.find_last_not_of('0') + 1, std::string::npos);
  if (s.back() == '.') {
    s.pop_back();
  }
  return s;
}
