//
// Created by Zhengzhong You on 5/31/22.
//

#ifndef CVRP_BBNODE_HPP
#define CVRP_BBNODE_HPP

#include <vector>
#include <unordered_map>
#include "MACRO.hpp"
#include "CUTS.hpp"
#include <deque>
#include <Eigen/Sparse>
#include "solver.hpp"

class CVRP;

class BBNODE;

struct Bucket;

using sparseRowMatrixXd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using RowVectorXT = Eigen::Matrix<size_t, 1, Eigen::Dynamic>;
using RowVectorXd = Eigen::Matrix<double, 1, Eigen::Dynamic>;
using VectorXi = Eigen::VectorXi;

struct BidirLinkList {
  //这个是选择完br以后，更新自己孩子的父辈，和自己的孩子
  BidirLinkList *PNode{nullptr}, *LNode{nullptr}, *RNode{nullptr};

  //偶数位置是column的lp位置，偶数位置是column的访问边的次数
  //此处还是不用int*比较好，大量的时间在此处的优化意义不大
  //而且访存的速度几乎一样，而我实现转移数据的操作并没有很多
  //vector 和 point的区别还是在于维度上面，维度多了以后就少用vector了
  //这里允许vector接受多个重复的值，如果一个col多次访问某边的话，但是这些重复的值是连续的（取决于我的实现方式）
  std::unordered_map<int, std::vector<int>> Edge2Cols;

  //优先去把col都重新整理下，成功整理到col_pool_in_memory里面，然后开始根据这个整理写map
  //当become为parent的时候一定是已经完全访问完成了
  void becomeParent(BBNODE *node, CVRP *cvrp);//当触发这个条件的时候，记得要把columns修改下

  void giveBirth(BidirLinkList *&lnode, BidirLinkList *&rnode);

  void deleteSelf();

  explicit BidirLinkList(BidirLinkList *pnode) {
    PNode = pnode;
  }

  ~BidirLinkList() = default;
};
//脚手架：可以先试验cut_index 和 col_index 和 一个col的边的unordered_map
//方案为：另外构建一个双向的二叉链表
//当子节点选择与父节点失去联系的时候，父节点同时也失去子节点的联系，同时父节点会检查自身是不是失去了所有的子节点，如果是就delete自己。
//同时子节点继承父节点的所有已有的值，失去联系的选择在于发生了列减少或者（和）行减少的时候
//程序退出的时候，记住了主动把所有的内存都delete掉

class BBNODE {
 public:

  SOLVER solver{};

  BidirLinkList *Ptr{};

  size_t *IdxCols{};

  int SizeEnuColPool{};
  int validSize{};

  RowVectorXT IdxColsInEnuColPool;
  RowVectorXd Cost4ColsInEnuColPool;

  std::deque<sparseRowMatrixXd> MatInEnu;
//  sparseRowMatrixXd MatInEnu;
//  std::vector<int> LPRow2MatRow;//use vector is fine!
//  std::unordered_map<int, int> MatRow2LPRow;// only keep the rows that are in the LP
//  int MaxUsedRowInMat{};

///  增加的这几个，需要考虑到 一点就是 需要写进去 子节点 中。

  //初始起就分配MAX_COL_INDEX这么些内存，这个是一一的映射，后面访问还要注意改内存，要知道哪些是在in memory里面取的
  //哪些是在in——pricing里面取的，父节点前的就是在memory里面的，其他的都是自己的
  int NumParentCols{};
  int NumParentColsInLPSols{};

  //下面两个变量不用特意去初始化
  std::vector<std::pair<size_t, double>> Idx4LPSolsInColPool;

  bool if_Int{};//解到cg收敛 int 一定是会被terminate的
  bool if_terminated{};//如果node 的Val已经大于ub了，那么就认为node已经被terminated了
  int TreeLevel{}, Idx{}, SizeAllocatedMem{};
  double Val{};
  int *EdgeTail{}, *EdgeHead{};
  double *EdgeVal{};
  int NumEdges{};
  int NumRowInBucketGraph{};
  double LastGap{1};
  bool *Deleted_ColsInEnuPool{};
//#ifdef NominalBranchingInEnu
//  bool *NoCompatibleColsInEnuPool{};
//#endif
  std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher>
      map_col_pool;// update when new mat is generated! //edge sense
//  std::vector<std::pair<int, double>> optSolInIdxCol; do not use this

  //cuts management:
  std::vector<RCC> RCCs;
  std::vector<R1C> R1Cs;
  std::vector<R1C_multi> R1Cs_multi;
  std::vector<BrC> BrCs;

#ifdef useM_dynamically
  double t4oneLP{};
  double geo_r_star{1};
  double c{1};
  //if 0,0 then, use the old ratio to do the estimation
  std::unordered_map<std::pair<int, int>, std::tuple<double, double, int>, PairHasher> obj_change;
  double l_r_ratio{1};
  void updateState(double new_value, double &old_value, int n);
  void calculateR_star(double lift, double &new_r_star, CVRP *cvrp);
#endif

//  std::vector<int> VBasis;
//  std::vector<int> CBasis;

  explicit BBNODE(int num, int p_col, CVRP *cvrp);

#ifdef readEnumerationTrees
  BBNODE(int num);
#endif

  BBNODE(BBNODE *node, BidirLinkList *ptr, int p_col, int idx, const BrC &bf, int NumBucketsPerVertex);

  BBNODE(BBNODE *node, int NumBucketsPerVertex, int NumCol, const bool *if_use_arc);

  BBNODE(BBNODE *node, int p_col, int idx, const BrC &bf);

  ~BBNODE();

  void allocateMem(int num);

  void freeMem() const;

  Bucket **AllForwardBuckets{};

  int NumForwardBucketArcs{};

  int NumForwardJumpArcs{};

  Bucket **AllBackwardBuckets{};

  int NumBackwardBucketArcs{};

  int NumBackwardJumpArcs{};

#ifdef NominalBranchingInEnu
  BBNODE(BBNODE *node, int p_col, const std::vector<int> &BranchingColSet, int idx, const NBrC &nbf);
  std::vector<NBrC> NBrCs;
#endif
};

#endif //CVRP_BBNODE_HPP
