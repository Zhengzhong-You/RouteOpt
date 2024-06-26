
#ifndef VRPTW_HPP_
#define VRPTW_HPP_
#include "CVRP.hpp"

#ifdef SOLVER_VRPTW

struct vrptwVertexInfo;

class VRPTW : public CVRP {
 public:

 public:
  explicit VRPTW(const InstanceData &instanceData);

  void checkSolutionFeasibleByCapacity(bool &feasible) override;

  void cleanColumnsCapInfeasible(BbNode *node) override;

  void popArcGraph(BbNode *node) override;

  double transformCost(double x) override;

  void setResourceInBucketGraph() override;

  void getLowerBoundofMinimumNumberCars() override;

  ~VRPTW() override;
};

struct double_same_Tolerance {
  bool operator()(const double &v1, const double &v2) const {
	return !(abs(v1 - v2) < TOLERANCE);
  }
};

struct vrptwVertexInfo {
  double e_time{};
  double l_time{};
  double s_time{};
  vrptwVertexInfo(double et, double lt, double st) : e_time(et), l_time(lt), s_time(st) {}
  vrptwVertexInfo() = default;
};
#endif
#endif //VRPTW_HPP_
