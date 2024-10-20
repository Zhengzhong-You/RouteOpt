#ifndef VRPTW_HPP_
#define VRPTW_HPP_
#include "cvrp.hpp"

class VRPTW : public CVRP {
public:
    std::vector<std::vector<double> > travel_time_matrix{};
#if  SOLVER_VRPTW == 1
 bool if_force_keep_rcc{};
#endif

public:
    explicit VRPTW(const InstanceData &instanceData);

    void tryGetTravelTimeMatrix(const std::string &file_name);

    void checkSolutionFeasibleByCapacity(bool &feasible) override;

    void cleanColumnsCapInfeasible(BbNode *node) override;

    void popArcGraph(BbNode *node) override;

    double transformCost(double x) override;

    void setResourceInBucketGraph() override;

    void getLowerBoundofMinimumNumberCars() override;

    ~VRPTW() override;
};
#endif //VRPTW_HPP_
