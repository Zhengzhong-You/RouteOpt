//
// Created by You, Zhengzhong on 5/6/24.
//


#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include "dynamics.hpp"
#include "branching_node.hpp"
#include "branching.hpp"

using namespace std;
using namespace chrono;
using namespace boost::math::tools;
using namespace std::placeholders;
using namespace boost::math;

CVRP *Dynamics::cvrp{};
BbNode *Dynamics::node{};
double Dynamics::alpha{};
double Dynamics::est_m{};
double Dynamics::f{};
double Dynamics::opt_k{};
double Dynamics::solve_cg_time{};
bool Dynamics::is_use_more_k{};
std::vector<std::pair<std::pair<double, int>, std::pair<double, int> > > Dynamics::r_star_depth{};
double Dynamics::solve_a_node_time{};
double Dynamics::t_for_one_lp_b4{};
double Dynamics::c_b4{};
int Dynamics::n_nodes_b4{};
double Dynamics::t_for_one_lp_enu{};
double Dynamics::c_enu{};
int Dynamics::n_nodes_enu{};
double Dynamics::alpha_heuristic{};
double Dynamics::est_m_heuristic{};
double Dynamics::t_for_one_heuristic_b4{};
double Dynamics::t_for_one_heuristic_enu{};
double Dynamics::r_best{1};

double getLBTk(double omega, double alpha, double B, double k);

double getUBTk(double omega, double alpha, double B, double k);

double getRealTk(double omega, double alpha, double B, double k);

void getRange(double omega,
              double alpha,
              double B,
              double target,
              double k1,
              double est_m,
              int &lb_k,
              int &ub_k);

void Dynamics::init(CVRP *pr_cvrp) {
    Dynamics::cvrp = pr_cvrp;
}

void Dynamics::updateNode(BbNode *pr_node) {
    node = pr_node;
}

void Dynamics::updateState(double new_value, double &old_value, int n) {
    if (n == 0) {
        old_value = new_value;
    } else {
        if (new_value / old_value < TOLERANCE) old_value /= double(n + 1) / n; //sometimes a no improve will happen
        else old_value = old_value * pow(new_value, 1.0 / (n + 1)) * pow(1.0 / old_value, 1.0 / (n + 1));
    }
}

void Dynamics::updateStateWithWeights(double new_value, double &old_value, int n) {
    if (n == 0) {
        old_value = new_value;
    } else {
        if (new_value / old_value < TOLERANCE) old_value /= double(n + 1) / n; //sometimes a no improve will happen
        else old_value = sqrt(old_value * new_value);
    }
}

void Dynamics::updateStateAverage(double new_value, double &old_value, int n) {
    old_value = (old_value * n + new_value) / (n + 1);
}

void Dynamics::calculateRStar(double lift) {
    if (lift < TOLERANCE)cout << "lift= " << lift << endl;
    lift = max(lift, TOLERANCE);
    if (node->getTreeLevel() == 0) throw runtime_error("calculateRStar: root node");
    auto &edge = node->getBrCs().back().edge;
    auto dir = node->getBrCs().back().br_dir;

    if (r_star_depth.size() <= node->getTreeLevel()) r_star_depth.resize(node->getTreeLevel() + 1);
    auto &r_star = r_star_depth[node->getTreeLevel()];

    auto &recordings = dir ? r_star.second.second : r_star.first.second;
    auto &increase = dir ? r_star.second.first : r_star.first.first;

    auto revised_lift = lift / (1 - alpha / (node->getDynamicKNode() + 1));
    updateState(revised_lift, increase, recordings);
    ++recordings;

    double new_r_star;
    if (r_star.second.second == 0 || r_star.first.second == 0) {
        new_r_star = numeric_limits<float>::max();
        for (auto &r: r_star_depth) {
            if (r.first.second == 0 || r.second.second == 0) continue;
            new_r_star = min(new_r_star, sqrt(r.first.first * r.second.first));
        }
        if (new_r_star == numeric_limits<float>::max()) new_r_star = revised_lift;
    } else {
        new_r_star = sqrt(r_star.first.first * r_star.second.first);
    }
    r_best = new_r_star * R_DISCOUNT;
}

void Dynamics::evaluateM1() {
    if (est_m == 0) {
        double max_num = 0.9 * Config::ML_BranchPhase0;
        est_m = min(FIX_M + 0., max_num);
        cout << "est_m for learning= " << est_m << endl;
        alpha = est_m / double(Config::ML_BranchPhase0);
        verbose_call(
            cout << "m= " << est_m << " alpha= " << alpha << endl;
        )
    }
}

#ifdef DYNAMICS_TYPE
void Dynamics::giveDynamicK(int &num, bool if_lp) {
  if(!if_lp) throw runtime_error("giveDynamicK in this dynamic type does not support: if_lp==false");
  if (node->getTreeLevel() == 0) {
	num = Config::ML_BranchPhase0;
	ml_call(MASTER_VALVE_ML, {
	  num = min(Config::ML_BranchPhase0, int(BaseBranching::current_branching_info.branch_pair.size()));
	  opt_k = num;
	}
	)
	verbose_call(cout << "we directly use the max n= " << num << endl;
	)
	return;
  }
  double beta = (BaseBranching::ub - node->getCurrentNodeVal() - f) / r_best;
  if (beta > 16) {
	num = Config::ML_BranchPhase0;
  } else {
	double B = DYNAMICS_TYPE == 1 ? pow(2, beta) : beta;
	cout << "B= " << B << endl;
	num = max(min(int(DYNAMICS_COEFFICIENT * B), int(Config::ML_BranchPhase0)), 1);
  }
  ml_call(MASTER_VALVE_ML,
		  num = min(num, int(BaseBranching::current_branching_info.branch_pair.size()));
  )
  opt_k = num;
  verbose_call(cout << "opt_k= " << opt_k << " in linear dynamics" << endl)
}
#else
void Dynamics::giveDynamicK(int &num, bool if_lp) {
    auto beg = high_resolution_clock::now();
    int max_num_heuristic = cvrp->getIfInEnuState() ? Config::BranPhase2InEnu : Config::BranPhase2;
    if (!if_lp) {
        double ratio = cvrp->getIfInEnuState()
                           ? (t_for_one_heuristic_enu / t_for_one_lp_enu)
                           : (t_for_one_heuristic_b4
                              / t_for_one_lp_b4);
        if (ratio < RATIO_HEURISTIC_LP) {
            num = max_num_heuristic;
            return;
        }
        est_m_heuristic = ABS_HEURISTIC_NO_MORE;
        alpha_heuristic = min(double(opt_k) / double(est_m), MAX_HEURISTIC_ALPHA);
        cout << "est_m_heuristic= " << est_m_heuristic << " alpha_heuristic= " << alpha_heuristic << endl;
    }
    if (node->getTreeLevel() == 0) {
        if (if_lp) {
            num = Config::ML_BranchPhase0;
            ml_call(MASTER_VALVE_ML,
                    num = min(num,
                        int(BaseBranching::current_branching_info.branch_pair.size()));
            )
        } else {
            num = min(max_num_heuristic, (int) BaseBranching::current_branching_info.branch_pair.size());
            num = min(num, ABS_HEURISTIC_NO_MORE);
        }
        if (if_lp) opt_k = num;
        verbose_call(cout << "we directly use the max n= " << num << endl;
        )
        return;
    }
    double beta;
    if (r_best < TOLERANCE) goto MAX_CANDIDATES;
    beta = max((BaseBranching::ub - node->getCurrentNodeVal() - f) / r_best, 0.);
    if (beta > 16) {
    MAX_CANDIDATES:
        if (if_lp) {
            num = (int) est_m;
            ml_call(MASTER_VALVE_ML,
                    num = min(num, int(BaseBranching::current_branching_info.branch_pair.size()));
            )
            opt_k = num;
        } else {
            num = min((int) est_m_heuristic, (int) BaseBranching::current_branching_info.branch_pair.size());
            num = min(num, ABS_HEURISTIC_NO_MORE);
        }
        verbose_call(cout << "we directly use the full m= " << num << endl;)
        return;
    }

    double omega;
    double t_enu = if_lp ? t_for_one_lp_enu : t_for_one_heuristic_enu;
    double t_b4 = if_lp ? t_for_one_lp_b4 : t_for_one_heuristic_b4;
    if (cvrp->getIfInEnuState()) {
        if (t_enu < TOLERANCE) goto MAX_CANDIDATES;
        omega = max(c_enu / t_enu, 1.);
    } else {
        if (t_b4 < TOLERANCE) goto MAX_CANDIDATES;
        omega = max(c_b4 / t_b4, 1.);
    }
    double B = pow(2, beta);

    double al = if_lp ? alpha : alpha_heuristic;
    double e_m = if_lp ? est_m : est_m_heuristic;

    double alpha_beta = al * beta * log(sqrt(2));
    double k1 = sqrt(alpha_beta * (alpha_beta + 2 * al + 2 * omega - 2)) + alpha_beta + al - 1;
    int k1_up = (int) max(min(ceil(k1), double(e_m)), 1.);
    int k1_down = (int) max(min(floor(k1), double(e_m)), 1.);
    verbose_call(cout << "alpha= " << al << " beta= " << beta << " omega= " << omega)

    if (k1_up == k1_down) {
        k1 = k1_up;
    } else {
        double k1_up_val = getLBTk(omega, al, B, k1_up);
        double k1_down_val = getLBTk(omega, al, B, k1_down);
        if (k1_up_val < k1_down_val) {
            k1 = k1_up;
        } else {
            k1 = k1_down;
        }
    }
    verbose_call(cout << " k1= " << k1 << endl)

    double ub_tk = getUBTk(omega, al, B, k1);
    int lb_k, ub_k;
    getRange(omega, al, B, ub_tk, k1, e_m, lb_k, ub_k);
    verbose_call(cout << "lb_k= " << lb_k << " ub_k= " << ub_k << endl)
    if (lb_k == ub_k) {
        num = lb_k;
    } else {
        double real_tk = getRealTk(omega, al, B, k1);
        getRange(omega, al, B, real_tk, k1, e_m, lb_k, ub_k);
        verbose_call(cout << "revise: lb_k= " << lb_k << " ub_k= " << ub_k << endl)
        if (lb_k == ub_k) {
            num = lb_k;
        } else {
            vector<double> tmp;
            for (int i = lb_k; i <= ub_k; ++i) {
                tmp.emplace_back(getRealTk(omega, al, B, i));
            }
            auto it = min_element(tmp.begin(), tmp.end()) - tmp.begin() + lb_k;
            num = (int) it;
        }
    }
    if (is_use_more_k) {
        num *= 2;
    }
    num = min(num, int(e_m));
    if (if_lp) {
        if (BaseBranching::ub - node->getCurrentNodeVal() < f && f > TOLERANCE)num = 2; //at least have 2 trials
        ml_call(MASTER_VALVE_ML,
                num = min(num, int(BaseBranching::current_branching_info.branch_pair.size()));
        )
    } else {
        num = min(num, int(BaseBranching::current_branching_info.branch_pair.size()));
    }
    if (!if_lp)num = min(num, ABS_HEURISTIC_NO_MORE);
    if (if_lp) opt_k = num;
    verbose_call(cout << "opt_k= " << num << endl)
    auto end = high_resolution_clock::now();
    cout << "giveDynamicK spent= " << duration<double>(end - beg).count() << "s" << endl;
}
#endif

void Dynamics::getAverageT4LPNHeuristic(double average_t, bool if_lp) {
    auto &input_t_b4 = if_lp ? t_for_one_lp_b4 : t_for_one_heuristic_b4;
    auto &input_t_enu = if_lp ? t_for_one_lp_enu : t_for_one_heuristic_enu;

    if (cvrp->getIfInEnuState()) {
        updateStateAverage(average_t, input_t_enu, n_nodes_enu);
    } else {
        updateStateAverage(average_t, input_t_b4, n_nodes_b4);
    }
}

void Dynamics::copyDynamicData4EachNode(const Dynamics &d, Dynamics &d2) {
    d2.k_node = d.k_node;
}

void Dynamics::calculateF(double eps, double node_value) {
    double old_f = f;
    if (node_value == 0) {
        if (eps < ENUMERATION_TIME_CONCERN_LIMIT) {
            double tmp_f = cvrp->getMaxEnumerationSuccessGap() * BaseBranching::ub;
            f = max(tmp_f, f);
        } else return;
    } else {
        if (eps < c_enu * MIP_TIME_ENU_TIME_RATIO || c_enu == 0) {
            double tmp_f = BaseBranching::ub - node_value;
            f = max(tmp_f, f);
        } else return;
    }
    if (f > old_f)cout << "f has been updated to " << f << endl;
}

double getLBTk(double omega, double alpha, double B, double k) {
    return (omega + k) * pow(B, 1 / (1 - alpha / (k + 1)));
}

double getUBTk(double omega, double alpha, double B, double k) {
    double exponent = alpha / (1.0 - alpha) * log(B);
    double omega_term = boost::math::gamma_p(k, exponent);
    return (omega + k) * k * pow(B, 1 / (1 - alpha)) * pow(exponent, -k) * omega_term;
}

double getRealTk(double omega, double alpha, double B, double k) {
    auto integrand = [&](double x) {
        return k * std::pow(x, k - 1) * std::pow(B, 1.0 / (alpha * x + 1 - alpha));
    };

    const int n = 100;

    double a = 0.0;
    double b = 1.0;
    double h = (b - a) / n;
    double integral = integrand(a) + integrand(b);

    for (int i = 1; i < n; i += 2) {
        integral += 4 * integrand(a + i * h);
    }
    for (int i = 2; i < n - 1; i += 2) {
        integral += 2 * integrand(a + i * h);
    }

    integral *= h / 3.0;
    return integral * (omega + k);
}

void getRange(double omega,
              double alpha,
              double B,
              double target,
              double k1,
              double est_m,
              int &lb_k,
              int &ub_k) {
    std::function<double(double)> bound_equation =
            [omega, alpha, B, target](double k) {
        return getLBTk(omega, alpha, B, k) - target;
    };

    std::array<std::pair<double, double>, 2> search_ranges = {
        std::make_pair(1., k1),
        std::make_pair(k1, est_m)
    };

    boost::uintmax_t max_iterations = 1000; // Maximum number of iterations

    int cnt = 0;
    for (const auto &range: search_ranges) {
        try {
            auto result = boost::math::tools::toms748_solve(
                bound_equation,
                range.first,
                range.second,
                boost::math::tools::eps_tolerance<double>(6),
                max_iterations
            );
            if (cnt == 0) {
                lb_k = int(ceil(result.second));
            } else {
                ub_k = int(floor(result.first));
            }
        } catch (const std::exception &e) {
            if (cnt == 0) {
                lb_k = 1;
            } else {
                ub_k = int(est_m);
            }
        }
        ++cnt;
    }
}
