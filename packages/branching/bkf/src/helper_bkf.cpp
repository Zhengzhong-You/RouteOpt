/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <cmath>
#include <functional>

#include "helper_bkf.hpp"

#include "route_opt_macro.hpp"

#include "bkf_macro.hpp"

namespace RouteOpt::Branching::BKF {
    double getLBTk(double omega, double alpha, double B, double k) {
        return (omega + k) * std::pow(B, 1 / (1 - alpha / (k + 1)));
    }

    double getUBTk(double omega, double alpha, double B, double k) {
        double exponent = alpha / (1.0 - alpha) * std::log(B);
        double omega_term = boost::math::gamma_p(k, exponent);
        return (omega + k) * k * std::pow(B, 1 / (1 - alpha)) * std::pow(exponent, -k) * omega_term;
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

        int cnt = 0;
        boost::uintmax_t iter = MAX_ITERATION;
        for (const auto &range: search_ranges) {
            try {
                auto result = boost::math::tools::toms748_solve(
                    bound_equation,
                    range.first,
                    range.second,
                    boost::math::tools::eps_tolerance<double>(BOOST_TOLERANCE_BIT),
                    iter
                );
                if (cnt == 0) {
                    lb_k = static_cast<int>(std::ceil(result.second));
                } else {
                    ub_k = static_cast<int>(std::floor(result.first));
                }
            } catch ([[maybe_unused]] const std::exception &e) {
                if (cnt == 0) {
                    lb_k = 1;
                } else {
                    ub_k = static_cast<int>(est_m);
                }
            }
            ++cnt;
        }
    }

    void updateState(double new_value, double &old_value, int n) {
        if (n == 0) {
            old_value = new_value;
        } else {
            if (new_value / old_value < TOLERANCE) old_value /= static_cast<double>(n + 1) / n;
                //sometimes a no improve will happen
            else old_value = old_value * std::pow(new_value, 1.0 / (n + 1)) * std::pow(1.0 / old_value, 1.0 / (n + 1));
        }
    }

    void updateStateWithWeights(double new_value, double &old_value, int n) {
        if (n == 0) {
            old_value = new_value;
        } else {
            if (new_value / old_value < TOLERANCE) old_value /= static_cast<double>(n + 1) / n;
                //sometimes a no improve will happen
            else old_value = std::sqrt(old_value * new_value);
        }
    }

    void updateStateAverage(double new_value, double &old_value, int n) {
        old_value = (old_value * n + new_value) / (n + 1);
    }
}
