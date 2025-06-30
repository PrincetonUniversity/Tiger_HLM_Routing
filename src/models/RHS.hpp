#ifndef RHS_HPP
#define RHS_HPP

#include <vector>
#include <algorithm>  // for std::max
#include <cmath>      // for std::pow
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

struct RHS {
    const std::vector<double>& runoff_series;
    const std::vector<double>& y_p_series;
    const double A_h;
    const double lambda_1;
    const double invtau;

    RHS(const std::vector<double>& runoff_series,
        const std::vector<double>& y_p_series,
        const double A_h,
        const double lambda_1,
        const double invtau
    );

    void operator()(const double q, double& dQdt, const double t) const;
};

// RK45 stepper type alias
typedef runge_kutta_dopri5<double, double, double, double, vector_space_algebra, default_operations, never_resizer> stepper_type;
extern stepper_type::controlled_type stepper;

#endif // RHS_HPP
