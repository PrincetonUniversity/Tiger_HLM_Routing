#ifndef RHS_HPP
#define RHS_HPP

#include <vector>
#include <algorithm> // for std::max
#include <cmath>     // for std::pow
#include <boost/numeric/odeint.hpp> // assuming you're using Boost odeint

using namespace boost::numeric::odeint;

/**
 * @brief Right-hand side function for the nonlinear routing ODE.
 * This function computes the rate of change of discharge (dQdt) based on the current discharge (q),
 * runoff, inflow from parent nodes, hillslope area, and routing parameters.
 * @param q Current discharge (m^3/s).
 * @param dQdt Rate of change of discharge (m^3/s^2).
 * @param t Current time (in minutes).
 * @note The function assumes that the inputs are provided in hourly resolution and converts them to the
 * appropriate units for the calculations. The runoff is converted from mm/h to m/min. (will change by adding index for resolution)
 */

struct RHS {
    const float* runoff_series; // runoff, hourly
    const std::vector<double>& y_p_series; // inflow from parent nodes, hourly (for now)
    const double A_h; // hillslope area in m^2
    const double lambda_1; // exponent
    const double invtau; // constant term

    RHS(const float* runoff_series,
        const std::vector<double>& y_p_series,
        const double A_h, 
        const double lambda_1, 
        const double invtau
    )
    : runoff_series(runoff_series),
      y_p_series(y_p_series),
      A_h(A_h),
      lambda_1(lambda_1),
      invtau(invtau) {}

    inline void operator()(const double q, double& dQdt, const double t) const {
        // t is in minutes, inputs are in hours 
        size_t hour_idx = static_cast<size_t>(t / 60);  // assumes t in minutes 

        // runoff_series is hourly 
        double runoff = runoff_series[hour_idx] * (0.001 / 60); // mm/h to m/min

        // y_p is hourly — 60 minutes per step
        double y_p = y_p_series[hour_idx];

        // Nonlinear routing ODE
        double q_safe = std::max(q, 1e-8); // avoid division by zero
        dQdt = invtau * pow(q_safe, lambda_1) * (-q_safe + (runoff * A_h / 60.0) + y_p);
    }
};

// RK45 Solver typedef
typedef runge_kutta_dopri5<
    double, double, double, double,
    vector_space_algebra, default_operations, never_resizer
> stepper_type;

#endif // RHS_HPP
