#include "RHS.hpp"

RHS::RHS(const std::vector<double>& runoff_series,
         const std::vector<double>& y_p_series,
         const double A_h,
         const double lambda_1,
         const double invtau)
    : runoff_series(runoff_series),
      y_p_series(y_p_series),
      A_h(A_h),
      lambda_1(lambda_1),
      invtau(invtau) {}

void RHS::operator()(const double q, double& dQdt, const double t) const {
    size_t idx = static_cast<size_t>(t);
    double runoff = (idx < runoff_series.size()) ? runoff_series[idx] : 0.0;

    size_t hour_idx = idx / 60;
    double y_p = (hour_idx < y_p_series.size()) ? y_p_series[hour_idx] : 0.0;

    double q_safe = std::max(q, 1e-8);
    dQdt = invtau * std::pow(q_safe, lambda_1) * (-q + (runoff * A_h / 60.0) + y_p);
}

// Define the global stepper instance
stepper_type::controlled_type stepper = make_controlled(1E-9, 1E-6, stepper_type());
