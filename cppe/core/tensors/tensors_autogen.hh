#pragma once

#include <Eigen/Core>
#include <functional>
#include <vector>

namespace libcppe {

namespace tensors {
Eigen::VectorXd T0(double x, double y, double z);
Eigen::VectorXd T1(double x, double y, double z);
Eigen::VectorXd T2(double x, double y, double z);
Eigen::VectorXd T3(double x, double y, double z);
Eigen::VectorXd T4(double x, double y, double z);
Eigen::VectorXd T5(double x, double y, double z);
Eigen::VectorXd T6(double x, double y, double z);

Eigen::VectorXd T0_damp_thole(double x, double y, double z, double a);
Eigen::VectorXd T1_damp_thole(double x, double y, double z, double a);
Eigen::VectorXd T2_damp_thole(double x, double y, double z, double a);
Eigen::VectorXd T3_damp_thole(double x, double y, double z, double a);
Eigen::VectorXd T4_damp_thole(double x, double y, double z, double a);
Eigen::VectorXd T5_damp_thole(double x, double y, double z, double a);
Eigen::VectorXd T6_damp_thole(double x, double y, double z, double a);

using Tfun      = std::function<Eigen::VectorXd(double, double, double)>;
using Tfun_damp = std::function<Eigen::VectorXd(double, double, double, double)>;

static const std::vector<Tfun> T = {
      Tfun{T0}, Tfun{T1}, Tfun{T2}, Tfun{T3}, Tfun{T4}, Tfun{T5}, Tfun{T6},
};

static const std::vector<Tfun_damp> T_damp_thole = {
      Tfun_damp{T0_damp_thole}, Tfun_damp{T1_damp_thole}, Tfun_damp{T2_damp_thole}, Tfun_damp{T3_damp_thole},
      Tfun_damp{T4_damp_thole}, Tfun_damp{T5_damp_thole}, Tfun_damp{T6_damp_thole},
};

}  // namespace tensors

}  // namespace libcppe