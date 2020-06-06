#pragma once

#include <Eigen/Core>
#include <functional>
#include <vector>

namespace libcppe {

namespace tensors {
Eigen::VectorXd T0(const Eigen::Vector3d&);
Eigen::VectorXd T1(const Eigen::Vector3d&);
Eigen::VectorXd T2(const Eigen::Vector3d&);
Eigen::VectorXd T3(const Eigen::Vector3d&);
Eigen::VectorXd T4(const Eigen::Vector3d&);
Eigen::VectorXd T5(const Eigen::Vector3d&);
Eigen::VectorXd T6(const Eigen::Vector3d&);

Eigen::VectorXd T0_damp_thole(const Eigen::Vector3d&, double a);
Eigen::VectorXd T1_damp_thole(const Eigen::Vector3d&, double a);
Eigen::VectorXd T2_damp_thole(const Eigen::Vector3d&, double a);
Eigen::VectorXd T3_damp_thole(const Eigen::Vector3d&, double a);
Eigen::VectorXd T4_damp_thole(const Eigen::Vector3d&, double a);
Eigen::VectorXd T5_damp_thole(const Eigen::Vector3d&, double a);
Eigen::VectorXd T6_damp_thole(const Eigen::Vector3d&, double a);

using Tfun      = std::function<Eigen::VectorXd(const Eigen::Vector3d&)>;
using Tfun_damp = std::function<Eigen::VectorXd(const Eigen::Vector3d&, double)>;

static const std::vector<Tfun> T = {
      Tfun{T0}, Tfun{T1}, Tfun{T2}, Tfun{T3}, Tfun{T4}, Tfun{T5}, Tfun{T6},
};

static const std::vector<Tfun_damp> T_damp_thole = {
      Tfun_damp{T0_damp_thole}, Tfun_damp{T1_damp_thole}, Tfun_damp{T2_damp_thole},
      Tfun_damp{T3_damp_thole}, Tfun_damp{T4_damp_thole}, Tfun_damp{T5_damp_thole},
      Tfun_damp{T6_damp_thole},
};

}  // namespace tensors

}  // namespace libcppe