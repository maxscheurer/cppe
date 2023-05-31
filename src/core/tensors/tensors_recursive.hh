#pragma once
#include <Eigen/Core>
#include <vector>

namespace libcppe {
namespace tensors_recursive {
double T(const Eigen::Vector3d& Rij, int x, int y, int z,
         std::vector<Eigen::MatrixXi>& Cijn, double damping_factor = 0.0);

Eigen::VectorXd T_recursive(int k, const Eigen::Vector3d& Rij,
                            double damping_factor = 0.0);

std::vector<double> thole_screening_factors(double v, int k);

std::vector<Eigen::MatrixXi> Tk_coefficients(int max_order);

int xyz2idx(int x, int y, int z);
}  // namespace tensors_recursive
}  // namespace libcppe