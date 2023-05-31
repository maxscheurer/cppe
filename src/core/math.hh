#pragma once

#include <Eigen/Core>
#include <vector>

namespace libcppe {

double factorial(int n);

void make_df(unsigned k, std::vector<double>& df);

int trinom(int i, int j, int k);

std::vector<double> symmetry_factors(unsigned k);

std::vector<double> prefactors(unsigned k);

std::vector<double> prefactors_nuclei(unsigned k);

int multipole_components(int k);

template <class T>
Eigen::Matrix3d triangle_to_mat(T);

Eigen::VectorXd mat_to_triangle(Eigen::Matrix3d);

}  // namespace libcppe
