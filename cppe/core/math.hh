#pragma once

#include <Eigen/Core>
#include <vector>

namespace libcppe {

Eigen::Vector3d smat_vec(Eigen::VectorXd mat, Eigen::Vector3d vec,
                         double alpha);

Eigen::VectorXd multipole_derivative(int k, int l, Eigen::Vector3d Rji,
                                     Eigen::VectorXd Mkj,
                                     std::vector<Eigen::MatrixXi> &Tk_coeffs);

double T(Eigen::Vector3d Rij, int x, int y, int z,
         std::vector<Eigen::MatrixXi> &Cijn);

Eigen::VectorXd Tk_tensor(int k, Eigen::Vector3d Rij,
                          std::vector<Eigen::MatrixXi> &Tk_coeffs);

std::vector<Eigen::MatrixXi> Tk_coefficients(int max_order);

int xyz2idx(int x, int y, int z);

// TODO: return type should be unsigned long
double factorial(int n);

void make_df(unsigned k, std::vector<double> &df);

int trinom(int i, int j, int k);

std::vector<double> symmetry_factors(unsigned k);

std::vector<double> prefactors(unsigned k);

std::vector<double> prefactors_nuclei(unsigned k);

int multipole_components(int k);

}  // namespace libcppe