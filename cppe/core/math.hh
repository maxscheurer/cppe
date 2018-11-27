#ifndef INCLUDE_LIBCPPE_CPPE_CORE_MATH_H
#define INCLUDE_LIBCPPE_CPPE_CORE_MATH_H

#include <armadillo>

namespace libcppe {

arma::vec smat_vec(arma::vec mat, arma::vec vec, double alpha);

arma::vec multipole_derivative(int k, int l, arma::vec Rji, arma::vec Mkj,
                               arma::Cube<int> &Tk_coeffs);

double T(arma::vec Rij, int x, int y, int z, arma::Cube<int> &Cijn);

arma::vec Tk_tensor(int k, arma::vec Rij, arma::Cube<int> &Tk_coeffs);

arma::Cube<int> Tk_coefficients(int max_order);

int xyz2idx(int x, int y, int z);

double factorial(int n);

void make_df(unsigned k, std::vector<double> &df);

int trinom(int i, int j, int k);

void symmetry_factors(unsigned k, std::vector<double> &pf);

void prefactors(unsigned k, std::vector<double> &pf);

void prefactors_nuclei(unsigned k, std::vector<double> &pf);

int multipole_components(int k);

} // namespace libcppe

#endif // INCLUDE_LIBCPPE_CPPE_CORE_MATH_H