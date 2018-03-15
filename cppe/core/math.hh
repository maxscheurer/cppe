#ifndef INCLUDE_LIBCPPE_CPPE_CORE_MATH_H
#define INCLUDE_LIBCPPE_CPPE_CORE_MATH_H

#include <armadillo>

namespace libcppe {
  

double T(arma::vec Rij, int x, int y, int z, arma::Cube<int>& Cijn);

arma::Cube<int> Tk_coefficients(int max_order);
  
double factorial(int n);

void make_df(unsigned k, std::vector<double> &df);

int trinom(int i, int j, int k);

void symmetry_factors(unsigned k, std::vector<double> &pf);

void prefactors(unsigned k, std::vector<double> &pf);

void prefactors_nuclei(unsigned k, std::vector<double> &pf);

} // namespace libcppe

#endif //INCLUDE_LIBCPPE_CPPE_CORE_MATH_H