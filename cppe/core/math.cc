#include "math.hh"

#include <cassert>

namespace libcppe {

template <class T>
Eigen::Matrix3d triangle_to_mat(T y) {
  Eigen::Matrix3d ret;
  ret(0, 0) = y[0];
  ret(0, 1) = ret(1, 0) = y[1];
  ret(0, 2) = ret(2, 0) = y[2];
  ret(1, 1)             = y[3];
  ret(1, 2) = ret(2, 1) = y[4];
  ret(2, 2)             = y[5];
  return ret;
}

template Eigen::Matrix3d triangle_to_mat(std::vector<double>);
template Eigen::Matrix3d triangle_to_mat(Eigen::VectorXd);

Eigen::VectorXd mat_to_triangle(Eigen::Matrix3d m) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(6);
  ret[0]              = m(0, 0);
  ret[1]              = m(0, 1);
  ret[2]              = m(0, 2);
  ret[3]              = m(1, 1);
  ret[4]              = m(1, 2);
  ret[5]              = m(2, 2);
  return ret;
}

double factorial(int n) {
  if (n < 2) return 1.0;
  double x = 1.0;
  for (int i = 2; i <= n; i++) x *= double(i);
  return x;
}

void make_df(unsigned k, std::vector<double>& df) {
  double f = -1.0;
  for (unsigned i = 1; i <= k; i++) {
    f *= double(i) / double(2 * i - 1);
  }
  for (unsigned ix = 0, ixyz = 0; ix <= k; ix++)
    for (unsigned iy = 0; iy <= ix; iy++, ixyz++) {
      unsigned kx = k - ix, ky = ix - iy, kz = k - kx - ky;
      double g = factorial(kx) * factorial(ky) * factorial(kz);
      df[ixyz] = f / g;
    }
}

int trinom(int i, int j, int k) {
  return factorial(i + j + k) / (factorial(i) * factorial(j) * factorial(k));
}

std::vector<double> symmetry_factors(unsigned k) {
  std::vector<double> pf(multipole_components(k));
  int x, y, z, idx;
  idx = 0;
  for (x = k; x > -1; x--) {
    for (y = k; y > -1; y--) {
      for (z = k; z > -1; z--) {
        if (x + y + z != k) continue;
        pf[idx] = static_cast<double>(trinom(x, y, z));
        idx++;
      }
    }
  }
  return pf;
}

std::vector<double> prefactors(unsigned k) {
  double taylor = -1.0 / factorial(k);
  // changed signs here because electron charges are included downstream
  // (integral library)
  std::vector<double> pref = symmetry_factors(k);
  for (size_t i = 0; i < pref.size(); i++) {
    pref[i] *= taylor;
  }
  return pref;
}

std::vector<double> prefactors_nuclei(unsigned k) {
  double taylor;
  if (k % 2 == 0) {
    taylor = 1.0 / factorial(k);
  } else if (k % 2 != 0) {
    taylor = -1.0 / factorial(k);
  }

  std::vector<double> pf = symmetry_factors(k);
  for (size_t i = 0; i < pf.size(); i++) {
    pf[i] *= taylor;
  }
  return pf;
}

int multipole_components(int k) { return (k + 1) * (k + 2) / 2; }

}  // namespace libcppe
