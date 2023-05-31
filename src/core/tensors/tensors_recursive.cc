#include "tensors_recursive.hh"
#include "../math.hh"

namespace libcppe {
namespace tensors_recursive {
Eigen::VectorXd T_recursive(int k, const Eigen::Vector3d& Rij, double damping_factor) {
  auto Tk_coeffs = Tk_coefficients(k + 1);
  int x, y, z;
  Eigen::VectorXd Tk(multipole_components(k));
  int idx = 0;
  for (x = k; x > -1; x--) {
    for (y = k; y > -1; y--) {
      for (z = k; z > -1; z--) {
        if (x + y + z != k) continue;
        Tk[idx] = T(Rij, x, y, z, Tk_coeffs, damping_factor);
        idx++;
      }
    }
  }
  return Tk;
}

// x y z to index in contiguous array
int xyz2idx(int x, int y, int z) {
  int idx = 0;
  int k   = x + y + z;
  for (int a = k; a > -1; a--) {
    for (int b = k; b > -1; b--) {
      for (int c = k; c > -1; c--) {
        if (a + b + c != k) continue;
        if (a != x || b != y || c != z) {
          idx++;
        } else {
          return idx;
        }
      }
    }
  }
  return -1;
}

double T(const Eigen::Vector3d& Rij, int x, int y, int z,
         std::vector<Eigen::MatrixXi>& Cijn, double damping_factor) {
  double t = 0.0;
  double R = Rij.norm();
  double Cx, Cy, Cz;
  int k = x + y + z;

  // Thole Damping
  std::vector<double> scr_facs;
  if (damping_factor != 0.0) {
    if (k > 3) {
      throw std::runtime_error("Thole damping only implemented up to third order.");
    }
    scr_facs = thole_screening_factors(R * damping_factor, k);
  }

  for (size_t l = 0; l <= x; l++) {
    Cx = Cijn[0](x, l) * pow((Rij(0) / R), l);
    for (size_t m = 0; m <= y; m++) {
      Cy = Cx * Cijn[l + x](y, m) * pow((Rij(1) / R), m);
      for (size_t n = 0; n <= z; n++) {
        Cz     = Cy * Cijn[l + x + m + y](z, n) * pow((Rij(2) / R), n);
        int kk = l + m + n;
        // Thole damping
        if (scr_facs.size()) {
          if (kk == 0) {
            if (k == 0) {
              Cz *= scr_facs[0];
            } else if (k == 2) {
              Cz *= scr_facs[1];
            }
          } else if (kk == 1) {
            if (k == 1) {
              Cz *= scr_facs[1];
            } else if (k == 3) {
              Cz *= scr_facs[2];
            }
          } else if (kk == 2) {
            if (k == 2) Cz *= scr_facs[2];
          } else if (kk == 3) {
            if (k == 3) Cz *= scr_facs[3];
          }
        }
        t += Cz;
      }
    }
  }
  t /= pow(R, x + y + z + 1);
  return t;
}

std::vector<double> thole_screening_factors(double v, int k) {
  std::vector<double> ret;
  // screening factors for potential, field, field gradient, field Hessian
  double f_v, f_E, f_T, f_D = 1.0;
  if (k >= 0) {
    f_v = 1.0 - (0.5 * v + 1.0) * std::exp(-v);
    ret.push_back(f_v);
  }
  if (k >= 1) {
    f_E = f_v - (0.5 * std::pow(v, 2) + 0.5 * v) * std::exp(-v);
    ret.push_back(f_E);
  }
  if (k >= 2) {
    f_T = f_E - 1.0 / 6.0 * std::pow(v, 3) * std::exp(-v);
    ret.push_back(f_T);
  }
  if (k >= 3) {
    f_D = f_T - 1.0 / 30.0 * std::pow(v, 4) * std::exp(-v);
    ret.push_back(f_D);
  }
  return ret;
}

std::vector<Eigen::MatrixXi> Tk_coefficients(int max_order) {
  int maxi = 2 * max_order + 3;
  std::vector<Eigen::MatrixXi> Cijn;
  for (int n = 0; n < maxi; ++n) {
    int k;
    Eigen::MatrixXi mat = Eigen::MatrixXi::Zero(max_order + 2, max_order + 2);
    mat(0, 0)           = 1;
    if ((n + 1) % 2 == 0) {
      Cijn.push_back(mat);
      continue;
    }
    for (size_t i = 1; i <= max_order + 1; i++) {
      if (i % 2 != 0) {
        k = i - 1;
      } else if (i % 2 == 0) {
        k = i;
      }
      for (size_t j = 0; j <= i; j++) {
        if ((i + j) % 2 != 0) continue;
        if (j == 0) {
          mat(i, j) = mat(i - 1, j + 1);
        } else if (j != i) {
          mat(i, j) = (j + 1) * mat(i - 1, j + 1);
          mat(i, j) = mat(i, j) - ((n + 1) + k) * mat(i - 1, j - 1);
          k         = k + 2;
        } else if (j == i) {
          mat(i, j) = -((n + 1) + k) * mat(i - 1, j - 1);
        }
        mat(j, i) = mat(i, j);
      }
    }
    Cijn.push_back(mat);
  }
  return Cijn;
}

}  // namespace tensors_recursive
}  // namespace libcppe