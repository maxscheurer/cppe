#include "math.hh"

#include <cassert>

namespace libcppe {

// this only works for the contraction of polarizability/interaction tensors
// with vectors of size 3
// Eigen::Vector3d smat_vec(const Eigen::VectorXd& mat, const Eigen::Vector3d& vec,
//                          double alpha) {
//   Eigen::Vector3d result;
//   // expect upper triangle be provided
//   result[0] = mat[0] * vec[0] + mat[1] * vec[1] + mat[2] * vec[2];
//   result[1] = mat[1] * vec[0] + mat[3] * vec[1] + mat[4] * vec[2];
//   result[2] = mat[2] * vec[0] + mat[4] * vec[1] + mat[5] * vec[2];
//   result *= alpha;
//   return result;
// }

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

Eigen::VectorXd Tk_tensor(int k, const Eigen::Vector3d& Rij,
                          std::vector<Eigen::MatrixXi>& Tk_coeffs, double damping_factor,
                          double alpha_i, double alpha_j) {
  int x, y, z;
  Eigen::VectorXd Tk(multipole_components(k));
  int idx = 0;
  for (x = k; x > -1; x--) {
    for (y = k; y > -1; y--) {
      for (z = k; z > -1; z--) {
        if (x + y + z != k) continue;
        Tk[idx] = T(Rij, x, y, z, Tk_coeffs, damping_factor, alpha_i, alpha_j);
        idx++;
      }
    }
  }
  return Tk;
}

// TODO: there must be a way to make this more efficient...
// only 1st derivative supported
Eigen::VectorXd multipole_derivative(int k, int l, const Eigen::Vector3d& Rji,
                                     Eigen::VectorXd Mkj,
                                     std::vector<Eigen::MatrixXi>& Tk_coeffs,
                                     double damping_factor, double alpha_i,
                                     double alpha_j) {
  if (l > 1) throw std::runtime_error("Only 1st derivatives supported for multipoles");
  Eigen::VectorXd Fi = Eigen::VectorXd::Zero(3);

  double taylor;
  if ((k + l) % 2 == 0) {
    taylor = 1.0 / factorial(k);
  } else if ((k + l) % 2 != 0) {
    taylor = -1.0 / factorial(k);
  }

  int i, j, m, x, y, z;
  double symfac;
  // std::cout << "mul k = " << k << std::endl;
  // std::cout << "l = " << l << std::endl;
  Eigen::VectorXd Tk = Tk_tensor(k + l, Rji, Tk_coeffs, damping_factor, alpha_i, alpha_j);
  for (x = k + l; x > -1; x--) {
    for (y = k + l; y > -1; y--) {
      for (z = k + l; z > -1; z--) {
        if (x + y + z != k + l) continue;
        i = xyz2idx(x, y, z);
        for (int a = x; a > -1; a--) {
          for (int b = y; b > -1; b--) {
            for (int c = z; c > -1; c--) {
              if (a + b + c != k) continue;
              j      = xyz2idx(a, b, c);
              m      = xyz2idx(x - a, y - b, z - c);
              symfac = trinom(a, b, c);
              Fi(m) += taylor * symfac * Tk(i) * Mkj(j);
            }
          }
        }
      }
    }
  }
  return Fi;
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
         std::vector<Eigen::MatrixXi>& Cijn, double damping_factor, double alpha_i,
         double alpha_j) {
  double t = 0.0;
  double R = Rij.norm();
  double Cx, Cy, Cz;
  int k = x + y + z;

  // Thole Damping
  std::vector<double> scr_facs;
  if (damping_factor != 0.0) {
    if (alpha_i <= 0.0 || alpha_j <= 0.0) {
      throw std::runtime_error(
            "Thole damping only valid for non-zero"
            " isotropic polarizabilities.");
    }
    if (k > 3) {
      throw std::runtime_error("Thole damping only implemented up to third order.");
    }
    // Molecular Simulation, 32:6, 471-484, DOI: 10.1080/08927020600631270
    // v = factor * u, with u = R / (alpha_i * alpha_j)**(1/6)
    double v = damping_factor * R / std::pow(alpha_i * alpha_j, 1.0 / 6.0);
    scr_facs = thole_screening_factors(v, k);
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
