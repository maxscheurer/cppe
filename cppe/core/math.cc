#include "math.hh"

#include <cassert>

namespace libcppe {
  
// this only works for the contraction of polarizability/interaction tensors with vectors of size 3
// TODO: need to make this faster
arma::vec smat_vec(arma::vec mat, arma::vec vec, double alpha) {
  assert(mat.n_elem == 6);
  assert(vec.n_elem == 3);
  arma::vec result(3, arma::fill::zeros);
  // expect upper triangle be provided
  result[0] = mat[0]*vec[0] + mat[1]*vec[1] + mat[2]*vec[2];
  result[1] = mat[1]*vec[0] + mat[3]*vec[1] + mat[4]*vec[2];
  result[2] = mat[2]*vec[0] + mat[4]*vec[1] + mat[5]*vec[2];  
  result *= alpha;
  return result;
}

// TODO: add option for damping
arma::vec Tk_tensor(int k, arma::vec Rij, arma::Cube<int>& Tk_coeffs) {
  int x,y,z;
  arma::vec Tk(multipole_components(k));
  int idx = 0;
  for (x = k; x > -1; x--) {
    for (y = k; y > -1; y--) {
      for (z = k; z > -1; z--) {
        if (x + y + z != k) continue;
        Tk[idx] = T(Rij, x,y,z, Tk_coeffs);
        idx++;
      }
    }
  }
  return Tk;
}

// TODO: there must be a way to make this more efficient...
// only 1st derivative supported
arma::vec multipole_derivative(int k, int l, arma::vec Rji, arma::vec Mkj, arma::Cube<int>& Tk_coeffs) {
  if (l > 1) throw std::runtime_error("Only 1st derivatives supported for multipoles");
  arma::vec Fi(3, arma::fill::zeros);
  
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
  arma::vec Tk = Tk_tensor(k + l, Rji, Tk_coeffs);
  for (x = k+l; x > -1; x--) {
    for (y = k+l; y > -1; y--) {
      for (z = k+l; z > -1; z--) {
        if (x + y + z != k+l) continue;
        i = xyz2idx(x, y, z);
        for (int a = x; a > -1; a--) {
          for (int b = y; b > -1; b--) {
            for (int c = z; c > -1; c--) {
              if (a + b + c != k) continue;
              j = xyz2idx(a, b, c);
              m = xyz2idx(x-a, y-b, z-c);
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
  int k = x + y + z;
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
}

double T(arma::vec Rij, int x, int y, int z, arma::Cube<int>& Cijn) {
  double t = 0.0;
  double R = arma::norm(Rij);
  double Cx, Cy, Cz;
  for (size_t l = 0; l <= x; l++) {
    Cx = Cijn.slice(0)(x,l) * pow((Rij(0) / R), l);
    for (size_t m = 0; m <= y; m++) {
      Cy = Cx * Cijn.slice(l+x)(y,m) * pow((Rij(1) / R), m);
      for (size_t n = 0; n <= z; n++) {
        Cz = Cy * Cijn.slice(l+x+m+y)(z,n) * pow((Rij(2) / R), n);
        t += Cz;
      }
    }
  }
  t /= pow(R, x+y+z+1);
  return t;
}


arma::Cube<int> Tk_coefficients(int max_order) {
  int maxi = 2*max_order+3;
  arma::Cube<int> Cijn(max_order+2, max_order+2, maxi, 
                  arma::fill::zeros);
  
  // TODO: I'm so bad with armadillo currently...
  for (int n = 0; n < maxi; ++n) {
    int k;
    arma::Mat<int>& mat = Cijn.slice(n);
    mat(0,0) = 1;
    if ((n+1) % 2 == 0) continue;
    for (size_t i = 1; i <= max_order+1; i++) {
      if (i % 2 != 0) {
        k = i - 1;
      } else if(i % 2 == 0) {
        k = i;
      }
      for (size_t j = 0; j <= i; j++) {
        if ((i+j) % 2 != 0) continue;
        if (j == 0) {
          mat(i,j) = mat(i-1,j+1);
        } else if (j != i) {
          mat(i,j) = (j + 1) * mat(i-1,j+1);
          mat(i,j) = mat(i,j) - ((n+1) + k) * mat(i-1,j-1);
          k = k + 2;
        } else if (j == i) {
          mat(i,j) = - ((n+1) + k) * mat(i-1,j-1);
        }
        mat(j,i) = mat(i,j);
      }
    }
  }
  return Cijn;
}

  
double factorial(int n) {
    if(n < 2) return 1.0;
    double x = 1.0;
    for(int i = 2; i <= n; i++) x *= double(i);
    return x;
}


void make_df(unsigned k, std::vector<double> &df) {

    double f = -1.0;
    for(unsigned i = 1; i <= k; i++) {
        f *= double(i) / double(2*i-1);
    }
    for(unsigned ix = 0, ixyz = 0; ix <= k; ix++)
    for(unsigned iy = 0; iy <= ix; iy++, ixyz++) {
        unsigned kx = k-ix, ky = ix-iy, kz = k-kx-ky;
        double g = factorial(kx) * factorial(ky) * factorial(kz);
        df[ixyz] = f / g;
    }
}

int trinom(int i, int j, int k) {
  return factorial(i+j+k) / (factorial(i) * factorial(j) * factorial(k));
}

void symmetry_factors(unsigned k, std::vector<double> &pf) {
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
}

void prefactors(unsigned k, std::vector<double> &pf) {
  double taylor = -1.0 / factorial(k);;
  // changed signs here because electron charges are included downstream (integral library)
  
  symmetry_factors(k, pf);
  for (size_t i = 0; i < pf.size(); i++) {
    pf[i] *= taylor;
  }
      
}

void prefactors_nuclei(unsigned k, std::vector<double> &pf) {
  double taylor;
  if (k % 2 == 0) {
    taylor = 1.0/factorial(k);
  } else if (k % 2 != 0) {
    taylor = -1.0/factorial(k);
  }
  
  symmetry_factors(k, pf);
  for (size_t i = 0; i < pf.size(); i++) {
    pf[i] *= taylor;
  }
      
}

int multipole_components(int k) {
  return (k+1)*(k+2)/2;
}


} // namespace libcppe