#include "math.hh"

namespace libcppe {
  

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
  // std::cout << Cijn << std::endl;
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
  // changed signs here because electron charges are included downstream
  
  
  symmetry_factors(k, pf);
  for (size_t i = 0; i < pf.size(); i++) {
    pf[i] *= taylor;
  }
      
}

void prefactors_nuclei(unsigned k, std::vector<double> &pf) {
  double taylor;
  // changed signs here because electron charges are included downstream
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

} // namespace libcppe