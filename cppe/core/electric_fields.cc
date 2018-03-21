#include <iomanip>

#include "electric_fields.hh"
#include "math.hh"

namespace libcppe {
  
void NuclearFields::compute(arma::vec& nuc_fields, bool damp_core) {
  if (damp_core) {
    throw std::runtime_error("damping not implemented");
  }
  arma::Cube<int> Tk_coeffs = Tk_coefficients(5);
  size_t site_counter = 0;
  // TODO: this could be parallelized
  for (auto& potential : m_potentials) {
    if (!potential.is_polarizable()) continue;
    arma::vec site_position = potential.get_site_position();
    for (auto& atom : m_mol) {
      arma::vec core_position = atom.get_pos();
      arma::vec diff = site_position-core_position;
      arma::vec Tms = Tk_tensor(1, diff, Tk_coeffs);
      nuc_fields(site_counter) -=  atom.charge * Tms(0);
      nuc_fields(site_counter+1) -=  atom.charge * Tms(1);
      nuc_fields(site_counter+2) -=  atom.charge * Tms(2);
    }
    site_counter += 3;
  }
}

void MultipoleFields::compute(arma::vec& mult_fields, bool damp) {
  if (damp) {
    throw std::runtime_error("damping not implemented");
  }
  arma::Cube<int> Tk_coeffs = Tk_coefficients(5);
  // Field at site of potential1 caused by all other sites
  size_t site_counter = 0;
  for (auto& potential1 : m_potentials) {
    if (!potential1.is_polarizable()) continue;
    
    // std::cout << "Calculating field on site " << potential1.index << std::endl;
    for (auto& potential2 : m_potentials) {
      if (potential1.index == potential2.index) continue;
      if (potential1.excludes_site(potential2.index)) continue;
      arma::vec diff = potential1.get_site_position()-potential2.get_site_position();
      // std::cout << "-- created by site " << potential2.index << std::endl;
      for (auto& mul : potential2.get_multipoles()) {
        // TODO: exclude zero value multipoles here...
        // int non_zeros = std::count_if( mul.get_values().begin(), mul.get_values().end(), [](double val){return abs(val) > 0.0;} );
        // std::cout << "non-zeros: " << non_zeros << std::endl; 
        // if (non_zeros == 0) continue;
        arma::vec Fi = multipole_derivative(mul.m_k, 1, diff, mul.get_values(), Tk_coeffs);
        mult_fields(site_counter) += Fi(0);
        mult_fields(site_counter+1) += Fi(1);
        mult_fields(site_counter+2) += Fi(2);
      }
    }
    site_counter += 3;
  }
}


void InducedMoments::compute(arma::vec& total_fields, arma::vec& induced_moments, bool make_guess) {
  std::cout << "run induced moments" << std::endl;
  arma::Cube<int> Tk_coeffs = Tk_coefficients(5);
  // guess
  if (make_guess) {
    size_t site_counter = 0;
    for (auto& pot : m_potentials) {
      if(!pot.is_polarizable()) continue;
      arma::vec res = smat_vec(pot.get_polarizabilities()[0].get_values(), total_fields.subvec(site_counter, site_counter+2), 1.0);
      induced_moments.subvec(site_counter, site_counter+2) = res;
      site_counter += 3;
    }
  }
  // std::cout << "induced mom. guess" << std::endl;
  // induced_moments.raw_print(std::cout << std::setprecision(10));
  int thresh = 50;
  int iteration = 0;
  bool converged = false;
  double norm;
  
  int l, m;
  arma::vec Ftmp(3, arma::fill::zeros);
  arma::vec M1tmp(3, arma::fill::zeros);
  
  // iterations
  while (!converged) {
    if (iteration >= thresh) break;
    norm = 0.0;
    l = 0;
    for (auto& pot1 : m_potentials) {
      if(!pot1.is_polarizable()) continue;
      m = 0;
      Ftmp.fill(0.0);
      for (auto& pot2 : m_potentials) {
        if(!pot2.is_polarizable()) continue;
        if (pot1.excludes_site(pot2.index) || pot1.index == pot2.index) {
          m += 3;
          continue;
        }
        arma::vec diff = pot2.get_site_position()-pot1.get_site_position();
        arma::vec T2 = Tk_tensor(2, diff, Tk_coeffs);
        // += ???
        Ftmp += smat_vec(T2, induced_moments.subvec(m, m+2), 1.0);
        m += 3;
      }
      // keep value to calculate residual
      M1tmp = induced_moments.subvec(l, l+2);
      Ftmp += total_fields.subvec(l, l+2);
      induced_moments.subvec(l, l+2) = smat_vec(pot1.get_polarizabilities()[0].get_values(), Ftmp, 1.0);
      M1tmp = induced_moments.subvec(l, l+2) - M1tmp;
      norm += arma::norm(M1tmp);
      l += 3;
    }
    std::cout << iteration << " -- Norm: " << norm << std::endl;
    // calculate based on iteration
    if (norm < 1e-8) converged = true;
    iteration++;
  }
  
  if (!converged) {
    throw std::runtime_error("Failed to converge induced moments.");
  }
  // TODO: warning if induced moments > 1 a.u.
}
  

} // namespace libcppe