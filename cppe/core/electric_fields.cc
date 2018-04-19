#include <iomanip>

#include "electric_fields.hh"
#include "math.hh"

namespace libcppe {

void NuclearFields::compute(arma::vec& nuc_fields, bool damp_core) {
  if (damp_core) {
    throw std::runtime_error("damping not implemented");
  }
  arma::Cube<int> Tk_coeffs = Tk_coefficients(5);

  #pragma omp parallel for firstprivate(Tk_coeffs)
  for (size_t i = 0; i < m_n_polsites; i++) {
    size_t site_counter = 3*i;
    Potential& potential = m_polsites[i];
    arma::vec site_position = potential.get_site_position();
    for (auto& atom : m_mol) {
      arma::vec core_position = atom.get_pos();
      arma::vec diff = site_position-core_position;
      arma::vec Tms = Tk_tensor(1, diff, Tk_coeffs);
      nuc_fields(site_counter) -=  atom.charge * Tms(0);
      nuc_fields(site_counter+1) -=  atom.charge * Tms(1);
      nuc_fields(site_counter+2) -=  atom.charge * Tms(2);
    }
  }
  // for (auto& potential : m_potentials) {
  //   if (!potential.is_polarizable()) continue;
  //   arma::vec site_position = potential.get_site_position();
  //   for (auto& atom : m_mol) {
  //     arma::vec core_position = atom.get_pos();
  //     arma::vec diff = site_position-core_position;
  //     arma::vec Tms = Tk_tensor(1, diff, Tk_coeffs);
  //     nuc_fields(site_counter) -=  atom.charge * Tms(0);
  //     nuc_fields(site_counter+1) -=  atom.charge * Tms(1);
  //     nuc_fields(site_counter+2) -=  atom.charge * Tms(2);
  //   }
  //   site_counter += 3;
  // }
}

void MultipoleFields::compute(arma::vec& mult_fields, bool damp) {
  if (damp) {
    throw std::runtime_error("damping not implemented");
  }
  arma::Cube<int> Tk_coeffs = Tk_coefficients(5);
  // Field at site of potential1 caused by all other sites
  // size_t site_counter = 0;
  #pragma omp parallel for firstprivate(Tk_coeffs)
  for (size_t i = 0; i < m_n_polsites; i++) {
    size_t site_counter = 3*i;
    Potential& potential1 = m_polsites[i];
    for (size_t j = 0; j < m_n_polsites; j++) {
      if (i == j) continue;
      Potential& potential2 = m_polsites[j];
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
  }
  // for (auto& potential1 : m_potentials) {
  //   if (!potential1.is_polarizable()) continue;
  //
  //   // std::cout << "Calculating field on site " << potential1.index << std::endl;
  //   for (auto& potential2 : m_potentials) {
  //     if (potential1.index == potential2.index) continue;
  //     if (potential1.excludes_site(potential2.index)) continue;
  //     arma::vec diff = potential1.get_site_position()-potential2.get_site_position();
  //     // std::cout << "-- created by site " << potential2.index << std::endl;
  //     for (auto& mul : potential2.get_multipoles()) {
  //       // TODO: exclude zero value multipoles here...
  //       // int non_zeros = std::count_if( mul.get_values().begin(), mul.get_values().end(), [](double val){return abs(val) > 0.0;} );
  //       // std::cout << "non-zeros: " << non_zeros << std::endl;
  //       // if (non_zeros == 0) continue;
  //       arma::vec Fi = multipole_derivative(mul.m_k, 1, diff, mul.get_values(), Tk_coeffs);
  //       mult_fields(site_counter) += Fi(0);
  //       mult_fields(site_counter+1) += Fi(1);
  //       mult_fields(site_counter+2) += Fi(2);
  //     }
  //   }
  //   site_counter += 3;
  // }
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
  int max_iter = 50;
  double norm_thresh = 1e-8;
  int iteration = 0;
  bool converged = false;
  double norm = 0.0;

  // int l, m;
  // arma::vec Ftmp(3, arma::fill::zeros);
  // arma::vec M1tmp(3, arma::fill::zeros);

  bool diis = false;
  int diis_maxvec = 10;
  std::vector<arma::vec> diis_prev_moments;
  std::vector<arma::vec> diis_residuals;
  arma::vec diis_old_moments = induced_moments;


  // iterations
  while (!converged) {
    if (iteration >= max_iter) break;
    if (norm <= 1.0 && iteration > 1 && !diis) {
      std::cout << "--- Turning on DIIS. ---" << std::endl;
      diis = true;
    }

    norm = 0.0;
    #pragma omp parallel for reduction(+:norm)
    for (int i = 0; i < m_n_polsites; ++i) {
      arma::vec Ftmp(3, arma::fill::zeros);
      arma::vec M1tmp(3);
      int l = i*3;
      Potential& pot1 = m_polsites[i];
      for (int j = 0; j < m_n_polsites; ++j) {
        int m = 3*j;
        Potential& pot2 = m_polsites[j];
        if (pot1.excludes_site(pot2.index) || i == j) {
          continue;
        }
        arma::vec diff = pot2.get_site_position()-pot1.get_site_position();
        arma::vec T2 = Tk_tensor(2, diff, Tk_coeffs);
        Ftmp += smat_vec(T2, induced_moments.subvec(m, m+2), 1.0);
      }
      // keep value to calculate residual
      M1tmp = induced_moments.subvec(l, l+2);
      Ftmp += total_fields.subvec(l, l+2);
      induced_moments.subvec(l, l+2) = smat_vec(pot1.get_polarizabilities()[0].get_values(), Ftmp, 1.0);
      // Calculate the residual
      M1tmp = induced_moments.subvec(l, l+2) - M1tmp;
      norm += arma::norm(M1tmp);
    }

    diis_prev_moments.push_back(induced_moments);
    if (diis_prev_moments.size() > diis_maxvec) {
      diis_prev_moments.erase(diis_prev_moments.begin());
    }
    diis_residuals.push_back(induced_moments - diis_old_moments);
    if (diis_residuals.size() > diis_maxvec) {
      diis_residuals.erase(diis_residuals.begin());
    }

    if (diis_residuals.size() > 2 && diis) {
      int diis_size = diis_residuals.size() + 1;
      arma::mat B(diis_size, diis_size, arma::fill::zeros);
      for (size_t i = 1; i < diis_size; i++) {
        B(i,0) = -1.0;
        B(0,i) = -1.0;
      }
      for (size_t i = 1; i < diis_size; i++) {
        for (size_t j = 1; j < diis_size; j++) {
          B(i,j) = arma::dot(diis_residuals[i-1], diis_residuals[j-1]);
          B(j,i) = B(i,j);
        }
      }
      // std::cout << "B-matrix" << std::endl;
      // std::cout << B << std::endl;
      arma::vec rhs(diis_size, arma::fill::zeros);
      rhs(0) = -1.0;

      arma::vec weights = arma::solve(B, rhs);
      // std::cout << weights << std::endl;
      induced_moments.fill(0.0);
      for (size_t i = 0; i < diis_size - 1; i++) {
        induced_moments += weights[i+1] * diis_prev_moments[i];
      }
    }

    if (diis) {
      norm = arma::norm(induced_moments - diis_old_moments);
    }

    diis_old_moments = induced_moments;

    std::cout << iteration << std::setprecision(12) << "--- Norm: " << norm << std::endl;
    // calculate based on iteration
    if (norm < norm_thresh) converged = true;

    iteration++;
  }

  if (!converged) {
    throw std::runtime_error("Failed to converge induced moments.");
  }

  double nrm = 0.0;
  for (int j = 0; j < m_n_polsites; ++j) {
    int m = 3*j;
    nrm = arma::norm(induced_moments.subvec(m, m+2));
    if (nrm > 1.0) {
      int site = m_polsites[j].index;
      std::cout << "WARNING: Induced moment on site " << site << " is greater than 1 a.u.!" << std::endl;
    }
  }
}

// returns a vector of potentials that are polarizable
std::vector<Potential> get_polarizable_sites(std::vector<Potential> potentials) {
  std::vector<Potential> result;
  for (auto p : potentials) {
    if (p.is_polarizable()) {
      result.push_back(p);
    }
  }
  return result;
}


} // namespace libcppe
