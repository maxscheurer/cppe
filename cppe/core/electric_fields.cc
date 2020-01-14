#include <Eigen/Dense>

#include "electric_fields.hh"
#include "math.hh"
#include <iomanip>

namespace libcppe {

Eigen::VectorXd NuclearFields::compute() {
  std::vector<Eigen::MatrixXi> Tk_coeffs = Tk_coefficients(5);
  Eigen::VectorXd nuc_fields             = Eigen::VectorXd::Zero(3 * m_n_polsites);
#pragma omp parallel for firstprivate(Tk_coeffs)
  for (size_t i = 0; i < m_n_polsites; i++) {
    size_t site_counter           = 3 * i;
    Potential& potential          = m_polsites[i];
    Eigen::Vector3d site_position = potential.get_site_position();
    for (auto& atom : m_mol) {
      Eigen::Vector3d core_position = atom.get_pos();
      Eigen::Vector3d diff          = site_position - core_position;
      Eigen::VectorXd Tms           = Tk_tensor(1, diff, Tk_coeffs);
      nuc_fields(site_counter) -= atom.charge * Tms(0);
      nuc_fields(site_counter + 1) -= atom.charge * Tms(1);
      nuc_fields(site_counter + 2) -= atom.charge * Tms(2);
    }
  }
  return nuc_fields;
}

Eigen::VectorXd MultipoleFields::compute() {
  std::vector<Eigen::MatrixXi> Tk_coeffs = Tk_coefficients(5);
  Eigen::VectorXd mult_fields            = Eigen::VectorXd::Zero(3 * m_n_polsites);
// Field at site of potential1 caused by all other sites (also non-polarizable
// sites!!!) size_t site_counter = 0;
#pragma omp parallel for firstprivate(Tk_coeffs)
  for (size_t i = 0; i < m_n_polsites; i++) {
    size_t site_counter   = 3 * i;
    Potential& potential1 = m_polsites[i];
    double alpha_i_isotropic =
          potential1.get_polarizability().get_isotropic_value();  // for damping
    for (size_t j = 0; j < m_potentials.size(); j++) {
      Potential& potential2 =
            m_potentials[j];  // all other multipoles create el. field at site i
      if (potential1.index == potential2.index) continue;
      if (potential1.excludes_site(potential2.index)) continue;
      Eigen::Vector3d diff =
            potential1.get_site_position() - potential2.get_site_position();

      // Thole damping: potential2 needs to be polarizable
      bool damp_enabled        = m_options.damp_multipole;
      double alpha_j_isotropic = 0.0;
      if (!potential2.is_polarizable()) {
        damp_enabled = false;
      } else {
        alpha_j_isotropic =
              potential2.get_polarizability().get_isotropic_value();  // for damping
      }

      // std::cout << "-- created by site " << potential2.index << std::endl;
      for (auto& mul : potential2.get_multipoles()) {
        // if (std::all_of(mul.get_values().begin(), mul.get_values().end(),
        //                 [](double v) { return std::abs(v) == 0.0; })) {
        //   continue;
        // }
        Eigen::VectorXd Fi;
        if (damp_enabled) {
          Fi = multipole_derivative(mul.m_k, 1, diff, mul.get_values_vec(), Tk_coeffs,
                                    m_options.damping_factor_multipole, alpha_i_isotropic,
                                    alpha_j_isotropic);
        } else {
          Fi = multipole_derivative(mul.m_k, 1, diff, mul.get_values_vec(), Tk_coeffs);
        }
        mult_fields(site_counter) += Fi(0);
        mult_fields(site_counter + 1) += Fi(1);
        mult_fields(site_counter + 2) += Fi(2);
      }
    }
  }
  return mult_fields;
}

void InducedMoments::compute(const Eigen::VectorXd& total_fields,
                             Eigen::VectorXd& induced_moments, bool make_guess) {
  m_printer("        Running solver for induced moments.");
  std::vector<Eigen::MatrixXi> Tk_coeffs = Tk_coefficients(5);
  // guess
  if (make_guess) {
    size_t site_counter = 0;
    for (auto& pot : m_potentials) {
      if (!pot.is_polarizable()) continue;
      induced_moments.segment<3>(site_counter) =
            pot.get_polarizability().get_matrix() * total_fields.segment<3>(site_counter);
      site_counter += 3;
    }
  }
  int max_iter           = m_options.maxiter;
  bool do_diis           = m_options.do_diis;
  double diis_start_norm = m_options.diis_start_norm;

  int iteration  = 0;
  bool converged = false;
  double norm    = 0.0;

  bool diis       = false;
  int diis_maxvec = 10;
  std::vector<Eigen::VectorXd> diis_prev_moments;
  std::vector<Eigen::VectorXd> diis_residuals;
  Eigen::VectorXd diis_old_moments = induced_moments;

  // iterations
  while (!converged) {
    if (iteration >= max_iter) break;
    if (norm <= diis_start_norm && iteration > 1 && !diis && do_diis) {
      m_printer("        --- Turning on DIIS. ---");
      diis = true;
    }

    norm = 0.0;
    // TODO: abstract matrix apply, generalized solver
#pragma omp parallel for reduction(+ : norm) firstprivate(Tk_coeffs)
    for (int i = 0; i < m_n_polsites; ++i) {
      Eigen::Vector3d Ftmp  = Eigen::Vector3d::Zero();
      Eigen::Vector3d M1tmp = Eigen::Vector3d::Zero();
      int l                 = i * 3;
      Potential& pot1       = m_polsites[i];
      for (int j = 0; j < m_n_polsites; ++j) {
        int m           = 3 * j;
        Potential& pot2 = m_polsites[j];
        if (pot1.excludes_site(pot2.index) || i == j) {
          continue;
        }
        Eigen::Vector3d diff = pot2.get_site_position() - pot1.get_site_position();
        Eigen::VectorXd T2;
        if (m_options.damp_induced) {
          Polarizability& alpha_i = pot1.get_polarizability();
          Polarizability& alpha_j = pot2.get_polarizability();
          T2 = Tk_tensor(2, diff, Tk_coeffs, m_options.damping_factor_induced,
                         alpha_i.get_isotropic_value(), alpha_j.get_isotropic_value());
        } else {
          T2 = Tk_tensor(2, diff, Tk_coeffs);
        }
        Eigen::Matrix3d T2m = triangle_to_mat(T2);
        Ftmp += T2m * induced_moments.segment<3>(m);
      }
      // keep value to calculate residual
      M1tmp = induced_moments.segment<3>(l);
      Ftmp += total_fields.segment<3>(l);
      induced_moments.segment<3>(l) = pot1.get_polarizability().get_matrix() * Ftmp;
      // Calculate the residual
      M1tmp = induced_moments.segment<3>(l) - M1tmp;
      norm += M1tmp.norm();
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
      int diis_size     = diis_residuals.size() + 1;
      Eigen::MatrixXd B = Eigen::MatrixXd::Zero(diis_size, diis_size);
      for (size_t i = 1; i < diis_size; i++) {
        B(i, 0) = -1.0;
        B(0, i) = -1.0;
      }
      for (size_t i = 1; i < diis_size; i++) {
        for (size_t j = 1; j < diis_size; j++) {
          B(i, j) = diis_residuals[i - 1].dot(diis_residuals[j - 1]);
          B(j, i) = B(i, j);
        }
      }
      Eigen::VectorXd rhs = Eigen::VectorXd::Zero(diis_size);
      rhs(0)              = -1.0;

      Eigen::VectorXd weights = B.colPivHouseholderQr().solve(rhs);
      induced_moments.fill(0.0);
      for (size_t i = 0; i < diis_size - 1; i++) {
        induced_moments += weights[i + 1] * diis_prev_moments[i];
      }
    }

    if (diis) {
      norm = (induced_moments - diis_old_moments).norm();
    }

    diis_old_moments = induced_moments;

    std::stringstream ss;
    ss.precision(12);
    ss << std::fixed << norm;
    m_printer(std::to_string(iteration) + "        --- Norm: " + ss.str());
    // calculate based on iteration
    if (norm < m_options.induced_thresh) converged = true;

    iteration++;
  }

  if (!converged) {
    throw std::runtime_error("Failed to converge induced moments.");
  }

  double nrm = 0.0;
  for (int j = 0; j < m_n_polsites; ++j) {
    int m = 3 * j;
    nrm   = (induced_moments.segment<3>(m)).norm();
    if (nrm > 1.0) {
      int site = m_polsites[j].index;
      m_printer("WARNING: Induced moment on site " + std::to_string(site) +
                " is greater than 1 a.u.!");
    }
  }
}

// returns a vector of potentials that have polarizabilities
std::vector<Potential> get_polarizable_sites(std::vector<Potential> potentials) {
  std::vector<Potential> result;
  std::copy_if(potentials.begin(), potentials.end(), std::back_inserter(result),
               [](Potential& p) { return p.is_polarizable(); });
  return result;
}

}  // namespace libcppe
