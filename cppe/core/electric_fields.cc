#include <Eigen/Dense>

#include "bmatrix.hh"
#include "electric_fields.hh"
#include "math.hh"
#include "tensors.hh"
#include <iomanip>

namespace libcppe {
using namespace tensors_recursive;

// TODO: there must be a way to make this more efficient...
// only 1st derivative supported
Eigen::VectorXd multipole_derivative(int k, int l, const Eigen::Vector3d& Rji,
                                     Eigen::VectorXd Mkj, double damping_factor) {
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

  Eigen::VectorXd Tk;
  if (damping_factor > 0.0) {
    Tk = tensors::T_damp_thole[k + l](Rji, damping_factor);
  } else {
    Tk = tensors::T[k + l](Rji);
  }
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

Eigen::VectorXd NuclearFields::compute() {
  Eigen::VectorXd nuc_fields = Eigen::VectorXd::Zero(3 * m_n_polsites);
#pragma omp parallel for
  for (size_t i = 0; i < m_n_polsites; i++) {
    size_t site_counter           = 3 * i;
    Potential& potential          = m_polsites[i];
    Eigen::Vector3d site_position = potential.get_site_position();
    for (auto& atom : m_mol) {
      Eigen::Vector3d core_position = atom.get_pos();
      Eigen::Vector3d diff          = site_position - core_position;
      Eigen::VectorXd Tms           = tensors::T1(diff);
      nuc_fields(site_counter) -= atom.charge * Tms(0);
      nuc_fields(site_counter + 1) -= atom.charge * Tms(1);
      nuc_fields(site_counter + 2) -= atom.charge * Tms(2);
    }
  }
  return nuc_fields;
}

Eigen::VectorXd MultipoleFields::compute() {
  Eigen::VectorXd mult_fields = Eigen::VectorXd::Zero(3 * m_n_polsites);
// Field at site of potential1 caused by all other sites (also non-polarizable
// sites!!!) size_t site_counter = 0;
#pragma omp parallel for
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
      for (auto& mul : potential2.get_multipoles()) {
        // if (std::all_of(mul.get_values().begin(), mul.get_values().end(),
        //                 [](double v) { return std::abs(v) == 0.0; })) {
        //   continue;
        // }
        Eigen::VectorXd Fi;
        // Molecular Simulation, 32:6, 471-484, DOI: 10.1080/08927020600631270
        double v = m_options.damping_factor_multipole /
                   std::pow(alpha_i_isotropic * alpha_j_isotropic, 1.0 / 6.0);
        if (damp_enabled) {
          Fi = multipole_derivative(mul.m_k, 1, diff, mul.get_values_vec(), v);
        } else {
          Fi = multipole_derivative(mul.m_k, 1, diff, mul.get_values_vec());
        }
        mult_fields(site_counter) += Fi(0);
        mult_fields(site_counter + 1) += Fi(1);
        mult_fields(site_counter + 2) += Fi(2);
      }
    }
  }
  return mult_fields;
}

Eigen::VectorXd InducedMoments::compute_cg(const Eigen::VectorXd& rhs) {
  // TODO: cleanup, print warning -> outside of solver
  // TODO: probably not entirely correct...
  BMatrix bmat(m_polsites, m_options);
  bool converged = false;

  Eigen::VectorXd x0 = bmat.compute_apply_diagonal(rhs);
  Eigen::VectorXd r0 = rhs - bmat.compute_apply(x0);
  Eigen::VectorXd z0 = bmat.compute_apply_diagonal(r0);
  Eigen::VectorXd p  = z0;

  Eigen::VectorXd x_k1, r_k1, z_k1;
  double alpha_k, beta_k;

  std::vector<Eigen::VectorXd> x{x0};
  std::vector<Eigen::VectorXd> r{r0};
  std::vector<Eigen::VectorXd> z{z0};
  for (int k = 0; k < m_options.maxiter; ++k) {
    Eigen::VectorXd Ap = bmat.compute_apply(p);
    alpha_k            = r[k].dot(z[k]) / p.dot(Ap);
    x_k1               = x[k] + alpha_k * p;
    x.push_back(x_k1);
    r_k1 = r[k] - alpha_k * Ap;
    r.push_back(r_k1);
    double rnorm = r_k1.norm();

    std::stringstream ss;
    ss.precision(12);
    ss << std::fixed << rnorm;
    m_printer(std::to_string(k) + " --- Norm: " + ss.str());
    if (rnorm < m_options.induced_thresh) {
      converged = true;
      break;
    }

    z_k1 = bmat.compute_apply_diagonal(r_k1);
    z.push_back(z_k1);
    beta_k = z_k1.dot(r_k1) / z[k].dot(r[k]);
    p      = z_k1 + beta_k * p;
  }
  if (!converged) {
    throw std::runtime_error("Failed to converge induced dipole moments.");
  }
  return x.back();
}

void InducedMoments::compute(const Eigen::VectorXd& total_fields,
                             Eigen::VectorXd& induced_moments, bool make_guess) {
  m_printer("        Running solver for induced moments.");
  BMatrix bmat{m_polsites, m_options};
  // guess
  if (make_guess) {
    induced_moments = bmat.compute_apply_diagonal(total_fields);
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

    induced_moments = bmat.compute_gauss_seidel_update(induced_moments, total_fields);
    norm            = (diis_old_moments - induced_moments).norm();

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
    throw std::runtime_error("Failed to converge induced dipole moments.");
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
