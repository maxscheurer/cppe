#include <Eigen/Dense>

#include "bmatrix.hh"
#include "electric_fields.hh"
#include "fmm/tree.hh"
#include "math.hh"
#include "tensors.hh"
#include <iomanip>

namespace libcppe {
using namespace tensors_recursive;

// TODO: probably retire at some point...
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
      Eigen::Vector3d core_position = atom.get_position();
      Eigen::Vector3d diff          = site_position - core_position;
      Eigen::VectorXd Tms           = tensors::T1(diff);
      nuc_fields(site_counter) -= atom.charge * Tms(0);
      nuc_fields(site_counter + 1) -= atom.charge * Tms(1);
      nuc_fields(site_counter + 2) -= atom.charge * Tms(2);
    }
  }
  return nuc_fields;
}

Eigen::MatrixXd NuclearFields::nuclear_gradient() {
  int natoms           = m_mol.size();
  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(3 * natoms, 3 * m_n_polsites);
#pragma omp parallel for
  for (size_t i = 0; i < m_n_polsites; i++) {
    size_t site_counter           = 3 * i;
    Potential& potential          = m_polsites[i];
    Eigen::Vector3d site_position = potential.get_site_position();
    for (int ai = 0; ai < natoms; ++ai) {
      auto& atom                         = m_mol[ai];
      Eigen::Vector3d core_position      = atom.get_position();
      Eigen::Vector3d diff               = site_position - core_position;
      Eigen::VectorXd Tms                = tensors::T2(diff);
      grad(3 * ai + 0, site_counter + 0) = atom.charge * Tms(0);
      grad(3 * ai + 0, site_counter + 1) = atom.charge * Tms(1);
      grad(3 * ai + 0, site_counter + 2) = atom.charge * Tms(2);

      grad(3 * ai + 1, site_counter + 0) = atom.charge * Tms(1);
      grad(3 * ai + 1, site_counter + 1) = atom.charge * Tms(3);
      grad(3 * ai + 1, site_counter + 2) = atom.charge * Tms(4);

      grad(3 * ai + 2, site_counter + 0) = atom.charge * Tms(2);
      grad(3 * ai + 2, site_counter + 1) = atom.charge * Tms(4);
      grad(3 * ai + 2, site_counter + 2) = atom.charge * Tms(5);
    }
  }
  return grad;
}

Eigen::VectorXd MultipoleFields::compute() {
  int n_sites  = m_potentials.size();
  int n_crit   = m_options.tree_ncrit;
  int order    = m_options.tree_expansion_order;
  double theta = m_options.theta;
  std::vector<double> charges(n_sites, 0.0);
  std::vector<double> dipoles(3 * n_sites, 0.0);
  std::vector<double> quadrupoles(6 * n_sites, 0.0);

  auto scheme       = m_options.summation_induced_fields;
  bool damp_enabled = m_options.damp_multipole;
  double damping    = 0.0;
  if (damp_enabled) {
    damping = m_options.damping_factor_multipole;
  }

  int max_order = 0;
  for (int i = 0; i < n_sites; ++i) {
    int max_multipole_order = m_potentials[i].max_multipole_order();
    if (max_multipole_order > max_order) {
      max_order = max_multipole_order;
    }
    if (max_multipole_order > 2) {
      throw std::runtime_error(
            "MultipoleFields::compute_tree only supported up to 2nd order.");
    }
    auto ms = m_potentials[i].get_multipoles();
    if (ms.size() == 0) continue;
    auto qs    = ms[0].get_values_vec();
    charges[i] = qs[0];
    if (max_multipole_order > 0) {
      auto mu            = ms[1].get_values_vec();
      dipoles[i * 3 + 0] = mu[0];
      dipoles[i * 3 + 1] = mu[1];
      dipoles[i * 3 + 2] = mu[2];
    }
    if (max_multipole_order > 1) {
      auto theta             = ms[2].get_values_vec();
      quadrupoles[i * 6 + 0] = 0.5 * theta[0];
      quadrupoles[i * 6 + 1] = theta[1];
      quadrupoles[i * 6 + 2] = theta[2];
      quadrupoles[i * 6 + 3] = 0.5 * theta[3];
      quadrupoles[i * 6 + 4] = theta[4];
      quadrupoles[i * 6 + 5] = 0.5 * theta[5];
    }
  }
  Eigen::VectorXd mult_fields = Eigen::VectorXd::Zero(3 * n_sites);
  // TODO: what if no charges are there?!
  std::shared_ptr<Tree<0, 3>> tree = build_shared_tree<0, 3>(
        m_potentials, charges.data(), n_crit, order, theta, damping);

  std::vector<double> fields0_v(3 * n_sites);
  if (scheme == "fmm") {
    tree->compute_field_fmm(fields0_v.data());
  } else {
    tree->compute_field_exact(fields0_v.data());
  }
  mult_fields += Eigen::Map<Eigen::VectorXd>(fields0_v.data(), fields0_v.size());

  if (max_order > 0) {
    std::shared_ptr<Tree<1, 3>> tree1 = build_shared_tree<1, 3>(
          m_potentials, dipoles.data(), n_crit, order, theta, damping);
    std::vector<double> fields1_v(3 * n_sites);
    if (scheme == "fmm") {
      tree1->compute_field_fmm(fields1_v.data());
    } else {
      tree1->compute_field_exact(fields1_v.data());
    }
    mult_fields -= Eigen::Map<Eigen::VectorXd>(fields1_v.data(), fields1_v.size());
  }
  if (max_order > 1) {
    std::shared_ptr<Tree<2, 3>> tree2 = build_shared_tree<2, 3>(
          m_potentials, quadrupoles.data(), n_crit, order, theta, damping);
    ;
    std::vector<double> fields2_v(3 * n_sites);
    if (scheme == "fmm") {
      tree2->compute_field_fmm(fields2_v.data());
    } else {
      tree2->compute_field_exact(fields2_v.data());
    }
    mult_fields += Eigen::Map<Eigen::VectorXd>(fields2_v.data(), fields2_v.size());
  }
  return mult_fields;
}

Eigen::VectorXd InducedMoments::compute(const Eigen::VectorXd& rhs, Eigen::VectorXd guess,
                                        bool make_guess) {
  // TODO: cleanup, print warning -> outside of solver
  BMatrix bmat(m_polsites, m_options);
  bool converged = false;

  Eigen::VectorXd x0;
  if (make_guess) {
    x0 = bmat.apply_diagonal_inverse(rhs);
  } else {
    x0 = guess;
  }
  Eigen::VectorXd r0 = rhs - bmat.apply(x0);
  Eigen::VectorXd z0 = bmat.apply_diagonal_inverse(r0);
  Eigen::VectorXd p  = z0;

  Eigen::VectorXd x_k1, r_k1, z_k1;
  double alpha_k, beta_k;

  std::vector<Eigen::VectorXd> x{x0};
  std::vector<Eigen::VectorXd> r{r0};
  std::vector<Eigen::VectorXd> z{z0};
  for (int k = 0; k < m_options.maxiter; ++k) {
    Eigen::VectorXd Ap = bmat.apply(p);
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
    z_k1 = bmat.apply_diagonal_inverse(r_k1);
    z.push_back(z_k1);
    beta_k = z_k1.dot(r_k1) / z[k].dot(r[k]);
    p      = z_k1 + beta_k * p;
  }
  if (!converged) {
    throw std::runtime_error("Failed to converge induced dipole moments.");
  }
  return x.back();
}

// returns a vector of potentials that have polarizabilities
std::vector<Potential> get_polarizable_sites(std::vector<Potential> potentials) {
  std::vector<Potential> result;
  std::copy_if(potentials.begin(), potentials.end(), std::back_inserter(result),
               [](Potential& p) { return p.is_polarizable(); });
  return result;
}

}  // namespace libcppe
