#include <cassert>
#include <set>

#include <Eigen/Core>

#include "math.hh"
#include "pot_manipulation.hh"

namespace libcppe {

bool sortbysec(const std::pair<int, double>& a, const std::pair<int, double>& b) {
  return (a.second < b.second);
}

std::vector<Potential> PotManipulator::manipulate(const PeOptions& pe_options) {
  if (pe_options.iso_pol) {
    for (Potential& pot : m_potentials) {
      pot.get_polarizability().make_isotropic();
    }
  }
  return manipulate_border(pe_options);
}

std::vector<Potential> PotManipulator::manipulate_border(const PeOptions& pe_options) {
  std::vector<Potential> new_potentials;
  std::set<int> changed_sites;
  if (!pe_options.pe_border) return m_potentials;  // do nothing

  for (Atom a : m_mol) {
    Eigen::Vector3d a_pos = a.get_pos();
    for (Potential pot : m_potentials) {
      if ((a_pos - pot.get_site_position()).norm() <= pe_options.border_rmin) {
        if (changed_sites.find(pot.index) == changed_sites.end()) {
          changed_sites.insert(pot.index);
        }
      }
    }
  }

  if (changed_sites.size() == 0) return m_potentials;  // do nothing

  if (pe_options.border_type == "remove") {
    for (Potential pot : m_potentials) {
      if (changed_sites.find(pot.index) == changed_sites.end()) {
        new_potentials.push_back(pot);
      } else {
        m_printer("     Removing all parameters on site: " + std::to_string(pot.index));
      }
    }
  } else if (pe_options.border_type == "redist") {
    m_printer("     Redistributing moments in " + std::to_string(pe_options.border_rmin) +
              " a.u. proximity.");
    int nredist = pe_options.border_nredist;
    if (nredist == -1) {
      m_printer("     Border sites will be redistributed to all other sites.");
      nredist = m_potentials.size() - changed_sites.size();
    }
    if (nredist > m_potentials.size()) {
      throw std::runtime_error("Cannot redistribute to more sites than available sites.");
    }
    // loop over all sites that are in rmin of core
    // redistribute their multipoles/polarizabilities
    for (auto site : changed_sites) {
      int redist_order = pe_options.border_redist_order;
      std::vector<std::pair<int, double>> neighbor_list;
      m_printer("     Redistributing site: " + std::to_string(site));
      // first element is pot.index, second the distance to site
      for (Potential pot : m_potentials) {
        if (changed_sites.find(pot.index) != changed_sites.end()) continue;
        Eigen::Vector3d dist_vec =
              m_potentials[site].get_site_position() - pot.get_site_position();
        double dist = dist_vec.norm();
        neighbor_list.push_back(std::pair<int, double>(pot.index, dist));
      }
      sort(neighbor_list.begin(), neighbor_list.end(), sortbysec);
      for (int k = 0; k < nredist; ++k) {
        Potential& pot = m_potentials[neighbor_list[k].first];
        if (pot.index == site) continue;
        m_printer("       to neighbor " + std::to_string(pot.index));
        // must have the same order of multipoles if we want to redist
        int m_idx = 0;
        for (Multipole& m : pot.get_multipoles()) {
          if (m.m_k >= redist_order) {
            m_idx++;
            continue;
          } else if (m.m_k != m_potentials[site].get_multipoles()[m_idx].m_k) {
            m_idx++;
            continue;
          } else {
            // std::cout << "Before: " << std::endl;
            // std::cout << m.get_values_vec() << std::endl;
            for (size_t i = 0; i < multipole_components(m.m_k); i++) {
              m.get_values()[i] +=
                    m_potentials[site].get_multipoles()[m_idx].get_values()[i] /
                    static_cast<double>(nredist);
            }
            // std::cout << "After: " << std::endl;
            // std::cout << m.get_values_vec() << std::endl;
            m_idx++;
          }
        }
        if (pe_options.border_redist_pol && pot.is_polarizable() &&
            m_potentials[site].is_polarizable()) {
          Polarizability& p = pot.get_polarizability();
          // TODO: what if a site that has been chosen nearest neighbor
          // has no multipole moments/polarizability of the respective order?!
          Eigen::Matrix3d other_pol =
                m_potentials[site].get_polarizability().get_matrix();
          p.get_matrix() += other_pol / static_cast<double>(nredist);
        }
      }
    }

    double total_charge = 0.0;
    for (Potential pot : m_potentials) {
      if (changed_sites.find(pot.index) == changed_sites.end()) {
        new_potentials.push_back(pot);
        for (auto& m : pot.get_multipoles()) {
          if (m.m_k == 0) {
            total_charge += m.get_values()[0];
          } else
            break;
        }
      } else {
        m_printer("     Removing all parameters on site: " + std::to_string(pot.index));
      }
    }
    m_printer("    Total charge of the classical region: " +
              std::to_string(total_charge));
  }

  if (new_potentials.size() > m_potentials.size()) {
    throw std::runtime_error(
          "Manipulated potentials cannot have more sites "
          "than original ones.");
  }
  return new_potentials;
}

}  // namespace libcppe
