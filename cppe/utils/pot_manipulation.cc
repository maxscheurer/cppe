#include <set>
#include <cassert>

#include "pot_manipulation.hh"
#include "../core/math.hh"

namespace libcppe {

bool sortbysec(const std::pair<int,int> &a,
                const std::pair<int,int> &b) {
      return (a.second < b.second);
}

std::vector<Potential> PotManipulator::manipulate(PeOptions& pe_options) {
  std::vector<Potential> new_potentials;
  std::set<int> changed_sites;
  if (!pe_options.pe_border) return m_potentials; // do nothing

  for (Atom a : m_mol) {
    arma::vec a_pos = a.get_pos();
    for (Potential pot : m_potentials) {
      if (arma::norm(a_pos-pot.get_site_position()) <= pe_options.border_options.rmin) {
        if (changed_sites.find(pot.index) == changed_sites.end()) {
          changed_sites.insert(pot.index);
        }
      }
    }
  }

  if (changed_sites.size() == 0) return m_potentials; // do nothing

  if (pe_options.border_options.border_type == BorderType::rem) {
    for (Potential pot : m_potentials) {
      if (changed_sites.find(pot.index) == changed_sites.end()) {
        new_potentials.push_back(pot);
      } else {
        std::cout << "Removing all parameters on site: " << pot.index << std::endl;
      }
    }
  } else if (pe_options.border_options.border_type == BorderType::redist) {
    std::cout << "redistributing moments" << std::endl;
    int nredist = pe_options.border_options.nredist;
    assert(nredist <= m_potentials.size());
    // loop over all sites that are in rmin of core
    std::cout << "changed sites: " << changed_sites.size() << std::endl;
    for (auto site : changed_sites) {
      int redist_order = pe_options.border_options.redist_order;
      std::vector<std::pair<int,double>> neighbor_list;
      std::cout << "site: " << site << std::endl;
      // first element is pot.index, second the distance to site
      for (Potential pot: m_potentials) {
        if (pot.index == site) continue;
        arma::vec dist_vec = m_potentials[site].get_site_position() - pot.get_site_position();
        double dist = arma::norm(dist_vec);
        neighbor_list.push_back(std::pair<int,double>(pot.index, dist));
      }
      sort(neighbor_list.begin(), neighbor_list.end(), sortbysec);
      std::cout << "sorted nbl for site: " << site << std::endl;
      for (int k = 0; k < nredist; ++k) {
        Potential& pot = m_potentials[neighbor_list[k].first];
        if (pot.index == site) continue;
        std::cout << "Redistributing to site " << pot.index << std::endl;
        // must have the same order of multipoles if we want to redist
        // we can warn the user later, but for now, we will stop the program
        assert(pot.get_multipoles().size() == m_potentials[site].get_multipoles().size());
        int m_idx = 0;
        for (Multipole& m : pot.get_multipoles()) {
          if (m.m_k >= redist_order) {
            m_idx++;
            continue;
          } else {
            std::cout << "Before: " << m.get_values_vec() << std::endl;
            // holy fuck... this is never gonna work...
            for (size_t i = 0; i < multipole_components(m.m_k); i++) {
               m.get_values()[i] += m_potentials[site].get_multipoles()[m_idx].get_values()[i] / nredist;
            }
            std::cout << "After: " << m.get_values_vec() << std::endl;
            m_idx++;
          }
        }
      }
    }

    for (Potential pot : m_potentials) {
      if (changed_sites.find(pot.index) == changed_sites.end()) {
        new_potentials.push_back(pot);
      } else {
        std::cout << "Removing all parameters on site: " << pot.index << std::endl;
      }
    }
  }
  assert(new_potentials.size() <= m_potentials.size());
  assert(new_potentials.size() != 0);
  return new_potentials;
}

} /* libcppe */
