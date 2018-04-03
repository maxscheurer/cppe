#include <set>
#include <cassert>

#include "pot_manipulation.hh"

namespace libcppe {

std::vector<Potential> PotManipulator::manipulate(PeOptions& pe_options) {
  std::vector<Potential> new_potentials;
  std::set<int> removed_sites;
  if (!pe_options.pe_border) return m_potentials; // do nothing

  if (pe_options.border_options.border_type == BorderType::rem) {
    for (Atom a : m_mol) {
      arma::vec a_pos = a.get_pos();
      for (Potential pot : m_potentials) {
        if (arma::norm(a_pos-pot.get_site_position()) <= pe_options.border_options.rmin) {
          if (removed_sites.find(pot.index) == removed_sites.end()) {
            removed_sites.insert(pot.index);
            std::cout << "Removing all parameters on site: " << pot.index << std::endl;
          }
        }
      }
    }
    for (Potential pot : m_potentials) {
      if (removed_sites.find(pot.index) == removed_sites.end()) {
        new_potentials.push_back(pot);
      }
    }
  } else if (pe_options.border_options.border_type == BorderType::redist) {
    throw std::runtime_error("Redistribution of border not implemented yet.");
  }

  assert(new_potentials.size() <= m_potentials.size());

  return new_potentials;
}

} /* libcppe */
