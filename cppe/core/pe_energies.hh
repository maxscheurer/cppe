#ifndef LIBCPPE_CORE_PE_ENERGIES_H
#define LIBCPPE_CORE_PE_ENERGIES_H

#include <map>
#include <string>
#include "../utils/string_utils.hh"

namespace libcppe {
  
// category and name
std::map<std::string, std::vector<std::string>> ens {
  {"Electrostatic", {"Electronic", "Nuclear", "Multipoles"}}, 
  {"Polarization", {"Electronic", "Nuclear", "Multipoles"}}
};

struct PeEnergyContribution {
  PeEnergyContribution(std::string cat, std::string name, double val) :
  m_category(cat), m_name(name), m_value(val) { }
  

	std::string m_category;
  std::string m_name;
  double m_value;
};

struct PeEnergy {
  
  // create PeEnergy with pre-defined energy names
  PeEnergy() {
    for (auto it : ens) {
	    for (auto nm : it.second) {
      	PeEnergyContribution contrib(it.first, nm, 0.0);	
        m_energies.push_back(contrib);
	    }
    }
  }
  
  double get(std::string energy_string) {
    std::vector<std::string> split_name = split(energy_string, '/');
    double energy = 0.0;
    // one category requested
    if (split_name.size() == 1) {
      for (auto& en : m_energies) {
        if (!en.m_category.compare(energy_string)) {
          energy += en.m_value;
        }
      }
      return energy;
    } else if (split_name.size() == 2) { // specific energy requested
      for (auto& en : m_energies) {
        if (!en.m_category.compare(split_name[0]) && !en.m_name.compare(split_name[1])) {
          energy = en.m_value;
        }
      }
      return energy;
    } else {
      throw std::runtime_error("Invalid energy name specified!");
    }
  }
  
  void set(std::string energy_string, double energy) {
    std::vector<std::string> split_name = split(energy_string, '/');
    if (split_name.size() != 2) {
      throw std::runtime_error("Can only set energy with explicit name!");
    }
    for (auto& en : m_energies) {
      if (!en.m_category.compare(split_name[0]) && !en.m_name.compare(split_name[1])) {
        en.m_value = energy;
      }
    }
  }
  
  double get_total_energy() {
    double energy = 0.0;
    for (auto& en : m_energies) {
      energy += en.m_value;
    }
    return energy;
  }
  
  private:
    std::vector<PeEnergyContribution> m_energies;
  
};

} /* libcppe */


#endif // LIBCPPE_CORE_PE_ENERGIES_H