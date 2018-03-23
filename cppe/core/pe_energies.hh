#ifndef LIBCPPE_CORE_PE_ENERGIES_H
#define LIBCPPE_CORE_PE_ENERGIES_H

#include <map>
#include <string>
#include "../utils/string_utils.hh"

namespace libcppe {

struct PeEnergyContribution {
  PeEnergyContribution(std::string cat, std::string name, double val) :
  m_category(cat), m_name(name), m_value(val) { }
  

	std::string m_category;
  std::string m_name;
  double m_value;
};

struct PeEnergy {
  
  PeEnergy();
  
  double get(std::string energy_string);
  
  void set(std::string energy_string, double energy);
  
  double get_total_energy();
  
  private:
    std::vector<PeEnergyContribution> m_energies;
  
};

} /* libcppe */


#endif // LIBCPPE_CORE_PE_ENERGIES_H