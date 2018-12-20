#ifndef LIBCPPE_CORE_PE_ENERGIES_H
#define LIBCPPE_CORE_PE_ENERGIES_H

#include "../utils/string_utils.hh"
#include <map>
#include <string>

namespace libcppe {

/*! PE Energy Contribution */
struct PeEnergyContribution {
  PeEnergyContribution(std::string cat, std::string name, double val)
      : m_category(cat), m_name(name), m_value(val) {}

  std::string m_category; /*!< category of the energy, either "Electrostatic"
                        or "Polarization"*/
  std::string m_name;     /*!< name of the energy,
                              "Electronic", "Nuclear", or "Multipoles" */
  double m_value;         //!< energy value
};

/*! PE Energy Container */
struct PeEnergy {
  PeEnergy();

  /*! returns energy contribution from given string */
  double get(std::string energy_string);

  /*!
     \brief sets the energy titled energy_string to energy
     \param energy_string: name of the energy contribution
     \param energy: value
     \return void
  */
  void set(std::string energy_string, double energy);

  /*! returns the total PE energy */
  double get_total_energy();

private:
  std::vector<PeEnergyContribution> m_energies;
};

} // namespace libcppe

#endif // LIBCPPE_CORE_PE_ENERGIES_H
