#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../utils/string_utils.hh"
#include "math.hh"
#include "potfile_reader.hh"

#define ang2bohr 1.8897261246

namespace libcppe {

struct Site {
  double x, y, z;
};

namespace {
inline bool file_exists(const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}
}  // unnamed namespace

PotfileReader::PotfileReader(std::string potfile_name) : m_potfile(potfile_name) {
  if (!file_exists(m_potfile)) {
    throw std::runtime_error("Potential file does not exist.");
  }
}

std::vector<Potential> PotfileReader::read() {
  std::ifstream infile(m_potfile);
  std::vector<Potential> potentials;

  std::string line;

  int num_sites = 0;
  std::string unit;
  std::vector<Site> sites;
  while (std::getline(infile, line)) {
    if (infile.eof() || infile.fail() || infile.bad()) {
      break;
    }
    if (line.find("@COORDINATES") != std::string::npos) {
      getline(infile, line);
      num_sites = stoi(split(line, ' ')[0]);
      // std::cout << "Number of sites: " << num_sites << std::endl;
      getline(infile, line);
      unit = reduce(line);

      double conversion;
      if (!unit.compare("AA")) {
        conversion = ang2bohr;
      } else if (!unit.compare("AU")) {
        conversion = 1.0;
      } else {
        throw std::runtime_error("Invalid unit for potential file.");
      }

      for (size_t i = 0; i < num_sites; i++) {
        Site site;
        getline(infile, line);
        std::vector<std::string> temp = split(reduce(line), ' ');
        std::string element           = temp[0];

        assert(temp.size() >= 4);
        site.x = stod(temp[1]) * conversion;
        site.y = stod(temp[2]) * conversion;
        site.z = stod(temp[3]) * conversion;
        sites.push_back(site);
        // create an empty potential for the site
        Potential p(site.x, site.y, site.z, i);
        potentials.push_back(p);
      }
    }
    if (line.find("ORDER") != std::string::npos) {
      std::vector<std::string> temp = split(reduce(line), ' ');
      // multipoles
      if (temp.size() == 2) {
        int order = stoi(temp[1]);
        getline(infile, line);
        int num_multipoles = stoi(line);
        int site_before    = -1;
        for (size_t n_mul = 0; n_mul < num_multipoles; n_mul++) {
          getline(infile, line);
          temp         = split(reduce(line), ' ');
          int site_num = stoi(temp[0]) - 1;

          // fill up the array if values were not defined for all sites
          if (site_num != site_before + 1) {
            int diff = site_num - site_before;
            for (size_t d = 1; d < diff; d++) {
              Site site = sites[site_before + d];
              Multipole mul(order);
              for (size_t vl = 1; vl <= multipole_components(order); vl++) {
                mul.add_value(0.0);
              }
              potentials[site_before + d].add_multipole(mul);
            }
          }

          Site site = sites[site_num];
          Multipole mul(order);
          for (size_t vl = 1; vl <= multipole_components(order); vl++) {
            mul.add_value(stod(temp[vl]));
          }
          mul.remove_trace();
          potentials[site_num].add_multipole(mul);
          site_before = site_num;

          // check if multipoles at the end of the list are missing
          if ((n_mul == num_multipoles - 1) && site_num != (num_sites - 1)) {
            int diff = num_sites - site_num;
            for (size_t d = 1; d < diff; d++) {
              Site site = sites[site_num + d];
              Multipole mul(order);
              for (size_t vl = 1; vl <= multipole_components(order); vl++) {
                mul.add_value(0.0);
              }
              potentials[site_num + d].add_multipole(mul);
            }
          }
        }
      } else if (temp.size() == 3) {  // polarizabilities
        int order1 = stoi(temp[1]);
        int order2 = stoi(temp[2]);
        if (order1 != 1 || order2 != 1) {
          throw std::runtime_error(
                "Only dipole-dipole polarizabilities "
                "are currently supported.");
        }
        getline(infile, line);
        int num_polarizabilities = stoi(line);
        for (size_t n_pol = 0; n_pol < num_polarizabilities; n_pol++) {
          getline(infile, line);
          temp         = split(reduce(line), ' ');
          int site_num = stoi(temp[0]) - 1;
          Site site    = sites[site_num];
          // std::cout << site.x << " " << site.y << " " << site.z << " " <<
          // site_num + 1 << std::endl;
          std::vector<double> pol_tmp;
          for (size_t vl = 1; vl <= multipole_components(order1 + order2); vl++) {
            pol_tmp.push_back(stod(temp[vl]));
          }
          Polarizability pol{pol_tmp};
          if (potentials[site_num].is_polarizable()) {
            throw std::runtime_error("Potential can only have one polarizability.");
          }
          potentials[site_num].set_polarizability(pol);
        }
      } else {  // unhandled
        throw std::runtime_error("Invalid number in potfile ORDER.");
      }
    }
    if (line.find("EXCLISTS") != std::string::npos) {
      getline(infile, line);
      int num_excl = stoi(split(line, ' ')[0]);
      std::vector<std::string> temp;
      for (size_t i = 0; i < num_excl; i++) {
        getline(infile, line);
        temp         = split(reduce(line), ' ');
        int site_num = stoi(temp[0]) - 1;
        int excl_site_number;
        int counter = 0;
        for (auto s : temp) {
          if (counter == 0) {
            counter++;
            continue;
          }
          int excl_site_number = stoi(s) - 1;
          if (excl_site_number < 0) {
            counter++;
            continue;
          }
          potentials[site_num].add_exclusion(excl_site_number);
          counter++;
        }
      }
    }
  }
  infile.close();

  // DEBUG
  // int sc = 0;
  // for (auto& pot : potentials) {
  //   std::cout << sc << std::endl;
  //   for (auto ex : pot.get_exclusions()) {
  //     std::cout << ex << " ";
  //   }
  //   std::cout << "Potential at site " << sc << std::endl;
  //   for (auto& mul : pot) {
  //     std::cout << "    k = " << mul.m_k << ", x= " << mul.m_x << ", y= " <<
  //     mul.m_y << ", z= " << mul.m_z << std::endl; std::cout << "       "; for
  //     (auto& val : mul.get_values()) {
  //       std::cout << val << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  //   sc++;
  // }
  //

  return potentials;
}

}  // namespace libcppe
