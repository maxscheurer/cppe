#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "potfile_reader.hh"
#include "string_utils.hh"
#include <libpe/libcppe/cppe/core/multipole.hh>

#define ang2bohr 1.8897261246

namespace libcppe {
  
struct Site {
  double x, y, z;
};

std::vector<int> mul_vals{1, 3, 6, 10, 15, 21};
  
namespace {
  inline bool file_exists (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
  }
} // unnamed namespace 
  
PotfileReader::PotfileReader(std::string potfile_name) : 
  m_potfile(potfile_name) {
    
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
  while( std::getline(infile, line) ) {
    if (infile.eof() || infile.fail() || infile.bad() ) {
      break;
    }
    if (line.find("@COORDINATES") != std::string::npos ) {
      getline(infile, line);
      num_sites = stoi(split(line, ' ')[0]);
      std::cout << "Number of sites: " << num_sites << std::endl; 
      getline(infile, line);
      unit = reduce(line);
      std::cout << "Unit: " << unit << std::endl;
      
      for (size_t i = 0; i < num_sites; i++) {
        Site site;
        getline(infile, line);
        std::vector<std::string> temp = split(reduce(line), ' ');
        std::string element = temp[0];
        std::cout << "Element: " << element << std::endl;
        
        double conversion;
        if (!unit.compare("AA")) {
          conversion = ang2bohr;
        } else if (!unit.compare("AU"))  {
          conversion = 1.0;
        } else {
          throw std::runtime_error("Invalid unit for potential file.");
        }
        assert(temp.size() == 4);
        std::cout << "-----" << std::endl;
        std::cout << temp[1] << std::endl;
        std::cout << temp[2] << std::endl;
        std::cout << temp[3] << std::endl;
        std::cout << "-----" << std::endl;
        site.x = stod(temp[1]) * conversion;
        site.y = stod(temp[2]) * conversion;
        site.z = stod(temp[3]) * conversion;
        sites.push_back(site);
        std::vector<Multipole> v_mul;
        potentials.push_back(v_mul);
      }
    }
    if (line.find("ORDER") != std::string::npos) {
      std::vector<std::string> temp = split(reduce(line), ' ');
      // multipoles
      if (temp.size() == 2) {
        int order = stoi(temp[1]);
        getline(infile, line);
        int num_multipoles = stoi(line);
        int site_before = -1;
        for (size_t n_mul = 0; n_mul < num_multipoles; n_mul++) {
          getline(infile, line);
          temp = split(reduce(line), ' ');
          int site_num = stoi(temp[0]) - 1;
          
          // fill up the array if values were not defined for all sites
          if (site_num != site_before+1) {
            int diff = site_num - site_before;
            for (size_t d = 1; d < diff; d++) {
              Site site = sites[site_before + d];
              Multipole mul(order, site.x, site.y, site.z);
              for (size_t vl = 1; vl <= mul_vals[order]; vl++) {
                mul.add_value(0.0);
              }
              potentials[site_before + d].push_back(mul);
            }
          }
          // ----------------------------------
          Site site = sites[site_num];
          Multipole mul(order, site.x, site.y, site.z);
          for (size_t vl = 1; vl <= mul_vals[order]; vl++) {
            mul.add_value(stod(temp[vl]));
          }
          potentials[site_num].push_back(mul);
          site_before = site_num;
          
          // check if multipoles at the end of the list are missing
          if ((n_mul == num_multipoles-1) && site_num != (num_sites-1)) {
            int diff = num_sites - site_num;
            for (size_t d = 1; d < diff; d++) {
              Site site = sites[site_num + d];
              Multipole mul(order, site.x, site.y, site.z);
              for (size_t vl = 1; vl <= mul_vals[order]; vl++) {
                mul.add_value(0.0);
              }
              potentials[site_num + d].push_back(mul);
            }
          }
        }
      } else if (temp.size() == 3) { // polarizabilities
        // read polarizabilities when we have the appropriate integrals
      } else { // unhandled
        throw std::runtime_error("Invalid number in potfile ORDER.");
      }
    }
  }
  infile.close();
  
  // DEBUG
  int sc = 0;
  for (auto& pot : potentials) {
    std::cout << "Potential at site " << sc << std::endl;
    for (auto& mul : pot) {
      std::cout << "    k = " << mul.m_k << ", x= " << mul.m_x << ", y= " << mul.m_y << ", z= " << mul.m_z << std::endl;
      std::cout << "       ";
      for (auto& val : mul.get_values()) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }
    sc++;
  }
  // 
  
  return potentials;
}

  
} // namespace libcppe