#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "potential.hh"

namespace libcppe {

class PotfileReader {
 private:
  std::string m_potfile;

 public:
  PotfileReader(std::string potfile_name);
  ~PotfileReader(){};
  std::vector<Potential> read();
};

}  // namespace libcppe
