#ifndef LIBCPPE_CORE_POTFILE_READER_H
#define LIBCPPE_CORE_POTFILE_READER_H

#include <iostream>
#include <string>
#include <vector>

#include "../core/multipole.hh"

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

#endif  // LIBCPPE_CORE_POTFILE_READER_H