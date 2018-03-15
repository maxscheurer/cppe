#ifndef LIBCPPE_CORE_POTFILE_READER_H
#define LIBCPPE_CORE_POTFILE_READER_H

#include <iostream>
#include <string>
#include <vector>

namespace libcppe {

class Multipole;
using Potential = std::vector<Multipole>;

class PotfileReader {
private:
  std::string m_potfile;

public:
  PotfileReader (std::string potfile_name);
  ~PotfileReader () {};
  std::vector<Potential> read();
};

} // namespace libcppe

#endif //LIBCPPE_CORE_POTFILE_READER_H