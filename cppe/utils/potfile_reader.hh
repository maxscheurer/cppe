#include <iostream>
#include <string>
#include <vector>

class PotfileReader {
private:
  std::string m_potfile;

public:
  PotfileReader (std::string potfile_name) : m_potfile(potfile_name) {};
  ~PotfileReader () {};
};