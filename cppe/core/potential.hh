#include <vector>

#include "multipole.hh"


namespace libcppe {

class Potential {
private:
  /* data */

public:
  Potential();
  ~Potential() {};
  
  std::vector<Multipole> m_multipoles;
  
};

}