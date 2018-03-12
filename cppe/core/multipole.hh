#include <vector>

namespace libcppe {
  
class Multipole {
private:
  std::vector<double> m_values;

public:
  Multipole (unsigned k, double x, double y, double z) :
    m_k(k), m_x(x), m_y(y), m_z(z)
    {
      
    };
  ~Multipole () {};
  
  void add_value(double val) {
    m_values.push_back(val);
  }
  
  std::vector<double>& get_values() {return m_values;}
  unsigned m_k;
  double m_x, m_y, m_z;
  
};


}