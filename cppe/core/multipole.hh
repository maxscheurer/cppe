namespace libcppe {
  
class Multipole {
private:
  /* data */

public:
  Multipole (unsigned k, double x, double y, double z) :
    m_k(k), m_x(x), m_y(y), m_z(z)
    {
    };
  ~Multipole () {};
  
  unsigned m_k;
  double m_x, m_y, m_z;
  
};


}