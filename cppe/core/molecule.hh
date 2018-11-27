#ifndef LIBCPPE_CPPE_CORE_MOLECULE_H
#define LIBCPPE_CPPE_CORE_MOLECULE_H

#include <armadillo>

namespace libcppe {

struct Atom {
  int atomic_number;
  int charge;
  double m_x, m_y, m_z;
  Atom(int an) : atomic_number(an) { charge = an; }
  Atom(int an, double x, double y, double z)
      : atomic_number(an), m_x(x), m_y(y), m_z(z) {
    charge = an;
  }

  arma::vec get_pos() {
    arma::vec pos(3);
    pos[0] = m_x;
    pos[1] = m_y;
    pos[2] = m_z;
    return pos;
  }
};

// Molecule is a slightly decorated std::vector, more features to come
struct Molecule : std::vector<Atom> {

  // TODO: ugly, probably highly stupid code?
  arma::vec get_atom_position(int atom) {
    if (this->size() <= atom) {
      throw std::out_of_range("Not enough atoms in Molecule.");
    }
    return (*this)[atom].get_pos();
  }

  ~Molecule() = default;
  Molecule &operator=(const Molecule &) = default;
};

} // namespace libcppe
#endif // LIBCPPE_CPPE_CORE_MOLECULE_H
