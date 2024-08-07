#include "tree.hh"
#include "../math.hh"
#include "calculate.hh"
#include "field_m1.hh"
#include "utils.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

Cell::Cell(double x, double y, double z, double r, size_t parent, size_t order,
           size_t level, size_t ncrit) {
  this->x      = x;
  this->y      = y;
  this->z      = z;
  this->r      = r;
  this->rmax   = sqrt(3 * r * r);
  this->parent = parent;
  this->level  = level;
  this->child.resize(8, 0);
  this->leaf.resize(ncrit, 0);
  this->nleaf  = 0;
  this->nchild = 0;
}

Cell::~Cell() {
#ifdef FMMLIBDEBUG
  std::cout << "Destructor of Cell called" << std::endl;
#endif
}

Cell::Cell(const Cell& other) {
  this->x      = other.x;
  this->y      = other.y;
  this->z      = other.z;
  this->r      = other.r;
  this->rmax   = other.rmax;
  this->parent = other.parent;
  this->level  = other.level;
  this->child  = other.child;
  std::copy(other.leaf.begin(), other.leaf.end(), std::back_inserter(this->leaf));
  std::copy(other.child.begin(), other.child.end(), std::back_inserter(this->child));
  this->nleaf  = other.nleaf;
  this->nchild = other.nchild;
}

Cell::Cell(Cell&& other) {
  this->x      = other.x;
  this->y      = other.y;
  this->z      = other.z;
  this->r      = other.r;
  this->rmax   = other.rmax;
  this->parent = other.parent;
  this->level  = other.level;
  this->child  = other.child;
  this->M      = other.M;
  this->L      = other.L;
  this->leaf   = other.leaf;
  this->nleaf  = other.nleaf;
  this->nchild = other.nchild;
  other.leaf.clear();
  other.child.clear();
}

void printTreeParticles(std::vector<Cell>& cells, size_t cell, size_t depth) {
  for (auto i = 0; i < depth; i++) {
    std::cout << "         ";
  }
  std::cout << cell << " (" << cells[cell].x << "," << cells[cell].y << ","
            << cells[cell].z << ") : (";
  size_t nchild = 0;
  for (auto octant = 0; octant < 8; octant++) {
    if (cells[cell].nchild & (1 << octant)) {
      nchild += 1;
    }
  }

  if (nchild == 0) {
    for (auto i = 0; i < cells[cell].nleaf; i++) {
      std::cout << cells[cell].leaf[i];
      if (i != (cells[cell].nleaf - 1)) {
        std::cout << ",";
      }
    }
  }
  std::cout << ")" << std::endl;
  for (auto octant = 0; octant < 8; octant++) {
    if (cells[cell].nchild & (1 << octant)) {
      printTreeParticles(cells, cells[cell].child[octant], depth + 1);
    }
  }
}

void add_child(std::vector<Cell>& cells, int octant, size_t p, size_t ncrit,
               size_t order) {
  int c = cells.size();
  // Do not change octant to size_t - otherwise the calculation
  // of x, y, z position through bit masking is *not* correct.
  double r      = cells[p].r / 2.0;
  double x      = cells[p].x + r * ((octant & 1) * 2 - 1);
  double y      = cells[p].y + r * ((octant & 2) - 1);
  double z      = cells[p].z + r * ((octant & 4) / 2 - 1);
  size_t parent = p;
  size_t level  = cells[p].level + 1;
  cells.push_back(Cell(x, y, z, r, parent, order, level, ncrit));
  cells[p].child[octant] = c;
  cells[c].nleaf         = 0;
  cells[p].nchild        = (cells[p].nchild | (1 << octant));
}

void split_cell(std::vector<Cell>& cells, std::vector<Particle>& particles, size_t p,
                size_t ncrit, size_t order) {
  size_t l, c;
  // Do not change octant to size_t - otherwise the calculation
  // of x, y, z position in add_child is not correct!
  int octant;
  for (auto i = 0; i < cells[p].leaf.size(); i++) {
    l      = cells[p].leaf[i];
    octant = (particles[l].r[0] > cells[p].x) + ((particles[l].r[1] > cells[p].y) << 1) +
             ((particles[l].r[2] > cells[p].z) << 2);

    if (!((cells[p].nchild) & (1 << octant))) {
      add_child(cells, octant, p, ncrit, order);
    }
    c                             = cells[p].child[octant];
    cells[c].leaf[cells[c].nleaf] = l;
    cells[c].nleaf += 1;
    if (cells[c].nleaf >= ncrit) {
      split_cell(cells, particles, c, ncrit, order);
    }
  }
}

template <int m_order, int osize>
std::shared_ptr<Tree<m_order, osize>> build_shared_tree(
      std::vector<Potential>& potentials, double* S, size_t ncrit, size_t order,
      double theta, double damping) {
  int sourcesize = multipole_components(m_order);
  int nparticles = potentials.size();
  std::vector<Particle> particles(nparticles);
  bool damping_enabled = damping > 0.0;
  for (auto i = 0; i < nparticles; i++) {
    particles[i].r          = potentials[i].ptr_position();
    particles[i].S          = &S[sourcesize * i];
    particles[i].exclusions = potentials[i].get_exclusions();
    if (damping_enabled && potentials[i].is_polarizable()) {
      particles[i].alpha = potentials[i].get_polarizability().get_isotropic_value();
    }
  }

  // Now create cells list
  std::vector<Cell> cells;
  size_t curr;
  int octant;

  // Compute average position
  double xavg = 0;
  double yavg = 0;
  double zavg = 0;
  for (auto i = 0; i < particles.size(); i++) {
    xavg += particles[i].r[0];
    yavg += particles[i].r[1];
    zavg += particles[i].r[2];
  }

  xavg /= particles.size();
  yavg /= particles.size();
  zavg /= particles.size();
#ifdef FMMLIBDEBUG
  std::cout << "Building Tree: Avg pos = (" << xavg << ", " << yavg << ", " << zavg << ")"
            << std::endl;
#endif
  double xmax = 0;
  double ymax = 0;
  double zmax = 0;

  for (auto i = 0; i < particles.size(); i++) {
    double x = std::abs(particles[i].r[0] - xavg);
    double y = std::abs(particles[i].r[1] - yavg);
    double z = std::abs(particles[i].r[2] - zavg);

    if (x > xmax) xmax = x;
    if (y > ymax) ymax = y;
    if (z > zmax) zmax = z;
  }

  double r =
        (xmax > ymax ? (xmax > zmax ? xmax : zmax) : (ymax > zmax ? ymax : zmax)) * 1.001;
  auto root = Cell(xavg, yavg, zavg, r, 0, order, 0, ncrit);

  cells.push_back(root);
  for (auto i = 0; i < particles.size(); i++) {
    curr = 0;
    while (cells[curr].nleaf >= ncrit) {
      cells[curr].nleaf += 1;
      octant = (particles[i].r[0] > cells[curr].x) +
               ((particles[i].r[1] > cells[curr].y) << 1) +
               ((particles[i].r[2] > cells[curr].z) << 2);
      if (!(cells[curr].nchild & (1 << octant))) {
        add_child(cells, octant, curr, ncrit, order);
      }
      curr = cells[curr].child[octant];
    }
    cells[curr].leaf[cells[curr].nleaf] = i;
    cells[curr].nleaf += 1;
    if (cells[curr].nleaf >= ncrit) {
      split_cell(cells, particles, curr, ncrit, order);
    }
  }

  // Now create tree object, and set properties.
  // Choosing a very simple data type here.
  std::shared_ptr<Tree<m_order, osize>> tree = std::make_shared<Tree<m_order, osize>>();
  tree->theta                                = theta;
  tree->ncrit                                = ncrit;
  tree->order                                = order;
  tree->cells                                = cells;
  tree->particles                            = particles;
  tree->damping                              = damping;

  // Create interaction lists, and sort M2L list for cache efficiency.
  interact_dehnen_lazy<m_order, osize>(0, 0, tree->cells, particles, theta, order, ncrit,
                                       tree->M2L_list, tree->P2P_list);
  std::sort(tree->M2L_list.begin(), tree->M2L_list.end(),
            [](std::pair<size_t, size_t>& left, std::pair<size_t, size_t>& right) {
              return left.first < right.first;
            });

  // Create memory into which each cell can point for the multipole arrays.
  tree->M.resize(tree->cells.size() * Msize(order, m_order), 0.0);
  tree->L.resize(tree->cells.size() * Lsize(order, m_order), 0.0);
  for (auto i = 0; i < tree->cells.size(); i++) {
    tree->cells[i].M = &tree->M[i * Msize(order, m_order)];
    tree->cells[i].L = &tree->L[i * Lsize(order, m_order)];
  }
  return tree;
}

template <int m_order, int osize>
std::shared_ptr<Tree<m_order, osize>> build_shared_tree(
      double* pos, double* S, size_t nparticles, size_t ncrit, size_t order, double theta,
      std::vector<std::vector<int>> exclusion_lists) {
  int sourcesize = multipole_components(m_order);
  // Create particles list for convenience
  std::vector<Particle> particles(nparticles);
  for (auto i = 0; i < nparticles; i++) {
    particles[i].r          = &pos[3 * i];
    particles[i].S          = &S[sourcesize * i];
    particles[i].exclusions = exclusion_lists[i];
  }

  // Now create cells list
  std::vector<Cell> cells;
  size_t curr;
  int octant;

  // Compute average position
  double xavg = 0;
  double yavg = 0;
  double zavg = 0;
  for (auto i = 0; i < particles.size(); i++) {
    xavg += particles[i].r[0];
    yavg += particles[i].r[1];
    zavg += particles[i].r[2];
  }

  xavg /= particles.size();
  yavg /= particles.size();
  zavg /= particles.size();
#ifdef FMMLIBDEBUG
  std::cout << "Building Tree: Avg pos = (" << xavg << ", " << yavg << ", " << zavg << ")"
            << std::endl;
#endif
  double xmax = 0;
  double ymax = 0;
  double zmax = 0;

  for (auto i = 0; i < particles.size(); i++) {
    double x = std::abs(particles[i].r[0] - xavg);
    double y = std::abs(particles[i].r[1] - yavg);
    double z = std::abs(particles[i].r[2] - zavg);

    if (x > xmax) xmax = x;
    if (y > ymax) ymax = y;
    if (z > zmax) zmax = z;
  }

  double r =
        (xmax > ymax ? (xmax > zmax ? xmax : zmax) : (ymax > zmax ? ymax : zmax)) * 1.001;
  auto root = Cell(xavg, yavg, zavg, r, 0, order, 0, ncrit);

  cells.push_back(root);
  for (auto i = 0; i < particles.size(); i++) {
    curr = 0;
    while (cells[curr].nleaf >= ncrit) {
      cells[curr].nleaf += 1;
      octant = (particles[i].r[0] > cells[curr].x) +
               ((particles[i].r[1] > cells[curr].y) << 1) +
               ((particles[i].r[2] > cells[curr].z) << 2);
      if (!(cells[curr].nchild & (1 << octant))) {
        add_child(cells, octant, curr, ncrit, order);
      }
      curr = cells[curr].child[octant];
    }
    cells[curr].leaf[cells[curr].nleaf] = i;
    cells[curr].nleaf += 1;
    if (cells[curr].nleaf >= ncrit) {
      split_cell(cells, particles, curr, ncrit, order);
    }
  }

  // Now create tree object, and set properties.
  // Choosing a very simple data type here.
  std::shared_ptr<Tree<m_order, osize>> tree = std::make_shared<Tree<m_order, osize>>();
  tree->theta                                = theta;
  tree->ncrit                                = ncrit;
  tree->order                                = order;
  tree->cells                                = cells;
  tree->particles                            = particles;

  // Create interaction lists, and sort M2L list for cache efficiency.
  interact_dehnen_lazy<m_order, osize>(0, 0, tree->cells, particles, theta, order, ncrit,
                                       tree->M2L_list, tree->P2P_list);
  std::sort(tree->M2L_list.begin(), tree->M2L_list.end(),
            [](std::pair<size_t, size_t>& left, std::pair<size_t, size_t>& right) {
              return left.first < right.first;
            });

  // Create memory into which each cell can point for the multipole arrays.
  tree->M.resize(tree->cells.size() * Msize(order, m_order), 0.0);
  tree->L.resize(tree->cells.size() * Lsize(order, m_order), 0.0);
  for (auto i = 0; i < tree->cells.size(); i++) {
    tree->cells[i].M = &tree->M[i * Msize(order, m_order)];
    tree->cells[i].L = &tree->L[i * Lsize(order, m_order)];
  }
  return tree;
}

template <int m_order, int osize>
void Tree<m_order, osize>::clear_M() {
  std::fill(M.begin(), M.end(), 0);
}

template <int m_order, int osize>
void Tree<m_order, osize>::clear_L() {
  std::fill(L.begin(), L.end(), 0);
}

template <int m_order, int osize>
void Tree<m_order, osize>::set_sources(double* S) {
  int sourcesize = multipole_components(m_order);
  clear_M();
  clear_L();
  for (auto i = 0; i < particles.size(); i++) {
    particles[i].S = &S[sourcesize * i];
  }
}

template <int m_order, int osize>
void Tree<m_order, osize>::compute_field_fmm(double* F) {
  // std::cout << "Computing FMM fields." << std::endl;
  for (auto i = 0; i < osize * particles.size(); i++) {
    F[i] = 0.0;
  }
  clear_M();
  clear_L();
#pragma omp parallel
  { evaluate_P2M<m_order, osize>(particles, cells, 0, ncrit, order); }

  evaluate_M2M<m_order, osize>(particles, cells, order);

#ifdef FMMLIBDEBUG
  M_sanity_check(cells);
#endif

#pragma omp parallel
  {
    evaluate_M2L_lazy<m_order, osize>(cells, M2L_list, order);
    evaluate_P2P_lazy<m_order, osize>(cells, particles, P2P_list, F, damping);
  }

  evaluate_L2L<m_order, osize>(cells, order);
#pragma omp parallel
  { evaluate_L2P<m_order, osize>(particles, cells, F, ncrit, order); }
}

template <int m_order, int osize>
void Tree<m_order, osize>::compute_field_exact(double* F) {
  if (damping > 0.0) {
    evaluate_direct_damping<m_order, osize>(particles, F, damping);
  } else {
    evaluate_direct<m_order, osize>(particles, F);
  }
}

template class Tree<0, 3>;
template class Tree<1, 3>;
template class Tree<2, 3>;
template std::shared_ptr<Tree<0, 3>> build_shared_tree<0, 3>(
      double* pos, double* S, size_t nparticles, size_t ncrit, size_t order, double theta,
      std::vector<std::vector<int>> exclusion_lists);
template std::shared_ptr<Tree<1, 3>> build_shared_tree<1, 3>(
      double* pos, double* S, size_t nparticles, size_t ncrit, size_t order, double theta,
      std::vector<std::vector<int>> exclusion_lists);
template std::shared_ptr<Tree<2, 3>> build_shared_tree<2, 3>(
      double* pos, double* S, size_t nparticles, size_t ncrit, size_t order, double theta,
      std::vector<std::vector<int>> exclusion_lists);

template std::shared_ptr<Tree<0, 3>> build_shared_tree<0, 3>(
      std::vector<Potential>& potentials, double* S, size_t ncrit, size_t order,
      double theta, double damping);
template std::shared_ptr<Tree<1, 3>> build_shared_tree<1, 3>(
      std::vector<Potential>& potentials, double* S, size_t ncrit, size_t order,
      double theta, double damping);
template std::shared_ptr<Tree<2, 3>> build_shared_tree<2, 3>(
      std::vector<Potential>& potentials, double* S, size_t ncrit, size_t order,
      double theta, double damping);