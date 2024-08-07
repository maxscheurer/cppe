#include "calculate.hh"
#include "../math.hh"
#include "field_m1.hh"

#include "tree.hh"
#include "utils.hh"
#include <cmath>
#include <iostream>
#include <stack>

#include "../tensors.hh"

using namespace libcppe;

void M_sanity_check(const std::vector<Cell>& cells) {
  double M0 = 0;
  for (auto c = 1; c < cells.size(); c++) {
    if (cells[c].nchild == 0) {
      M0 += cells[c].M[0];
    }
  }
  std::cout << "Cell 0 has M0 = " << cells[0].M[0] << std::endl;
  std::cout << "Other cells   = " << M0 << std::endl;
  if (std::abs((cells[0].M[0] - M0) / M0) > 10e-10) {
    throw std::runtime_error("M0 sanity check failed");
  }
}

template <int m_order, int osize>
void P2P_Cells(size_t A, size_t B, std::vector<Cell>& cells,
               std::vector<Particle>& particles, double* F) {
  // A - target
  // B - source
  for (auto p1 = 0; p1 < cells[A].nleaf; p1++) {
    size_t l1 = cells[A].leaf[p1];
    for (auto p2 = 0; p2 < cells[B].nleaf; p2++) {
      size_t l2 = cells[B].leaf[p2];
      if (l2 != l1 && !particles[l1].excludes_particle(l2)) {
        double dx = (particles[l1].r[0] - particles[l2].r[0]);
        double dy = (particles[l1].r[1] - particles[l2].r[1]);
        double dz = (particles[l1].r[2] - particles[l2].r[2]);
        P2P<m_order, osize>(dx, dy, dz, particles[l2].S, &F[osize * l1]);
      }
    }
  }
}

template <int m_order, int osize>
void P2P_Cells_damping(size_t A, size_t B, std::vector<Cell>& cells,
                       std::vector<Particle>& particles, double* F, double damping) {
  // A - target
  // B - source
  for (auto p1 = 0; p1 < cells[A].nleaf; p1++) {
    size_t l1 = cells[A].leaf[p1];
    for (auto p2 = 0; p2 < cells[B].nleaf; p2++) {
      size_t l2 = cells[B].leaf[p2];
      if (l2 != l1 && !particles[l1].excludes_particle(l2)) {
        double dx = (particles[l1].r[0] - particles[l2].r[0]);
        double dy = (particles[l1].r[1] - particles[l2].r[1]);
        double dz = (particles[l1].r[2] - particles[l2].r[2]);
        double a1 = particles[l1].alpha;
        double a2 = particles[l2].alpha;
        if (a1 > 0.0 && a2 > 0.0) {
          double v = damping / std::pow(a1 * a2, 1.0 / 6.0);
          P2P<m_order, osize>(dx, dy, dz, v, particles[l2].S, &F[osize * l1]);
        } else {
          P2P<m_order, osize>(dx, dy, dz, particles[l2].S, &F[osize * l1]);
        }
      }
    }
  }
}

template <int m_order, int osize>
void interact_dehnen_lazy(const size_t A, const size_t B, const std::vector<Cell>& cells,
                          const std::vector<Particle>& particles, const double theta,
                          const size_t order, const size_t ncrit,
                          std::vector<std::pair<size_t, size_t>>& M2L_list,
                          std::vector<std::pair<size_t, size_t>>& P2P_list) {
  const double dx = cells[A].x - cells[B].x;
  const double dy = cells[A].y - cells[B].y;
  const double dz = cells[A].z - cells[B].z;
  const double R  = sqrt(dx * dx + dy * dy + dz * dz);

  if (R * theta > (cells[A].rmax + cells[B].rmax)) {
    // if (cells[A].nleaf < ncrit && cells[B].nleaf < ncrit) {
    std::pair<size_t, size_t> m2l_pair = std::make_pair(B, A);
    M2L_list.push_back(m2l_pair);
    //}
  }

  else if (cells[A].nchild == 0 && cells[B].nchild == 0) {
    if (cells[B].nleaf >= ncrit) {
      std::pair<size_t, size_t> m2l_pair = std::make_pair(B, A);
      M2L_list.push_back(m2l_pair);
      M2L<m_order, osize>(dx, dy, dz, cells[B].M, cells[A].L, order);
    } else {
      // if (cells[A].nleaf < ncrit and cells[B].nleaf < ncrit) {
      std::pair<size_t, size_t> p2p_pair = std::make_pair(A, B);
      P2P_list.push_back(p2p_pair);
      //}
    }
  }

  else if (cells[B].nchild == 0 ||
           (cells[A].rmax >= cells[B].rmax && cells[A].nchild != 0)) {
    for (auto oa = 0; oa < 8; oa++) {
      // For all 8 children of A, if child exists
      if (cells[A].nchild & (1 << oa)) {
        int a = cells[A].child[oa];
        interact_dehnen_lazy<m_order, osize>(a, B, cells, particles, theta, order, ncrit,
                                             M2L_list, P2P_list);
      }
    }
  }

  else {
    for (auto ob = 0; ob < 8; ob++) {
      // for all 8 children of B, if child exists:
      if (cells[B].nchild & (1 << ob)) {
        int b = cells[B].child[ob];
        interact_dehnen_lazy<m_order, osize>(A, b, cells, particles, theta, order, ncrit,
                                             M2L_list, P2P_list);
      }
    }
  }
}

template <int m_order, int osize>
void evaluate_P2M(std::vector<Particle>& particles, std::vector<Cell>& cells, size_t cell,
                  size_t ncrit, size_t exporder) {
  int sourcesize = multipole_components(m_order);
#pragma omp for
  for (auto c = 0; c < cells.size(); c++) {
    // std::cout << "Cell " << c << std::endl;
    // std::cout << "  Msize = " << Msize(exporder, FMMGEN_SOURCEORDER) << std::endl;
    size_t msize = Msize(exporder, m_order);
    double* M    = new double[msize]();

    if (cells[c].nleaf < ncrit) {
      for (auto i = 0; i < cells[c].nleaf; i++) {
        size_t l = cells[c].leaf[i];
        // Walter dehnen's definition:
        // (-1)^m / m! (x_a - z_a)^m
        double dx = (cells[c].x - particles[l].r[0]);
        double dy = (cells[c].y - particles[l].r[1]);
        double dz = (cells[c].z - particles[l].r[2]);
        for (auto k = 0; k < sourcesize; k++) {
          // std::cout << particles[l].S[k] << std::endl;
          M[k] = particles[l].S[k];
        }
        M2M<m_order, osize>(dx, dy, dz, M, cells[c].M, exporder);
      }
    }
    delete[] M;
  }
}

template <int M, int osize>
void evaluate_M2M(std::vector<Particle>& particles, std::vector<Cell>& cells,
                  size_t exporder) {
  /*
  evaluate_M2M(particles, cells)

  This function evaluates the multipole to
  multipole kernel. It does this by working up the
  tree from the leaf nodes, which is possible
  by iterating backwards through the nodes because
  of the way the tree is constructed.
  */
  // #pragma omp for schedule(dynamic)
  // Can't currently go up the tree in parallel.
  // Needs to be recursive or summing not correct.

  // Dehnen definition:
  // M_m(z_p) = (z_p - z_c)^n / n! M_{m - n}
  for (auto c = cells.size() - 1; c > 0; c--) {
    size_t p  = cells[c].parent;
    double dx = (cells[p].x - cells[c].x);
    double dy = (cells[p].y - cells[c].y);
    double dz = (cells[p].z - cells[c].z);
    M2M<M, osize>(dx, dy, dz, cells[c].M, cells[p].M, exporder);
  }
}

template <int m_order, int osize>
void evaluate_M2L_lazy(std::vector<Cell>& cells,
                       std::vector<std::pair<size_t, size_t>>& M2L_list, size_t order) {
#pragma omp for
  for (auto i = 0; i < M2L_list.size(); i++) {
    size_t B = M2L_list[i].first;
    size_t A = M2L_list[i].second;
    // Dehnen definition:
    // F_n(z_B) = M_m(z_A) * D_{n+m} (z_B - z_A)

    // So here, we've reversed the order:
    // F_n(z_A) = M_m(z_B) * D_{n+m} (z_A - z_B)
    double dx = cells[A].x - cells[B].x;
    double dy = cells[A].y - cells[B].y;
    double dz = cells[A].z - cells[B].z;
    M2L<m_order, osize>(dx, dy, dz, cells[B].M, cells[A].L, order);
  }
}

template <int M, int osize>
void evaluate_P2P_lazy(std::vector<Cell>& cells, std::vector<Particle>& particles,
                       std::vector<std::pair<size_t, size_t>>& P2P_list, double* F,
                       double damping) {

  if (damping > 0.0) {
    #pragma omp for
    for (auto i = 0; i < P2P_list.size(); i++) {
      size_t A = P2P_list[i].first;
      size_t B = P2P_list[i].second;
      P2P_Cells_damping<M, osize>(A, B, cells, particles, F, damping);
    }
  } else {
    #pragma omp for
    for (auto i = 0; i < P2P_list.size(); i++) {
      size_t A = P2P_list[i].first;
      size_t B = P2P_list[i].second;
      P2P_Cells<M, osize>(A, B, cells, particles, F);
    }
  }
}

template <int M, int osize>
void evaluate_L2L(std::vector<Cell>& cells, size_t exporder) {
  // Can't currently go down the tree in parallel!
  // needs to be recursive or summing not correct.
  for (auto p = 0; p < cells.size(); p++) {
    for (auto octant = 0; octant < 8; octant++) {
      if (cells[p].nchild & (1 << octant)) {
        // for child c in cell p
        size_t c  = cells[p].child[octant];
        double dx = cells[c].x - cells[p].x;
        double dy = cells[c].y - cells[p].y;
        double dz = cells[c].z - cells[p].z;
        L2L<M, osize>(dx, dy, dz, cells[p].L, cells[c].L, exporder);
      }
    }
  }
}

template <int M, int osize>
void evaluate_L2P(std::vector<Particle>& particles, std::vector<Cell>& cells, double* F,
                  size_t ncrit, size_t exporder) {
#pragma omp for schedule(runtime)
  for (auto i = 0; i < cells.size(); i++) {
    if (cells[i].nleaf < ncrit) {
      for (auto p = 0; p < cells[i].nleaf; p++) {
        size_t k  = cells[i].leaf[p];
        double dx = particles[k].r[0] - cells[i].x;
        double dy = particles[k].r[1] - cells[i].y;
        double dz = particles[k].r[2] - cells[i].z;
        L2P<M, osize>(dx, dy, dz, cells[i].L, &F[osize * k], exporder);
      }
    }
  }
}

template <int M, int osize>
void evaluate_direct(std::vector<Particle>& particles, double* F) {
  int n = particles.size();
#pragma omp parallel for schedule(runtime)
  for (auto i = 0; i < n; i++) {
    for (auto j = 0; j < n; j++) {
      if (i != j && !particles[i].excludes_particle(j)) {
        double dx = particles[i].r[0] - particles[j].r[0];
        double dy = particles[i].r[1] - particles[j].r[1];
        double dz = particles[i].r[2] - particles[j].r[2];
        P2P<M, osize>(dx, dy, dz, particles[j].S, &F[osize * i]);
      }
    }
  }
}

template <int M, int osize>
void evaluate_direct_damping(std::vector<Particle>& particles, double* F,
                             double damping) {
  int n = particles.size();
#pragma omp parallel for schedule(runtime)
  for (auto i = 0; i < n; i++) {
    for (auto j = 0; j < n; j++) {
      if (i != j && !particles[i].excludes_particle(j)) {
        double dx = particles[i].r[0] - particles[j].r[0];
        double dy = particles[i].r[1] - particles[j].r[1];
        double dz = particles[i].r[2] - particles[j].r[2];
        double a1 = particles[i].alpha;
        double a2 = particles[j].alpha;
        if (std::abs(a1) > 0.0 && std::abs(a2) > 0.0) {
          double v = damping / std::pow(a1 * a2, 1.0 / 6.0);
          P2P<M, osize>(dx, dy, dz, v, particles[j].S, &F[osize * i]);
        } else {
          P2P<M, osize>(dx, dy, dz, particles[j].S, &F[osize * i]);
        }
      }
    }
  }
}

template <int M, int osize>
void evaluate_M2P_and_P2P(std::vector<Particle>& particles, unsigned int p,
                          unsigned int i, std::vector<Cell>& cells, double* F,
                          unsigned int n_crit, double theta, unsigned int exporder) {
  // For particle i, start at cell p
  double dx, dy, dz, r;
  int c;
  unsigned int j;
  // If cell p is not a leaf cell:
  if (cells[p].nleaf >= n_crit) {
    // Iterate through it's children
    for (unsigned int octant = 0; octant < 8; octant++) {
      // If a child exists in a given octant:
      if (cells[p].nchild & (1 << octant)) {
        // Get the child's index
        c = cells[p].child[octant];
        // Calculate the distance from the particle to the child cell
        dx = particles[i].r[0] - cells[c].x;
        dy = particles[i].r[1] - cells[c].y;
        dz = particles[i].r[2] - cells[c].z;
        r  = sqrt(dx * dx + dy * dy + dz * dz);
        // Apply the Barnes-Hut criterion:
        if (cells[c].r > theta * r) {
          // If the cell is 'near':
          evaluate_M2P_and_P2P<M, osize>(particles, c, i, cells, F, n_crit, theta,
                                         exporder);
        } else {
          // If the cell is 'far', calculate the potential and field
          // on the particle from the multipole expansion:
          M2P<M, osize>(dx, dy, dz, cells[c].M, &F[osize * i], exporder);
        }
      }
    }
  } else {
    // loop in leaf cell's particles
    for (unsigned int l = 0; l < (cells[p].nleaf); l++) {
      // Get the particle index:
      j = cells[p].leaf[l];
      if (i != j && !particles[i].excludes_particle(j)) {
        // Calculate the interparticle distance:
        dx = particles[i].r[0] - particles[j].r[0];
        dy = particles[i].r[1] - particles[j].r[1];
        dz = particles[i].r[2] - particles[j].r[2];
        // Compute the field:
        P2P<M, osize>(dx, dy, dz, particles[j].S, &F[osize * i]);
      }
    }
  }
}

#define INSTANTIATE(M_ORDER, OUTPUT_SIZE)                                                \
  template void evaluate_P2M<M_ORDER, OUTPUT_SIZE>(                                      \
        std::vector<Particle> & particles, std::vector<Cell> & cells, size_t cell,       \
        size_t ncrit, size_t exporder);                                                  \
                                                                                         \
  template void evaluate_M2M<M_ORDER, OUTPUT_SIZE>(                                      \
        std::vector<Particle> & particles, std::vector<Cell> & cells, size_t exporder);  \
  template void evaluate_L2L<M_ORDER, OUTPUT_SIZE>(std::vector<Cell> & cells,            \
                                                   size_t exporder);                     \
                                                                                         \
  template void evaluate_L2P<M_ORDER, OUTPUT_SIZE>(std::vector<Particle> & particles,    \
                                                   std::vector<Cell> & cells, double* F, \
                                                   size_t ncrit, size_t exporder);       \
                                                                                         \
  template void evaluate_direct<M_ORDER, OUTPUT_SIZE>(std::vector<Particle> & particles, \
                                                      double* F);                        \
  template void evaluate_direct_damping<M_ORDER, OUTPUT_SIZE>(                           \
        std::vector<Particle> & particles, double* F, double damping);                   \
                                                                                         \
  template void interact_dehnen_lazy<M_ORDER, OUTPUT_SIZE>(                              \
        const size_t A, const size_t B, const std::vector<Cell>& cells,                  \
        const std::vector<Particle>& particles, const double theta, const size_t order,  \
        const size_t ncrit, std::vector<std::pair<size_t, size_t>>& M2L_list,            \
        std::vector<std::pair<size_t, size_t>>& P2P_list);                               \
                                                                                         \
  template void P2P_Cells<M_ORDER, OUTPUT_SIZE>(                                         \
        size_t A, size_t B, std::vector<Cell> & cells,                                   \
        std::vector<Particle> & particles, double* F);                                   \
                                                                                         \
  template void P2P_Cells_damping<M_ORDER, OUTPUT_SIZE>(                                 \
        size_t A, size_t B, std::vector<Cell> & cells,                                   \
        std::vector<Particle> & particles, double* F, double dampign);                   \
                                                                                         \
  template void evaluate_M2L_lazy<M_ORDER, OUTPUT_SIZE>(                                 \
        std::vector<Cell> & cells, std::vector<std::pair<size_t, size_t>> & M2L_list,    \
        size_t order);                                                                   \
                                                                                         \
  template void evaluate_P2P_lazy<M_ORDER, OUTPUT_SIZE>(                                 \
        std::vector<Cell> & cells, std::vector<Particle> & particles,                    \
        std::vector<std::pair<size_t, size_t>> & P2P_list, double* F, double damping);   \
                                                                                         \
  template void evaluate_M2P_and_P2P<M_ORDER, OUTPUT_SIZE>(                              \
        std::vector<Particle> & particles, unsigned int p, unsigned int i,               \
        std::vector<Cell>& cells, double* F, unsigned int n_crit, double theta,          \
        unsigned int exporder);

INSTANTIATE(0, 3)
INSTANTIATE(1, 3)
INSTANTIATE(2, 3)
