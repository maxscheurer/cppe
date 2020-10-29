//############################################
//#
//# Functions for running the Barnes-Hut
//# method for gravitational source particles.
//#
//# (C) Ryan Pepper, 2018
//# University of Southampton, UK
//#
//#
//###########################################
#pragma once
#include "tree.hh"
#include "utils.hh"
#include <iostream>

void M_sanity_check(const std::vector<Cell>& cells);

template <int m_order, int osize>
void evaluate_P2M(std::vector<Particle>& particles, std::vector<Cell>& cells, size_t cell,
                  size_t ncrit, size_t exporder);

template <int m_order, int osize>
void evaluate_M2M(std::vector<Particle>& particles, std::vector<Cell>& cells,
                  size_t exporder);

template <int m_order, int osize>
void evaluate_L2L(std::vector<Cell>& cells, size_t exporder);

template <int m_order, int osize>
void evaluate_L2P(std::vector<Particle>& particles, std::vector<Cell>& cells, double* F,
                  size_t ncrit, size_t exporder);

template <int m_order, int osize>
void evaluate_direct(std::vector<Particle>& particles, double* F, size_t Nparticles);

template <int m_order, int osize>
void interact_dehnen(size_t A, size_t B, std::vector<Cell>& cells,
                     std::vector<Particle>& particles, double theta, size_t order,
                     size_t ncrit, double* F);

template <int m_order, int osize>
void interact_dehnen_lazy(const size_t A, const size_t B, const std::vector<Cell>& cells,
                          const std::vector<Particle>& particles, const double theta,
                          const size_t order, const size_t ncrit,
                          std::vector<std::pair<size_t, size_t>>& M2L_list,
                          std::vector<std::pair<size_t, size_t>>& P2P_list);

template <int m_order, int osize>
void P2P_Cells(size_t A, size_t B, std::vector<Cell>& cells,
               std::vector<Particle>& particles, double* F);

template <int m_order, int osize>
void evaluate_M2L_lazy(std::vector<Cell>& cells,
                       std::vector<std::pair<size_t, size_t>>& M2L_list, size_t order);

template <int m_order, int osize>
void evaluate_P2P_lazy(std::vector<Cell>& cells, std::vector<Particle>& particles,
                       std::vector<std::pair<size_t, size_t>>& P2P_list, double* F);

template <int m_order, int osize>
void evaluate_M2P_and_P2P(std::vector<Particle>& particles, unsigned int p,
                          unsigned int i, std::vector<Cell>& cells, double* F,
                          unsigned int n_crit, double theta, unsigned int exporder);
