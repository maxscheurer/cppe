#include "utils.hh"

size_t TriangleNumbers(size_t n) { return (n * (n + 1)) / 2; }

size_t Nterms(size_t p) {
  return ((p + 1) * (p + 2) * (p + 3)) / 6;
}

size_t Msize(size_t order, size_t source_order) {
  if (source_order == 0) {
    return Nterms(order);
  }
  return Nterms(order) - Nterms(source_order - 1);
}

size_t Lsize(size_t order, size_t source_order) { return Nterms(order - source_order); }