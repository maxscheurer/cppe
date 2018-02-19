cmake_minimum_required(VERSION 3.6)
project(pe_interface Fortran)

include_directories(
  ${PROJECT_BINARY_DIR}/../../external/gen1int/modules
  ${PROJECT_BINARY_DIR}/../../modules
  ${PROJECT_BINARY_DIR}/../../external/pelib/modules
)

add_library(pe_interface STATIC
  pe_interface.F90
)
