cmake_minimum_required(VERSION 3.6)
project(pelib_interface Fortran)

# include_directories(
#   ${PROJECT_BINARY_DIR}/../../external/gen1int/modules
#   ${PROJECT_BINARY_DIR}/../../modules
#   ${PROJECT_BINARY_DIR}/../../external/pelib/modules
# )

include_directories(${CMAKE_INSTALL_PREFIX}/modules)

add_library(pelib_interface STATIC
  pe_interface.F90
)
set_target_properties(pelib_interface PROPERTIES Fortran_FORMAT "FREE")
