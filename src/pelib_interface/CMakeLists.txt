cmake_minimum_required(VERSION 3.0.3)
project(pelib_interface Fortran)

include_directories(${CMAKE_INSTALL_PREFIX}/modules)

add_library(pelib_interface STATIC
  pe_interface.F90
)
set_target_properties(pelib_interface PROPERTIES Fortran_FORMAT "FREE")

target_link_libraries(pelib_interface "gen1int" "gen1int_interface" "lapack" "blas" "gfortran")

install(TARGETS pelib_interface EXPORT pelib_interface_install LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
