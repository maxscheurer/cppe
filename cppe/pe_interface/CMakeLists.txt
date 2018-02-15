cmake_minimum_required(VERSION 3.6)
project(pe_interface Fortran)

include_directories(
  ${CMAKE_DALTON_DIR}/DALTON/include
  ${CMAKE_DALTON_DIR}/build/modules
)

add_library(pe_interface STATIC
  pe_interface.F90
)
