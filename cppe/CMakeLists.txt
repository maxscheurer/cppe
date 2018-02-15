set(CMAKE_CXX_STANDARD 11)

add_library(cppe
  libcppe.cc
)

target_link_libraries(cppe pe_interface armadillo lapack arpack blas
dalton gen1int gen1int_interface pelib gfortran)
