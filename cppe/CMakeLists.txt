set(CMAKE_CXX_STANDARD 11)

# add_subdirectory(pe_interface)
# add_dependencies(pe_interface gen1int_interface pelib)

add_library(cppe
  libcppe.cc
  utils/potfile_reader.cc
  utils/string_utils.cc
  utils/pot_manipulation.cc
  core/cppe_state.cc
  core/multipole_expansion.cc
  core/electric_fields.cc
  core/math.cc
  core/pe_energies.cc
)

# target_link_libraries(cppe pe_interface armadillo lapack blas gen1int gen1int_interface pelib gfortran)
target_link_libraries(cppe armadillo)
