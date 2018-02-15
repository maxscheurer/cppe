SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
add_subdirectory(gen1int)
add_subdirectory(gen1int_interface)

set(INTEGRAL_LIBRARY "GEN1INT")
set(HOST_PROGRAM "DALTON")

set(PARENT_INCLUDE_DIR ${CMAKE_Fortran_MODULE_DIRECTORY})
add_subdirectory(pelib)
add_dependencies(pelib gen1int_interface gen1int)
