add_library(
    gen1int_interface
    gen1int_api.F90
    gen1int_cube.F90
    gen1int_matrix.F90
    gen1int_shell.F90
    quit.F90
    )

set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
    )

if(WITH_MPI)
  set(ENABLE_MPI ON)
  add_definitions(-DVAR_MPI)
  message("Using MPI for gen1int_interface.")
endif()


# add_external(gen1int)
include_directories(
  ${PROJECT_BINARY_DIR}/modules
  ${PROJECT_BINARY_DIR}/external/gen1int/modules
)
target_link_libraries(gen1int_interface gen1int)
add_dependencies(gen1int_interface gen1int)
