find_package(pelib)
if(${pelib_FOUND})
    add_library(pelib_external INTERFACE)
    message(STATUS "pelib found here: " ${pelib_DIR})
else()
    include(ExternalProject)
    message(STATUS "Suitable pelib could not be located, ${Magenta}Building pelib${ColourReset} instead.")
    ExternalProject_Add(
      pelib_external
      DEPENDS gen1int_external gen1int_interface
      GIT_REPOSITORY "https://gitlab.com/pe-software/pelib-public.git"
      GIT_TAG "7b045d7"
      CMAKE_ARGS -DINTEGRAL_LIBRARY=GEN1INT
                 -DHOST_PROGRAM=DALTON
                 -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                 -DPARENT_MODULE_DIR=${STAGED_INSTALL_PREFIX}/modules
    )
    set(pelib_DIR ${STAGED_INSTALL_PREFIX}/lib CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
    message(STATUS "pelib_DIR = " ${pelib_DIR})
endif()