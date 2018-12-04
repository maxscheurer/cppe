find_package(pybind11 CONFIG)
if(${pybind11_FOUND})
    add_library(pybind11_external INTERFACE)
    message(STATUS "Pybind11 found here: " ${pybind11_DIR})
else()
    include(ExternalProject)
    message(STATUS "Suitable pybind11 could not be located, ${Magenta}Building pybind11${ColourReset} instead.")
    ExternalProject_Add(
      pybind11_external
      GIT_REPOSITORY "https://github.com/pybind/pybind11.git"
      GIT_TAG "master"
      CMAKE_ARGS -DPYBIND11_TEST=OFF
                 -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    )
    set(pybind11_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/pybind11 CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
endif()