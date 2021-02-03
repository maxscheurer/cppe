find_package(pybind11 CONFIG)
if(${pybind11_FOUND})
    add_library(pybind11_external INTERFACE)
    message(STATUS "${Cyan}Found pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
    if(${CMAKE_INSIST_FIND_PACKAGE_pybind11})
        message(FATAL_ERROR "Suitable pybind11 could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable pybind11 could not be located, ${Magenta}Building pybind11${ColourReset} instead.")
    ExternalProject_Add(
      pybind11_external
      URL https://github.com/pybind/pybind11/archive/v2.6.1.tar.gz
      CMAKE_ARGS -DPYBIND11_TEST=OFF
                 -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                 -DCMAKE_CXX_STANDARD=${cppe_CXX_STANDARD}
    )
    set(pybind11_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/pybind11 CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
endif()
