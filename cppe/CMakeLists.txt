set(CMAKE_CXX_STANDARD 14)

# add_subdirectory(pe_interface)
# add_dependencies(pe_interface gen1int_interface pelib)
find_package(Armadillo REQUIRED)

if (ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(ENABLE_PYTHON_INTERFACE OFF)
set(PYTHON_INTERFACE_SOURCES "")
set(PYTHON_SHARED "")
if (${ENABLE_PYTHON_INTERFACE})
    message(STATUS "Building Python interface for cppe.")
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
        
          set(pybind11_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/pybind11 CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
        )
    endif()

    set(PYTHON_INTERFACE_SOURCES
        pycppe/export_molecule.cc
        pycppe/export_multipole.cc
        pycppe/export_options.cc
        pycppe/export_fields.cc
        pycppe/export_utils.cc
        pycppe/export_state.cc
        pycppe/pycppe.cc)
    set(PYTHON_SHARED SHARED)
endif()

add_library(cppe ${PYTHON_SHARED}
  libcppe.cc
  utils/potfile_reader.cc
  utils/string_utils.cc
  utils/pot_manipulation.cc
  core/cppe_state.cc
  core/multipole_expansion.cc
  core/electric_fields.cc
  core/math.cc
  core/pe_energies.cc
  ${PYTHON_INTERFACE_SOURCES}
)
# target_link_libraries(cppe pe_interface armadillo lapack blas gen1int gen1int_interface pelib gfortran)
target_link_libraries(cppe PUBLIC ${ARMADILLO_LIBRARIES})
target_include_directories(cppe PUBLIC ${ARMADILLO_INCLUDE_DIRS})

if(${ENABLE_PYTHON_INTERFACE})
    target_link_libraries(cppe PUBLIC pybind11::module)
    set_target_properties(cppe PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                             SUFFIX "${PYTHON_MODULE_EXTENSION}")
endif()

set(PN "cppe")
option(INSTALL_DEVEL_HEADERS "Install devel headers" ON)
if (INSTALL_DEVEL_HEADERS)
    install(DIRECTORY "core"
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN}
            FILES_MATCHING PATTERN "*.hh")
    install(DIRECTORY "utils"
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN}
            FILES_MATCHING PATTERN "*.hh")
    install(FILES libcppe.hh
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN})
endif()
install(TARGETS cppe
        EXPORT "${PN}Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

target_compile_definitions(cppe INTERFACE USING_${PN})
target_include_directories(cppe INTERFACE
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> ${ARMADILLO_INCLUDE_DIRS})

set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
message(STATUS ${CMAKECONFIG_INSTALL_DIR})
configure_package_config_file(cmake/${PN}Config.cmake.in
                         "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                         INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                            VERSION ${${PN}_VERSION}
                            COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
         ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
   DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT "${PN}Targets"
   NAMESPACE "${PN}::"
   DESTINATION ${CMAKECONFIG_INSTALL_DIR})
