cmake_minimum_required(VERSION 3.0.0)
project(cppe-core VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(ExternalProject)

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

set(PYTHON_INTERFACE_SOURCES "")
set(PYTHON_SHARED "")
if (${ENABLE_PYTHON_INTERFACE})
    message(STATUS "Building Python interface for cppe.")
    find_package(pybind11 CONFIG REQUIRED)
    # message(STATUS "Using pybind11: " ${pybind11_DIR})

    set(PYTHON_INTERFACE_SOURCES
        pycppe/export_molecule.cc
        pycppe/export_multipole.cc
        pycppe/export_options.cc
        pycppe/export_fields.cc
        pycppe/export_utils.cc
        pycppe/export_state.cc
        pycppe/export_math.cc
        pycppe/pycppe.cc)
    set(PYTHON_SHARED SHARED)
endif()

set(PELIB_LIBRARIES "")
set(PELIB_INTERFACE_SOURCES "")
if(${ENABLE_PELIB})
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    add_subdirectory(pelib_interface)
    set(PELIB_LIBRARIES "pelib"
                        "gen1int"
                        "pelib_interface"
                        "gen1int_interface"
                        "lapack"
                        "blas"
                        "gfortran")
    link_directories(${CMAKE_INSTALL_PREFIX}/lib)
    set(PELIB_INTERFACE_SOURCES pelib_interface/pe_interface.F90)
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
  ${PYTHON_INTERFACE_SOURCES} ${PELIB_INTERFACE_SOURCES}
)

# TODO: hack
find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(cppe PUBLIC ${PELIB_LIBRARIES})
target_include_directories(cppe PUBLIC ${EIGEN3_INCLUDE_DIRS})

if(${ENABLE_PYTHON_INTERFACE})
    target_link_libraries(cppe PUBLIC pybind11::module)
    set_target_properties(cppe PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                             SUFFIX "${PYTHON_MODULE_EXTENSION}")
endif()

set(PN "cppe")
#### Installation ####
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
$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
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

if(${ENABLE_PELIB})
    install(EXPORT pelib_interface_install DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()