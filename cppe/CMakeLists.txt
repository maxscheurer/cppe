cmake_minimum_required(VERSION 3.0.0)
project(cppe-core VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(ExternalProject)

# TODO: move to external/upstream
find_package(Armadillo)
if(${ARMADILLO_FOUND})
    add_library(armadillo_external INTERFACE) # Dummy
    set(ARMADILLO_INCLUDE_DIRS_IFC ${ARMADILLO_INCLUDE_DIRS})
else()
    message(STATUS "${CYAN}Armadillo not found on system, building Armadillo from source.")
    set(arma_tar "armadillo-9.200.5.tar.xz")
    set(CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
    ExternalProject_Add(armadillo_external
    	DOWNLOAD_COMMAND wget https://sourceforge.net/projects/arma/files/${arma_tar}
    	CONFIGURE_COMMAND tar -xf ../${arma_tar} -C ../armadillo_external  --strip-components=1
    	BUILD_IN_SOURCE 1
    	BUILD_COMMAND cmake ${CMAKE_ARGS} .
    	INSTALL_COMMAND make
	)

    set(Armadillo_DIR ${CMAKE_BINARY_DIR}/cppe/armadillo_external-prefix/src/armadillo_external)
    set(ARMADILLO_LIBRARIES ${Armadillo_DIR}/libarmadillo${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(ARMADILLO_INCLUDE_DIRS ${Armadillo_DIR}/include)
    set(ARMADILLO_INCLUDE_DIRS_IFC "")
endif()

if (ENABLE_OPENMP)
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
    find_package(pybind11 CONFIG REQUIRED NO_DEFAULT_PATH)
    message(STATUS "Using pybind11: " ${pybind11_DIR})

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

if(${ENABLE_PELIB})
    add_subdirectory(pelib_interface)
    # add_dependencies(pelib_interface gen1int_interface pelib_external gen1int_external)
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

add_dependencies(cppe armadillo_external)
# target_link_libraries(cppe pe_interface armadillo lapack blas gen1int gen1int_interface pelib gfortran)
target_link_libraries(cppe PUBLIC ${ARMADILLO_LIBRARIES})
target_include_directories(cppe PRIVATE ${ARMADILLO_INCLUDE_DIRS})

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
$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> ${ARMADILLO_INCLUDE_DIRS_IFC})

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
