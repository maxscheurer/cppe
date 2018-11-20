set(CMAKE_CXX_STANDARD 14)

# add_subdirectory(pe_interface)
# add_dependencies(pe_interface gen1int_interface pelib)
find_package(Armadillo REQUIRED)

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
target_link_libraries(cppe ${ARMADILLO_LIBRARIES})
target_include_directories(cppe PUBLIC ${ARMADILLO_INCLUDE_DIRS})

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(PN "cppe")

set(INSTALL_DEVEL_HEADERS ON)
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
