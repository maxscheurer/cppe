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

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(PN "cppe")

set(INSTALL_DEVEL_HEADERS ON)
set(src_prefix "core")
file(GLOB_RECURSE header_files ${src_prefix}/*.hh)
message(${CYAN} ${header_files})
install(FILES ${header_files} libcppe.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cppe)
if (INSTALL_DEVEL_HEADERS)
    install(DIRECTORY ${src_prefix}
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN}
            FILES_MATCHING PATTERN "*.h")
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
