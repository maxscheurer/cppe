# define object library for sources common to C library and Python bindings
add_library(cppe-objs
  OBJECT
    metadata.cc
  )

# create empty list of CPPE private header files
add_subdirectory(utils)
add_subdirectory(core)

# OpenMP
if(ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
  if(NOT TARGET OpenMP::OpenMP_CXX)
    message(WARNING "Suitable OpenMP could not be found. Building without OpenMP instead.")
  endif()
endif()

# report on compiler flags in use
message(STATUS "Configuring a ${CMAKE_BUILD_TYPE} build")
if(BUILD_SHARED_LIBS)
  message(STATUS "Building a SHARED library.")
else()
  message(STATUS "Building a STATIC library.")
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} _cmake_build_type_upper)

message(STATUS "C++ compiler flags")
message(STATUS "   From environment       : ${CMAKE_CXX_FLAGS}")
message(STATUS "   Build-type-specific    : ${CMAKE_CXX_FLAGS_${_cmake_build_type_upper}}")
message(STATUS "   Vectorization flag     : ${ARCH_FLAG}")
message(STATUS "   Project defaults       : ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} ${CPPE_CXX_FLAGS}")
message(STATUS "   User-appended          : ${EXTRA_CXXFLAGS}")
if(TARGET OpenMP::OpenMP_CXX AND ENABLE_OPENMP)
  message(STATUS "   OpenMP parallelization : ${OpenMP_CXX_FLAGS}")
endif()

# transform CPPE_CXX_FLAGS and EXTRA_CXXFLAGS to ;-separated lists
string(REPLACE " " ";" CPPE_CXX_FLAGS ${CPPE_CXX_FLAGS})
if(DEFINED EXTRA_CXXFLAGS)
  string(REPLACE " " ";" EXTRA_CXXFLAGS ${EXTRA_CXXFLAGS})
endif()

# Eigen3
find_package(Eigen3 3.4 CONFIG REQUIRED)

# flags to prepend
target_compile_options(cppe-objs
  BEFORE
  PRIVATE
    "${CPPE_CXX_FLAGS}"
  PUBLIC
    "${ARCH_FLAG}"
  )

# flags to append
target_compile_options(cppe-objs
  PRIVATE
    "${EXTRA_CXXFLAGS}"
  )

if(WIN32)
  target_compile_definitions(cppe-objs
    PUBLIC
      _USE_MATH_DEFINES
      _CRT_SECURE_NO_WARNINGS
      _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
    )
endif()

if(APPLE)
  target_link_options(cppe-objs
    PUBLIC
      "SHELL:-undefined dynamic_lookup"
    )
endif()

target_link_libraries(cppe-objs
  PUBLIC
    Eigen3::Eigen
    $<$<AND:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,$<BOOL:${ENABLE_OPENMP}>>:OpenMP::OpenMP_CXX>
  )

# standard locations of libraries, headers, etc.
include(GNUInstallDirs)

# build C++ library
# we left STATIC/SHARED unspecified such that the setting of BUILD_SHARED_LIBS is honored
add_library(cppe "")

target_compile_definitions(cppe
  INTERFACE
    USING_${PROJECT_NAME}
  )

target_include_directories(cppe
  INTERFACE
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

target_link_libraries(cppe
  PRIVATE
    cppe-objs
  )

set_target_properties(cppe
  PROPERTIES
    VERSION ${cppe_VERSION_MAJOR}
    SOVERSION ${cppe_VERSION_MAJOR}
    MACOSX_RPATH ON
    OUTPUT_NAME "cppe"
    EXPORT_NAME "cppe"
    PUBLIC_HEADER "libcppe.hh"
    PRIVATE_HEADER "metadata.hh"
  )

if(ENABLE_PYTHON_INTERFACE)
  add_subdirectory(python_iface)
endif()

install(
  TARGETS
    cppe
  EXPORT
    "${PROJECT_NAME}Targets"
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  PRIVATE_HEADER DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  )

include(CMakePackageConfigHelpers)

set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION
    ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY AnyNewerVersion
  )

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )

install(
  EXPORT
    "${PROJECT_NAME}Targets"
  NAMESPACE
    "${PROJECT_NAME}::"
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )
