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
string(TOUPPER ${CMAKE_BUILD_TYPE} _cmake_build_type_upper)

message(STATUS "C++ compiler flags")
message(STATUS "   From environment       : ${CMAKE_CXX_FLAGS}")
message(STATUS "   Build-type-specific    : ${CMAKE_CXX_FLAGS_${_cmake_build_type_upper}}")
message(STATUS "   Vectorization flag     : ${ARCH_FLAG}")
message(STATUS "   Project defaults       : ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} ${CPPE_CXX_FLAGS}")
message(STATUS "   User-appended          : ${EXTRA_CXXFLAGS}")
message(STATUS "   OpenMP parallelization : ${OpenMP_CXX_FLAGS}")

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
  PUBLIC
    "${CPPE_CXX_FLAGS}"
    "${ARCH_FLAG}"
  )

# flags to append
target_compile_options(cppe-objs
  PUBLIC
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

target_link_libraries(cppe-objs
  PUBLIC
    Eigen3::Eigen
    $<$<AND:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,$<BOOL:${ENABLE_OPENMP}>>:OpenMP::OpenMP_CXX>
  )

# standard locations of libraries, headers, etc.
include(GNUInstallDirs)

add_subdirectory(python_iface)
