# OpenMP
if(ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
  if(NOT TARGET OpenMP::OpenMP_CXX)
    message(WARNING "Suitable OpenMP could not be found. Building without OpenMP instead.")
  endif()
endif()

# define object library for sources
add_library(cppe-objs
  OBJECT
    ""
  )

# create empty list of CPPE private header files
add_subdirectory(utils)
add_subdirectory(core)

# report on compiler flags in use
message(STATUS "Configuring a ${CMAKE_BUILD_TYPE} build")
string(TOUPPER ${CMAKE_BUILD_TYPE} _cmake_build_type_upper)

message(STATUS "Compiler flags for ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "   From environment       : ${CMAKE_CXX_FLAGS}")
message(STATUS "   Build-type-specific    : ${CMAKE_CXX_FLAGS_${_cmake_build_type_upper}}")
message(STATUS "   Vectorization flag     : ${ARCH_FLAG}")
message(STATUS "   Project defaults       : ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} ${CPPE_CXX_FLAGS}")
message(STATUS "   User-appended          : ${EXTRA_CXXFLAGS}")
message(STATUS "   OpenMP parallelization : ${OpenMP_CXX_FLAGS}")

# transform CPPE_CXX_FLAGS and EXTRA_CXXFLAGS to ;-separated lists
if(DEFINED CPPE_CXX_FLAGS)
  string(REPLACE " " ";" CPPE_CXX_FLAGS ${CPPE_CXX_FLAGS})
endif()
if(DEFINED EXTRA_CXXFLAGS)
  string(REPLACE " " ";" EXTRA_CXXFLAGS ${EXTRA_CXXFLAGS})
endif()

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

# Eigen3
message(STATUS "Eigen3_ROOT : $ENV{Eigen3_ROOT}")
message(STATUS "Eigen3_DIR : $ENV{Eigen3_DIR}")
find_package(Eigen3 3.4 REQUIRED)
if(TARGET Eigen3::Eigen)
  message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")
endif()

target_link_libraries(cppe-objs
  PUBLIC
    Eigen3::Eigen
    $<$<AND:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,$<BOOL:${ENABLE_OPENMP}>>:OpenMP::OpenMP_CXX>
  )

add_subdirectory(python_iface)
