cmake_minimum_required(VERSION 3.0.0)
cmake_policy(VERSION 3.0.0)

# Prefer shared to static libs (cmake does the reverse by default)
# option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones" ON)

##########################################################################
# Project setup

project(cppe LANGUAGES CXX)
if(NOT CMAKE_REQUIRED_QUIET)
	message("\n# Configuring ${PROJECT_NAME}\n#")
endif()

# Export compile commands if cmake supports it.
if (CMAKE_VERSION VERSION_GREATER 3.5.0)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)

add_subdirectory(external/upstream)
include(ExternalProject)
ExternalProject_Add(cppe-core
   DEPENDS pybind11_external
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cppe
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
			  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
              -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
              -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
              -DCMAKE_CXX_STANDARD_REQUIRED=ON
              -DENABLE_OPENMP=${ENABLE_OPENMP}
			  -DENABLE_PYTHON_INTERFACE=${ENABLE_PYTHON_INTERFACE}
			  -Dpybind11_DIR=${pybind11_DIR}
              -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
   CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
   BUILD_ALWAYS 1)


install(DIRECTORY ${STAGED_INSTALL_PREFIX}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS
        PATTERN "*.pyc" EXCLUDE)

##########################################################################
# Finish off

if(NOT CMAKE_REQUIRED_QUIET)
	message("#\n# Configuring ${PROJECT_NAME} done\n")
endif()