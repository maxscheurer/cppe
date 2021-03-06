cmake_minimum_required(VERSION 3.0.0)
cmake_policy(VERSION 3.0.0)

project(cppe LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(psi4OptionsTools)
option_with_print(BUILD_SHARED_LIBS "Build shared libraries instead of static ones (CMake default is static)" ON)
option_with_print(ENABLE_OPENMP "Enables OpenMP parallelization" ON)
option_with_print(ENABLE_PYTHON_INTERFACE "Build the Python interface" OFF)
option_with_print(INSTALL_DEVEL_HEADERS "Install devel headers" OFF)
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization (with MSVC, it enables AVX2 instructions)" ON
                  "-xHost" "-march=native" "/arch:AVX2")
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
option_with_default(BUILD_FPIC "Compile static libraries with position independent code" ON)
option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed" /)
option_with_default(cppe_CXX_STANDARD "Specify C++ standard for core cppe" 14)

##########################################################################
# Project setup

if(NOT CMAKE_REQUIRED_QUIET)
	message("\n# Configuring ${PROJECT_NAME}\n#")
endif()

# Export compile commands if cmake supports it.
if (CMAKE_VERSION VERSION_GREATER 3.5.0)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)

set(pybind_ext "")
set(python_flags "")
if(${ENABLE_PYTHON_INTERFACE})
    set(Python_ADDITIONAL_VERSIONS 3.8 3.7 3.6)  # adjust with CMake minimum FindPythonInterp
    find_package(PythonLibsNew 3.6 REQUIRED)
    message(STATUS "${Cyan}Found Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE} (found version ${PYTHON_VERSION_STRING})")
    set(pybind_ext "pybind11_external")
    set(python_flags
              -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
              -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
              -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
        )
endif()

include(GNUInstallDirs)
include(autocmake_safeguards)
include(custom_color_messages)
include(ExternalProject)

add_subdirectory(external/upstream)

ExternalProject_Add(cppe-core
   DEPENDS ${pybind_ext}
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cppe
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
              -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
              -DCMAKE_CXX_STANDARD=${cppe_CXX_STANDARD}
              -DCMAKE_CXX_STANDARD_REQUIRED=ON
              -DEIGEN3_INCLUDE_DIRS=${PROJECT_SOURCE_DIR}/external/eigen3
              -DENABLE_OPENMP=${ENABLE_OPENMP}
              -DENABLE_PYTHON_INTERFACE=${ENABLE_PYTHON_INTERFACE}
              ${python_flags}
              -Dpybind11_DIR=${pybind11_DIR}
              -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
              -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
              -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
              -DINSTALL_DEVEL_HEADERS=${INSTALL_DEVEL_HEADERS}
              -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
              -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
              -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
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
