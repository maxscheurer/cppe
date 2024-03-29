# set minimum cmake version
cmake_minimum_required(VERSION 3.16...3.20)

# project name
project(cppe LANGUAGES CXX VERSION 0.3.2)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# if CMAKE_BUILD_TYPE undefined, we set it to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# options handling utilities
include(CMakeDependentOption)
# Macro for printing an option in a consistent manner
# Syntax: print_option(<option to print> <was specified>)
function(print_option variable default)
  if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    message(STATUS "Setting (unspecified) option ${variable}: ${default}")
  else()
    message(STATUS "Setting option ${variable}: ${${variable}}")
  endif()
endfunction()

# Wraps an option with default ON/OFF. Adds nice messaging to option()
# Syntax: option_with_print(<option name> <description> <default value>)
macro(option_with_print variable msge default)
  print_option(${variable} ${default})
  option(${variable} ${msge} ${default})
endmacro()

# Wraps an option with a default other than ON/OFF and prints it
# Syntax: option_with_default(<option name> <description> <default value>)
macro(option_with_default variable msge default)
  # NOTE
  #  1. Can't combine with above b/c CMake handles ON/OFF options specially
  #  2. CMake variables are always defined so need to further check for if
  #     they are the NULL string. This is also why we need the force
  print_option(${variable} "${default}")
  if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    set(${variable} "${default}" CACHE STRING ${msge} FORCE)
  endif()
endmacro()

option_with_print(ENABLE_OPENMP "Enables OpenMP parallelization" ON)

# included cmake modules
include(${PROJECT_SOURCE_DIR}/cmake/downloaded/autocmake_default_build_paths.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/downloaded/autocmake_safeguards.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/compiler_flags/CXXFlags.cmake)

add_subdirectory(src)

include(${PROJECT_SOURCE_DIR}/cmake/downloaded/autocmake_save_flags.cmake)
