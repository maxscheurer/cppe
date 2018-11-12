cmake_minimum_required(VERSION 3.0.0)
cmake_policy(VERSION 3.0.0)

# Prefer shared to static libs (cmake does the reverse by default)
# option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones" ON)

##########################################################################
# Project setup

project(cppe VERSION 0.0.0 LANGUAGES Fortran C CXX)
message("\n# Configuring ${PROJECT_NAME}\n#")

##########################################################################
# Setup hard and optional dependencies and find components

# TODO setup

# Export compile commands if cmake supports it.
if (CMAKE_VERSION VERSION_GREATER 3.5.0)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

##########################################################################
# Pass on to subdirectories
# add_subdirectory(external)
add_subdirectory(cppe)

# TODO Testing

##########################################################################
# Finish off
message("#\n# Configuring ${PROJECT_NAME} done\n")
