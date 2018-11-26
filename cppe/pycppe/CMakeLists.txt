include(FindPythonLibs)
find_package(pybind11 CONFIG)
if(${pybind11_FOUND})
    # get_property(_loc TARGET pybind11 PROPERTY LOCATION)
    # message(STATUS "${Cyan}Found pybind11${ColourReset}: ${_loc} (found version ${pybind11_VERSION})")
    add_library(pybind11_external INTERFACE)  # dummy
    message(STATUS "found pybind11 " ${pybind11_DIR})
else()
    include(ExternalProject)
    message(STATUS "Suitable pybind11 could not be located, ${Magenta}Building pybind11${ColourReset} instead.")
    ExternalProject_Add(
      pybind11_external
      GIT_REPOSITORY "https://github.com/pybind/pybind11.git"
      GIT_TAG "master"
      CMAKE_ARGS -DPYBIND11_TEST=OFF
                 -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
    
      set(pybind11_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/pybind11 CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
    )
endif()

message(STATUS "pybind dir " ${pybind11_DIR})
# pybind11_add_module(pycppe export_cppe.cc)

add_library(pycppe MODULE export_core.cc export_utils.cc pycppe.cc)

target_include_directories(pycppe PRIVATE ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(pycppe PRIVATE ${ARMADILLO_LIBRARIES} cppe)

target_link_libraries(pycppe PUBLIC pybind11::module)
set_target_properties(pycppe PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")