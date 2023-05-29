message(STATUS "Building Python interface for cppe.")

find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development)

# figure out where to put the Python module
if(NOT DEFINED PYMOD_INSTALL_FULLDIR)
  if(NOT WIN32)
    set(PYMOD_INSTALL_FULLDIR
          "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/cppe"
        CACHE STRING
          "Location under CMAKE_INSTALL_PREFIX to which Python modules are installed"
        FORCE
      )
  else()
    set(PYMOD_INSTALL_FULLDIR
          "Lib/site-packages/cppe"
        CACHE STRING
          "Location under CMAKE_INSTALL_PREFIX to which Python modules are installed"
        FORCE
      )
  endif()
endif()
message(STATUS "Setting PYMOD_INSTALL_FULLDIR: ${PYMOD_INSTALL_FULLDIR}")
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})

find_package(pybind11 2.6 CONFIG REQUIRED)

pybind11_add_module(pycppe
  MODULE
  THIN_LTO
    export_molecule.cc
    export_potential.cc
    export_fields.cc
    export_utils.cc
    export_state.cc
    export_math.cc
    export_tensors.cc
    export_fmm.cc
    export_cppe.cc
)

target_link_libraries(pycppe
  PRIVATE
    cppe
)

# handle RPATH
file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_RPATH "@loader_path/${_rel}")
else()
  set(_RPATH "\$ORIGIN/${_rel}")
endif()

set_target_properties(pycppe
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${_RPATH}${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH ON
    ARCHIVE_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
    LIBRARY_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
    RUNTIME_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  )

# we glob the Python files in this folder and let CMake add a rule such that
# the glob is repeated every time we rebuild.
# This is NOT RECOMMENDED by CMake
# (https://cmake.org/cmake/help/v3.16/command/file.html#filesystem) but you only
# live once!
file(
  GLOB
    _pys
  LIST_DIRECTORIES
    FALSE
  CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py
  )

# link the Python files under the build folder
foreach(_py IN LISTS _pys)
  get_filename_component(__py ${_py} NAME)
  file(
    CREATE_LINK
      ${_py}
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/${__py}
    COPY_ON_ERROR
    SYMBOLIC
    )
endforeach()

install(
  TARGETS
    pycppe
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )

install(
  FILES
    ${_pys}
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )

# install the Python tests
install(
  DIRECTORY
    ${PROJECT_SOURCE_DIR}/tests/
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}/tests
  PATTERN "CMakeLists.txt" EXCLUDE
  PATTERN "__pycache__" EXCLUDE
  )
