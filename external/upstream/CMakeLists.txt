if(${ENABLE_PYTHON_INTERFACE})
    add_subdirectory(pybind11)
endif()

if(${ENABLE_PELIB})
    file(MAKE_DIRECTORY ${STAGED_INSTALL_PREFIX})
    file(MAKE_DIRECTORY ${STAGED_INSTALL_PREFIX}/lib)
    file(MAKE_DIRECTORY ${STAGED_INSTALL_PREFIX}/modules)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${STAGED_INSTALL_PREFIX}/modules)
    add_subdirectory(gen1int)
    
    include(ExternalProject)
    ExternalProject_Add(
      gen1int_interface
      DEPENDS gen1int_external
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gen1int_interface
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                 -DCMAKE_Fortran_MODULE_DIRECTORY=${CMAKE_Fortran_MODULE_DIRECTORY}
    )
    add_subdirectory(pelib)
endif()