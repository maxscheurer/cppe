find_package(gen1int)
if(${gen1int_FOUND})
    add_library(gen1int_external INTERFACE)
    message(STATUS "gen1int found here: " ${gen1int_DIR})
else()
    include(ExternalProject)
    message(STATUS "Suitable gen1int could not be located, ${Magenta}Building gen1int${ColourReset} instead.")
    ExternalProject_Add(
      gen1int_external
      GIT_REPOSITORY "https://gitlab.com/bingao/gen1int.git"
      GIT_TAG "3529af7"
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      INSTALL_COMMAND cp libgen1int.a ${STAGED_INSTALL_PREFIX}/lib &&
                      mkdir -p ${STAGED_INSTALL_PREFIX}/modules &&
                      cp modules/gen1int.mod ${STAGED_INSTALL_PREFIX}/modules &&
                      cp modules/gen1int_geom.mod ${STAGED_INSTALL_PREFIX}/modules &&
                      cp modules/london_ao.mod ${STAGED_INSTALL_PREFIX}/modules
    )
    # set(gen1int_DIR ${STAGED_INSTALL_PREFIX}/lib CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
    # file(COPY  DESTINATION ${STAGED_INSTALL_PREFIX}/lib)
    # find_file(gen1int_libfile libgen1int.a)
    # message(FATAL_ERROR "gen1int_DIR = " ${gen1int_DIR})
endif()