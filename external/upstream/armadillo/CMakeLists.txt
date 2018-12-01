# TODO: find out required min. version

# find_package(Armadillo)
# 
# if(${ARMADILLO_FOUND})
#     add_library(armadillo_external INTERFACE) # Dummy
#     set(ARMADILLO_INCLUDE_DIRS_IFC ${ARMADILLO_INCLUDE_DIRS})
# else()
#     message(STATUS "${CYAN}Armadillo not found on system, building Armadillo from source.")
#     set(arma_tar "armadillo-9.200.5.tar.xz")
#     set(CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
#     "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
#     "-DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}"
#     "-DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}")
#     ExternalProject_Add(armadillo_external
# 	DOWNLOAD_COMMAND wget https://sourceforge.net/projects/arma/files/${arma_tar}
# 	CONFIGURE_COMMAND tar -xf ../${arma_tar} -C ../armadillo_external  --strip-components=1
# 	BUILD_IN_SOURCE 1
# 	BUILD_COMMAND cmake ${CMAKE_ARGS} .
# 	INSTALL_COMMAND make && make install
# 	)
# 
#     set(Armadillo_DIR ${CMAKE_BINARY_DIR}/cppe/armadillo_external-prefix/src/armadillo_external
#     ${STAGED_INSTALL_PREFIX}/share/cmake/pybind11 CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
#     # set(ARMADILLO_LIBRARIES ${Armadillo_DIR}/libarmadillo${CMAKE_SHARED_LIBRARY_SUFFIX})
#     # set(ARMADILLO_INCLUDE_DIRS ${Armadillo_DIR}/include)
#     set(ARMADILLO_INCLUDE_DIRS_IFC "")
# 
# endif()