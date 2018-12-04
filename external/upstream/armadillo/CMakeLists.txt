find_package(Armadillo)
if(${ARMADILLO_FOUND})
    add_library(armadillo_external INTERFACE) # Dummy
else()
    include(ExternalProject)
    message(STATUS "Armadillo not found on system, building Armadillo from source.")
    set(arma_tar "armadillo-9.200.5.tar.xz")
    set(CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
                   "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
                   "-DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}")
    ExternalProject_Add(armadillo_external
    	DOWNLOAD_COMMAND wget https://sourceforge.net/projects/arma/files/${arma_tar}
    	CONFIGURE_COMMAND tar -xf ../${arma_tar} -C ../armadillo_external  --strip-components=1
    	BUILD_IN_SOURCE 1
    	BUILD_COMMAND cmake ${CMAKE_ARGS} .
    	INSTALL_COMMAND make && make install
	)

    set(Armadillo_DIR ${STAGED_INSTALL_PREFIX}/share/Armadillo/CMake CACHE PATH "ArmadilloConfig.cmake directory" FORCE)
endif()