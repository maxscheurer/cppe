if(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  set(CPPE_CXX_FLAGS "/W3 /EHsc /bigobj")
  set(CMAKE_CXX_FLAGS_RELEASE "/O2")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2")
  set(CMAKE_CXX_FLAGS_DEBUG "/Od /W4")
endif()