target_sources(cppe
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tensors_autogen.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tensors_recursive.cc
  )

list(APPEND _headers
  ${CMAKE_CURRENT_SOURCE_DIR}/tensors_autogen.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/tensors_recursive.hh
  )

foreach(_header IN LISTS _headers)
  install(
    FILES
      ${_header}
    DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/core/tensors
    )
endforeach()
