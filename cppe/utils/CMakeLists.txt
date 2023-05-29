target_sources(cppe
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/string_utils.cc
  )

list(APPEND _headers
  ${CMAKE_CURRENT_SOURCE_DIR}/string_utils.hh
  )

foreach(_header IN LISTS _headers)
  install(
    FILES
      ${_header}
    DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils
    )
endforeach()
