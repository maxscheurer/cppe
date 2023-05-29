target_sources(cppe-objs
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/calculate.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/field_m0.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/field_m1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/field_m2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tree.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/utils.cc
  )

list(APPEND _headers
  ${CMAKE_CURRENT_SOURCE_DIR}/calculate.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/field_m0.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/field_m1.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/field_m2.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/tree.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/utils.hh
  )

foreach(_header IN LISTS _headers)
  install(
    FILES
      ${_header}
    DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/core/fmm
    )
endforeach()
