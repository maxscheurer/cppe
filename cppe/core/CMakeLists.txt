target_sources(cppe-objs
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/potfile_reader.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/pot_manipulation.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/cppe_state.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/multipole_expansion.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/electric_fields.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/math.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/bmatrix.cc
  )

list(APPEND _headers
  ${CMAKE_CURRENT_SOURCE_DIR}/bmatrix.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/cppe_state.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/electric_fields.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/math.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/molecule.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/multipole_expansion.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/pe_options.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/pot_manipulation.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/potential.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/potfile_reader.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/tensors.hh
  )

foreach(_header IN LISTS _headers)
  install(
    FILES
      ${_header}
    DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/core
    )
endforeach()

add_subdirectory(tensors)
add_subdirectory(fmm)
