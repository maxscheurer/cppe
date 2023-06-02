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

add_subdirectory(tensors)
add_subdirectory(fmm)
