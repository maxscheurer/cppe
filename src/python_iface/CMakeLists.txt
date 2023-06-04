find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule
  )

find_package(pybind11 2.6 CONFIG REQUIRED)

pybind11_add_module(pycppe
  MODULE
    export_molecule.cc
    export_potential.cc
    export_fields.cc
    export_utils.cc
    export_state.cc
    export_math.cc
    export_tensors.cc
    export_fmm.cc
    export_cppe.cc
)

target_link_libraries(pycppe
  PRIVATE
    cppe-objs
)

install(
  TARGETS
    pycppe
  LIBRARY DESTINATION
    cppe
  )
