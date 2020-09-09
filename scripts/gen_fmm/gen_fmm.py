import fmmgen

cse = True
atomic = True
precision = 'double'

order = 8
# order = 4
for k in range(3):
    print("Generating code for source order", k)
    # fmmgen.generate_code(order, f"pot_m{k}",
    #                     src_dir="../../cppe/core/fmm",
    #                     include_dir="../../cppe/core/fmm",
    #                     precision=precision,
    #                     CSE=cse,
    #                     cython=False,
    #                     potential=True,
    #                     field=False,
    #                     source_order=k,
    #                     atomic=atomic, minpow=11,
    #                     harmonic_derivs=True,
    #                     language='c++', name_prefix=f"pot_m{k}_",
    #                     write_defines=False)
    fmmgen.generate_code(order, f"field_m{k}",
                        src_dir="../../cppe/core/fmm",
                        include_dir="../../cppe/core/fmm",
                        precision=precision,
                        CSE=cse,
                        cython=False,
                        potential=False,
                        field=True,
                        source_order=k,
                        atomic=atomic, minpow=11,
                        harmonic_derivs=True,
                        language='c++', name_prefix=f"field_m{k}_",
                        write_defines=False, write_templates=True)