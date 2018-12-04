from cppe import cppe as pycppe
import numpy as np
p = pycppe.Potential(1.0, 2.0, 3.0, 0)
a = pycppe.PotfileReader("/Users/maxscheurer/testing/qc/pehf_no_oct.pot")

atom = pycppe.Atom(8, 0.0, 1.0, 0.0)

print(atom.get_position())

mol = pycppe.Molecule()
mol.append(atom)

sites = a.read()
sites[0].x = -10.0
for s in sites:
    print(s.x)

# potman = pycppe.PotManipulator(sites, mol)

options = pycppe.PeOptions()
options.potfile = "/Users/maxscheurer/testing/qc/pehf_no_oct.pot"
options.do_diis = False
# options.pe_border = True
# options.border_options.rmin = 3.4

state = pycppe.CppeState(options, mol)
state.calculate_static_energies_and_fields()

# generate random data
fields = np.zeros(30)
operator = np.random.rand(49).reshape(7, 7)
density = np.random.rand(49).reshape(7, 7)

state.update_induced_moments(fields, 0, False)
state.set_es_operator(operator)
state.update_energies(density)
ind_norm = np.linalg.norm(state.get_induced_moments().reshape((10, 3)), axis=1)
print(ind_norm)
state.print_summary()

# stand-alone induced moments
ind_moms = pycppe.InducedMoments(a.read(), options)
fields = np.random.rand(30)
ind = ind_moms.compute(fields, True)
ind_norm = np.linalg.norm(ind.reshape((10, 3)), axis=1)
print(ind_norm)

print(pycppe.prefactors(2))
