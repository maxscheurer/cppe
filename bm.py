import os
import numpy as np
from pyscf import gto
import cppe
import yaml


def print_callback(string):
    print(string)


shells = np.arange(10, 65, 5, dtype=int)

folder = "../fmm_benchmarks/"

ret = {}
for shell in shells[3:4]:
    res = {}
    xyzfile = os.path.join(folder, f"solvated_{shell}", "pna.xyz")
    with open(xyzfile, "r") as f:
        xyz = f.readlines()[2:]
        xyz = "\n".join([x.strip() for x in xyz])
    mol = gto.Mole(atom=xyz, basis="sto3g")
    mol.build()

    cppe_mol = cppe.Molecule()
    for z, coord in zip(mol.atom_charges(), mol.atom_coords()):
        cppe_mol.append(cppe.Atom(z, *coord))

    options = {
        "potfile": os.path.join(folder, f"solvated_{shell}",
                                f"solvated_{shell}.pot"),
        "induced_thresh": 1e-8,
    }
    cppe_state = cppe.CppeState(options, cppe_mol, print_callback)
    potentials = cppe_state.potentials
    npolsites = cppe_state.get_polarizable_site_number()
    print(f"{npolsites} polarizable sites.")

    cppe_state.calculate_static_energies_and_fields()
    # res["time_static"] = time_static

    static_fields = cppe_state.static_fields
    res["n_polsites"] = npolsites

    options['summation_induced_fields'] = "direct"
    indmom = cppe.InducedMoments(potentials, options)
    indmom_exact = indmom.compute_cg(static_fields.flatten())

    options['summation_induced_fields'] = "fmm"
    indmom = cppe.InducedMoments(potentials, options)
    indmom_fmm = indmom.compute_cg(static_fields.flatten())
    # res["time_fmm"] = time_fmm

    absdiff = np.abs(indmom_exact - indmom_fmm)
    rnorm = np.linalg.norm(indmom_exact - indmom_fmm)
    # print("rnorm", rnorm)
    # res["time_exact"] = time_exact
    res["rnorm"] = float(rnorm)
    res["max_absdiff"] = float(np.max(absdiff))

    ret[f"{shell}"] = res

with open("results.yml", "w") as f:
    yaml.safe_dump(ret, f)
