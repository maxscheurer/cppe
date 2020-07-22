import os
import numpy as np
from pyscf import gto
import cppe
import yaml
from timings import Timer


def print_callback(string):
    print(string)


shells = np.arange(10, 65, 5, dtype=int)
print(shells)
folder = "../fmm_benchmarks/"

ret = {}
for shell in shells:
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

    timer = Timer()

    with timer.record("static_field_energies"):
        cppe_state.calculate_static_energies_and_fields()
        print("Static fields done.")

    static_fields = cppe_state.static_fields
    res["n_polsites"] = npolsites

    with timer.record("direct"):
        options['summation_induced_fields'] = "direct"
        indmom = cppe.InducedMoments(potentials, options)
        indmom_exact = indmom.compute_cg(static_fields.flatten())

    # with timer.record("bh"):
    #     options['summation_induced_fields'] = "bh"
    #     indmom = cppe.InducedMoments(potentials, options)
    #     indmom_exact = indmom.compute_cg(static_fields.flatten())

    with timer.record("fmm"):
        options['summation_induced_fields'] = "fmm"
        indmom = cppe.InducedMoments(potentials, options)
        indmom_fmm = indmom.compute_cg(static_fields.flatten())

    absdiff = np.abs(indmom_exact - indmom_fmm)
    rnorm = np.linalg.norm(indmom_exact - indmom_fmm)
    res["rnorm"] = float(rnorm)
    res["max_absdiff"] = float(np.max(absdiff))
    for t in timer.tasks:
        res[t] = timer.total(t)
    ret[f"{shell}"] = res

with open("results.yml", "w") as f:
    yaml.safe_dump(ret, f)
