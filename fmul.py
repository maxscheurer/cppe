import os
import numpy as np
from pyscf import gto
import cppe
import yaml
from timings import Timer
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

colors = sns.color_palette("Set2")

def setup():
    # Setup matplotlib
    tex_premable = [
        r"\usepackage[T1]{fontenc}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage{lmodern}",
        r"\usepackage{amsmath}",
    ]
    pgf_with_rc_fonts = {
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": True,
        "text.latex.preamble": tex_premable,
        "pgf.rcfonts": False,
        "pgf.preamble": tex_premable,
    }
    matplotlib.rcParams.update(pgf_with_rc_fonts)

setup()
sns.set_context("paper")
# sns.set_palette("Set2")

def print_callback(string):
    print(string)


shells = np.arange(50, 55, 5, dtype=int)
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
    }
    cppe_state = cppe.CppeState(options, cppe_mol, print_callback)
    potentials = cppe_state.potentials
    timer = Timer()
    
    fmuls = cppe.MultipoleFields(potentials, options)
    print("multipole fields prepared")
    with timer.record("direct"):
        fs = fmuls.compute_tree()

    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
    for ii, theta in enumerate([0.3, 0.5, 0.7, 0.99]):
        for ax, exp_order in zip(axes, [3,5,7]):
            options = {
                "summation_induced_fields": "fmm",
                "theta": theta,
                "tree_expansion_order": exp_order,
            }

            fmuls_tree = cppe.MultipoleFields(potentials, options) 
            with timer.record(f"tree_{theta}_{exp_order}"):
                fs_tree = fmuls_tree.compute_tree()

            fs = fs.reshape(len(potentials), 3)
            fs_tree = fs_tree.reshape(fs.shape)
            err_mu_i = np.linalg.norm(fs - fs_tree, axis=1) / np.linalg.norm(fs, axis=1)

            sns.histplot(x=np.log10(err_mu_i), ax=ax, stat="probability", label=fr"$\theta = {theta}$", color=colors[ii])
    # ax2 = ax.twinx()
    # sns.histplot(x=np.log10(err_mu_i), ax=ax2, stat="density", cumulative=True,
    #              element="step", fill=False)
    # ax2.set_ylim([0, 2])
            ax.set_xlabel(r"$\log_{10} F^\mathrm{err}$")
            ax.set_xlim([-16, 0])
            ax.set_title(f"order = {exp_order}")
    print(timer.describe())

    plt.legend()
    plt.tight_layout()
    plt.show()