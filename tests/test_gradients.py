import unittest
import os

import numpy as np

from .cache import cache

import cppe
from cppe import PotfileReader, MultipoleExpansion, NuclearFields


prefactors_5p = np.array([1.0, -8.0, 8.0, -1.0]) / 12.0
multipliers_5p = [-2, -1, 1, 2]
prefactors_2p = [-0.5, 0.5]
multipliers_2p = [-1, 1]
stencils = {
    "2p": (multipliers_2p, prefactors_2p),
    "5p": (multipliers_5p, prefactors_5p)
}
coords_label = ["x", "y", "z"]


def print_callback(output):
    print("cb: ", output)


class TestGradients(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = "{}/potfiles/pna_6w.pot".format(dirname)
    potfile_iso_path = "{}/potfiles/pna_6w_isopol.pot".format(dirname)

    def test_compute_nuclear_interaction_energy_gradient(self):
        mol = cache.molecule["pna"]
        options = {"potfile": self.potfile_path}

        p = PotfileReader(self.potfile_path)
        potentials = p.read()

        atoms = [a.atomic_number for a in mol]
        coords = np.array([a.position for a in mol])

        grad_nuc = MultipoleExpansion(mol, potentials).nuclear_gradient()

        grad_nuc_fd = grad_nuclear_interaction_energy_fdiff(atoms, coords, potentials)
        np.testing.assert_allclose(grad_nuc_fd, grad_nuc, atol=1e-10)
        
        grad_field_fd = grad_nuclear_field_fdiff(atoms, coords, potentials)

        natoms = len(atoms)
        polsites = cppe.get_polarizable_sites(potentials)
        grad_field = NuclearFields(mol, potentials).nuclear_gradient().reshape((natoms, 3, len(polsites), 3))
        np.testing.assert_allclose(grad_field_fd, grad_field, atol=1e-10)

def grad_nuclear_field_fdiff(atoms, coords, potentials, step_au=1e-3):
    """
    Computes the finite difference gradient
    with a 5-point stencil
    """
    natoms = len(atoms)
    polsites = cppe.get_polarizable_sites(potentials)
    grad = np.zeros((natoms, 3, len(polsites), 3))
    for i in range(natoms):
        for c in range(3):
            print("Computing dE/d{} for atom {}.".format(coords_label[c], i))
            for f, p in zip(*stencils["5p"]):
                coords_p = coords.copy()
                coords_p[i, c] += f * step_au
                mol_p = cppe.Molecule()
                for z, coord in zip(atoms, coords_p):
                    mol_p.append(cppe.Atom(z, *coord))

                nf = NuclearFields(mol_p, potentials)
                en_pert = nf.compute().reshape(len(polsites), 3)
                grad[i, c] += p * en_pert / step_au
    return grad


def grad_nuclear_interaction_energy_fdiff(atoms, coords, potentials, step_au=1e-3):
    """
    Computes the finite difference gradient
    with a 5-point stencil
    """
    natoms = len(atoms)
    grad = np.zeros((natoms, 3))
    for i in range(natoms):
        for c in range(3):
            print("Computing dE/d{} for atom {}.".format(coords_label[c], i))
            for f, p in zip(*stencils["5p"]):
                coords_p = coords.copy()
                coords_p[i, c] += f * step_au
                mol_p = cppe.Molecule()
                for z, coord in zip(atoms, coords_p):
                    mol_p.append(cppe.Atom(z, *coord))

                mexp = MultipoleExpansion(mol_p, potentials)
                en_pert = mexp.interaction_energy()
                grad[i, c] += p * en_pert / step_au
    return grad
