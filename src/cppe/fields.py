from __future__ import annotations

import numpy as np
from numba import njit, prange


@njit(cache=True, parallel=True)
def _nuclear_fields_kernel(site_positions, atom_positions, atom_charges):
    n_sites = site_positions.shape[0]
    n_atoms = atom_positions.shape[0]
    fields = np.zeros((n_sites, 3), dtype=np.float64)

    for i in prange(n_sites):
        sx = site_positions[i, 0]
        sy = site_positions[i, 1]
        sz = site_positions[i, 2]
        fx = 0.0
        fy = 0.0
        fz = 0.0
        for j in range(n_atoms):
            dx = sx - atom_positions[j, 0]
            dy = sy - atom_positions[j, 1]
            dz = sz - atom_positions[j, 2]
            r2 = dx * dx + dy * dy + dz * dz
            inv_r = 1.0 / np.sqrt(r2)
            inv_r3 = inv_r / r2
            q = atom_charges[j]
            fx += q * dx * inv_r3
            fy += q * dy * inv_r3
            fz += q * dz * inv_r3
        fields[i, 0] = fx
        fields[i, 1] = fy
        fields[i, 2] = fz

    return fields.reshape(3 * n_sites)


@njit(cache=True, parallel=True)
def _nuclear_gradient_kernel(site_positions, atom_positions, atom_charges):
    n_sites = site_positions.shape[0]
    n_atoms = atom_positions.shape[0]
    grad = np.zeros((3 * n_atoms, 3 * n_sites), dtype=np.float64)

    for i in prange(n_sites):
        sx = site_positions[i, 0]
        sy = site_positions[i, 1]
        sz = site_positions[i, 2]
        col = 3 * i

        for ai in range(n_atoms):
            dx = sx - atom_positions[ai, 0]
            dy = sy - atom_positions[ai, 1]
            dz = sz - atom_positions[ai, 2]

            r2 = dx * dx + dy * dy + dz * dz
            inv_r = 1.0 / np.sqrt(r2)
            inv_r5 = (inv_r * inv_r * inv_r) / r2
            q = atom_charges[ai]

            txx = (3.0 * dx * dx - r2) * inv_r5
            txy = 3.0 * dx * dy * inv_r5
            txz = 3.0 * dx * dz * inv_r5
            tyy = (3.0 * dy * dy - r2) * inv_r5
            tyz = 3.0 * dy * dz * inv_r5
            tzz = (3.0 * dz * dz - r2) * inv_r5

            row = 3 * ai
            grad[row + 0, col + 0] = q * txx
            grad[row + 0, col + 1] = q * txy
            grad[row + 0, col + 2] = q * txz

            grad[row + 1, col + 0] = q * txy
            grad[row + 1, col + 1] = q * tyy
            grad[row + 1, col + 2] = q * tyz

            grad[row + 2, col + 0] = q * txz
            grad[row + 2, col + 1] = q * tyz
            grad[row + 2, col + 2] = q * tzz

    return grad


class NuclearFields:
    def __init__(self, molecule, potentials):
        self._molecule = molecule
        self._polsites = [p for p in potentials if p.is_polarizable]

    @staticmethod
    def _as_atom_arrays(molecule):
        n_atoms = len(molecule)
        atom_positions = np.empty((n_atoms, 3), dtype=np.float64)
        atom_charges = np.empty(n_atoms, dtype=np.float64)
        for i, atom in enumerate(molecule):
            atom_positions[i, 0] = atom.x
            atom_positions[i, 1] = atom.y
            atom_positions[i, 2] = atom.z
            atom_charges[i] = atom.charge
        return atom_positions, atom_charges

    @staticmethod
    def _as_site_positions(polsites):
        n_sites = len(polsites)
        site_positions = np.empty((n_sites, 3), dtype=np.float64)
        for i, site in enumerate(polsites):
            site_positions[i, 0] = site.x
            site_positions[i, 1] = site.y
            site_positions[i, 2] = site.z
        return site_positions

    def compute(self):
        site_positions = self._as_site_positions(self._polsites)
        atom_positions, atom_charges = self._as_atom_arrays(self._molecule)
        return _nuclear_fields_kernel(site_positions, atom_positions, atom_charges)

    def nuclear_gradient(self):
        site_positions = self._as_site_positions(self._polsites)
        atom_positions, atom_charges = self._as_atom_arrays(self._molecule)
        return _nuclear_gradient_kernel(site_positions, atom_positions, atom_charges)
