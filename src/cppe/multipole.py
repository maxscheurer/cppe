from __future__ import annotations

import numpy as np
from numba import njit, prange


@njit(cache=True, parallel=True)
def _interaction_energy_kernel(
    site_positions,
    atom_positions,
    atom_charges,
    charges,
    dipoles,
    quadrupoles,
    use_dipoles,
    use_quadrupoles,
):
    n_sites = site_positions.shape[0]
    n_atoms = atom_positions.shape[0]
    site_energies = np.zeros(n_sites, dtype=np.float64)

    for i in prange(n_sites):
        sx = site_positions[i, 0]
        sy = site_positions[i, 1]
        sz = site_positions[i, 2]

        q0 = charges[i]
        mux = dipoles[i, 0]
        muy = dipoles[i, 1]
        muz = dipoles[i, 2]

        qxx = quadrupoles[i, 0]
        qxy = quadrupoles[i, 1]
        qxz = quadrupoles[i, 2]
        qyy = quadrupoles[i, 3]
        qyz = quadrupoles[i, 4]
        qzz = quadrupoles[i, 5]

        acc = 0.0
        for a in range(n_atoms):
            dx = atom_positions[a, 0] - sx
            dy = atom_positions[a, 1] - sy
            dz = atom_positions[a, 2] - sz

            r2 = dx * dx + dy * dy + dz * dz
            inv_r = 1.0 / np.sqrt(r2)
            inv_r3 = inv_r / r2

            e = q0 * inv_r
            if use_dipoles:
                e += (mux * dx + muy * dy + muz * dz) * inv_r3

            if use_quadrupoles:
                inv_r5 = inv_r3 / r2
                txx = (3.0 * dx * dx - r2) * inv_r5
                txy = 3.0 * dx * dy * inv_r5
                txz = 3.0 * dx * dz * inv_r5
                tyy = (3.0 * dy * dy - r2) * inv_r5
                tyz = 3.0 * dy * dz * inv_r5
                tzz = (3.0 * dz * dz - r2) * inv_r5
                e += (
                    0.5 * qxx * txx
                    + qxy * txy
                    + qxz * txz
                    + 0.5 * qyy * tyy
                    + qyz * tyz
                    + 0.5 * qzz * tzz
                )

            acc += atom_charges[a] * e
        site_energies[i] = acc

    return site_energies.sum()


@njit(cache=True, parallel=True)
def _nuclear_gradient_kernel(
    site_positions,
    atom_positions,
    atom_charges,
    charges,
    dipoles,
    quadrupoles,
    has_dipole,
    has_quadrupole,
):
    n_atoms = atom_positions.shape[0]
    n_sites = site_positions.shape[0]
    grad = np.zeros((n_atoms, 3), dtype=np.float64)

    for ai in prange(n_atoms):
        ax = atom_positions[ai, 0]
        ay = atom_positions[ai, 1]
        az = atom_positions[ai, 2]
        q_atom = atom_charges[ai]

        gx = 0.0
        gy = 0.0
        gz = 0.0

        for j in range(n_sites):
            dx = ax - site_positions[j, 0]
            dy = ay - site_positions[j, 1]
            dz = az - site_positions[j, 2]

            r2 = dx * dx + dy * dy + dz * dz
            inv_r = 1.0 / np.sqrt(r2)
            inv_r3 = inv_r / r2

            q = charges[j]
            gx += q * dx * inv_r3
            gy += q * dy * inv_r3
            gz += q * dz * inv_r3

            if has_dipole[j]:
                mx = dipoles[j, 0]
                my = dipoles[j, 1]
                mz = dipoles[j, 2]
                mdotr = mx * dx + my * dy + mz * dz
                inv_r5 = inv_r3 / r2
                gx += 3.0 * dx * mdotr * inv_r5 - mx * inv_r3
                gy += 3.0 * dy * mdotr * inv_r5 - my * inv_r3
                gz += 3.0 * dz * mdotr * inv_r5 - mz * inv_r3

            if has_quadrupole[j]:
                qxx = quadrupoles[j, 0]
                qxy = quadrupoles[j, 1]
                qxz = quadrupoles[j, 2]
                qyy = quadrupoles[j, 3]
                qyz = quadrupoles[j, 4]
                qzz = quadrupoles[j, 5]
                qrx = qxx * dx + qxy * dy + qxz * dz
                qry = qxy * dx + qyy * dy + qyz * dz
                qrz = qxz * dx + qyz * dy + qzz * dz
                rqrr = dx * qrx + dy * qry + dz * qrz
                inv_r5 = inv_r3 / r2
                inv_r7 = inv_r5 / r2
                gx += 7.5 * dx * rqrr * inv_r7 - 3.0 * qrx * inv_r5
                gy += 7.5 * dy * rqrr * inv_r7 - 3.0 * qry * inv_r5
                gz += 7.5 * dz * rqrr * inv_r7 - 3.0 * qrz * inv_r5

        grad[ai, 0] = -q_atom * gx
        grad[ai, 1] = -q_atom * gy
        grad[ai, 2] = -q_atom * gz

    return grad


class MultipoleExpansion:
    def __init__(self, molecule, potentials):
        self._molecule = molecule
        self._potentials = potentials

        n_atoms = len(molecule)
        n_sites = len(potentials)
        self._atom_positions = np.empty((n_atoms, 3), dtype=np.float64)
        self._atom_charges = np.empty(n_atoms, dtype=np.float64)
        for i, atom in enumerate(molecule):
            self._atom_positions[i, 0] = atom.x
            self._atom_positions[i, 1] = atom.y
            self._atom_positions[i, 2] = atom.z
            self._atom_charges[i] = atom.charge

        self._site_positions = np.empty((n_sites, 3), dtype=np.float64)
        self._charges = np.zeros(n_sites, dtype=np.float64)
        self._dipoles = np.zeros((n_sites, 3), dtype=np.float64)
        self._quadrupoles = np.zeros((n_sites, 6), dtype=np.float64)

        self._max_order = 0
        self._has_dipole = np.zeros(n_sites, dtype=np.int8)
        self._has_quadrupole = np.zeros(n_sites, dtype=np.int8)
        for i, site in enumerate(potentials):
            self._site_positions[i, 0] = site.x
            self._site_positions[i, 1] = site.y
            self._site_positions[i, 2] = site.z
            for multipole in site.multipoles:
                k = multipole.k
                if k > self._max_order:
                    self._max_order = k
                values = np.asarray(multipole.values, dtype=np.float64)
                if k == 0:
                    self._charges[i] = values[0]
                elif k == 1:
                    self._has_dipole[i] = 1
                    self._dipoles[i, :] = values
                elif k == 2:
                    self._has_quadrupole[i] = 1
                    self._quadrupoles[i, :] = values
        if self._max_order > 2:
            raise NotImplementedError(
                "Python MultipoleExpansion currently supports multipoles up to order 2."
            )

    def interaction_energy(self):
        return _interaction_energy_kernel(
            self._site_positions,
            self._atom_positions,
            self._atom_charges,
            self._charges,
            self._dipoles,
            self._quadrupoles,
            self._max_order > 0,
            self._max_order > 1,
        )

    def nuclear_gradient(self):
        return _nuclear_gradient_kernel(
            self._site_positions,
            self._atom_positions,
            self._atom_charges,
            self._charges,
            self._dipoles,
            self._quadrupoles,
            self._has_dipole,
            self._has_quadrupole,
        )
