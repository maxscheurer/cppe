from __future__ import annotations

import numpy as np
from numba import njit, prange


@njit(cache=True, parallel=True)
def _multipole_fields_direct_kernel(
    positions,
    include_mask,
    charges,
    dipoles,
    quadrupoles,
    has_dipole,
    has_quadrupole,
):
    n_sites = positions.shape[0]
    out = np.zeros(3 * n_sites, dtype=np.float64)

    for i in prange(n_sites):
        ix = positions[i, 0]
        iy = positions[i, 1]
        iz = positions[i, 2]
        fx = 0.0
        fy = 0.0
        fz = 0.0

        for j in range(n_sites):
            if include_mask[i, j] == 0:
                continue

            dx = ix - positions[j, 0]
            dy = iy - positions[j, 1]
            dz = iz - positions[j, 2]

            r2 = dx * dx + dy * dy + dz * dz
            inv_r = 1.0 / np.sqrt(r2)
            inv_r3 = inv_r / r2

            q = charges[j]
            fx += q * dx * inv_r3
            fy += q * dy * inv_r3
            fz += q * dz * inv_r3

            if has_dipole[j]:
                mx = dipoles[j, 0]
                my = dipoles[j, 1]
                mz = dipoles[j, 2]
                mdotr = mx * dx + my * dy + mz * dz
                inv_r5 = inv_r3 / r2
                fx += 3.0 * dx * mdotr * inv_r5 - mx * inv_r3
                fy += 3.0 * dy * mdotr * inv_r5 - my * inv_r3
                fz += 3.0 * dz * mdotr * inv_r5 - mz * inv_r3

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
                fx += 7.5 * dx * rqrr * inv_r7 - 3.0 * qrx * inv_r5
                fy += 7.5 * dy * rqrr * inv_r7 - 3.0 * qry * inv_r5
                fz += 7.5 * dz * rqrr * inv_r7 - 3.0 * qrz * inv_r5

        i3 = 3 * i
        out[i3 + 0] = fx
        out[i3 + 1] = fy
        out[i3 + 2] = fz

    return out


class MultipoleFields:
    def __init__(self, potentials, options):
        self._potentials = list(potentials)
        self._options = dict(options)
        if self._options.get("summation_induced_fields", "direct") != "direct":
            raise NotImplementedError(
                "Python MultipoleFields currently supports summation_induced_fields='direct' only."
            )
        if self._options.get("damp_multipole", False):
            raise NotImplementedError(
                "Python MultipoleFields currently does not support damp_multipole=True."
            )

    def compute(self):
        n_sites = len(self._potentials)
        positions = np.empty((n_sites, 3), dtype=np.float64)
        include_mask = np.zeros((n_sites, n_sites), dtype=np.int8)
        charges = np.zeros(n_sites, dtype=np.float64)
        dipoles = np.zeros((n_sites, 3), dtype=np.float64)
        quadrupoles = np.zeros((n_sites, 6), dtype=np.float64)
        has_dipole = np.zeros(n_sites, dtype=np.int8)
        has_quadrupole = np.zeros(n_sites, dtype=np.int8)

        for i, target in enumerate(self._potentials):
            positions[i, 0] = target.x
            positions[i, 1] = target.y
            positions[i, 2] = target.z

            for j, source in enumerate(self._potentials):
                if i == j:
                    continue
                if target.excludes_site(source.index):
                    continue
                include_mask[i, j] = 1

            for multipole in target.multipoles:
                values = np.asarray(multipole.values, dtype=np.float64)
                if multipole.k == 0:
                    charges[i] = values[0]
                elif multipole.k == 1:
                    has_dipole[i] = 1
                    dipoles[i, 0] = values[0]
                    dipoles[i, 1] = values[1]
                    dipoles[i, 2] = values[2]
                elif multipole.k == 2:
                    has_quadrupole[i] = 1
                    quadrupoles[i, :] = values

        return _multipole_fields_direct_kernel(
            positions,
            include_mask,
            charges,
            dipoles,
            quadrupoles,
            has_dipole,
            has_quadrupole,
        )
