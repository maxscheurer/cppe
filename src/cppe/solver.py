from __future__ import annotations

import numpy as np
from numba import njit, prange


def _triangle_to_mat(values):
    ret = np.zeros((3, 3), dtype=np.float64)
    tri = np.triu_indices(3)
    ret[tri] = values
    ret[(tri[1], tri[0])] = values
    return ret


@njit(cache=True, parallel=True)
def _induced_fields_direct_kernel(positions, include_mask, induced_moments):
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

            dx = positions[j, 0] - ix
            dy = positions[j, 1] - iy
            dz = positions[j, 2] - iz

            r2 = dx * dx + dy * dy + dz * dz
            inv_r = 1.0 / np.sqrt(r2)
            inv_r5 = (inv_r * inv_r * inv_r) / r2

            txx = (3.0 * dx * dx - r2) * inv_r5
            txy = 3.0 * dx * dy * inv_r5
            txz = 3.0 * dx * dz * inv_r5
            tyy = (3.0 * dy * dy - r2) * inv_r5
            tyz = 3.0 * dy * dz * inv_r5
            tzz = (3.0 * dz * dz - r2) * inv_r5

            j3 = 3 * j
            sx = induced_moments[j3 + 0]
            sy = induced_moments[j3 + 1]
            sz = induced_moments[j3 + 2]

            fx -= txx * sx + txy * sy + txz * sz
            fy -= txy * sx + tyy * sy + tyz * sz
            fz -= txz * sx + tyz * sy + tzz * sz

        i3 = 3 * i
        out[i3 + 0] = fx
        out[i3 + 1] = fy
        out[i3 + 2] = fz

    return out


@njit(cache=True, parallel=True)
def _apply_diagonal_kernel(alpha_inverse, induced_moments):
    n_sites = alpha_inverse.shape[0]
    out = np.zeros(3 * n_sites, dtype=np.float64)
    for i in prange(n_sites):
        i3 = 3 * i
        sx = induced_moments[i3 + 0]
        sy = induced_moments[i3 + 1]
        sz = induced_moments[i3 + 2]
        a = alpha_inverse[i]
        out[i3 + 0] = a[0, 0] * sx + a[0, 1] * sy + a[0, 2] * sz
        out[i3 + 1] = a[1, 0] * sx + a[1, 1] * sy + a[1, 2] * sz
        out[i3 + 2] = a[2, 0] * sx + a[2, 1] * sy + a[2, 2] * sz
    return out


@njit(cache=True, parallel=True)
def _apply_diagonal_inverse_kernel(alpha, in_vector):
    n_sites = alpha.shape[0]
    out = np.zeros(3 * n_sites, dtype=np.float64)
    for i in prange(n_sites):
        i3 = 3 * i
        sx = in_vector[i3 + 0]
        sy = in_vector[i3 + 1]
        sz = in_vector[i3 + 2]
        a = alpha[i]
        out[i3 + 0] = a[0, 0] * sx + a[0, 1] * sy + a[0, 2] * sz
        out[i3 + 1] = a[1, 0] * sx + a[1, 1] * sy + a[1, 2] * sz
        out[i3 + 2] = a[2, 0] * sx + a[2, 1] * sy + a[2, 2] * sz
    return out


class BMatrix:
    def __init__(self, polsites, options):
        self._polsites = polsites
        self._options = dict(options)
        self._n_polsites = len(polsites)

        self._positions = np.empty((self._n_polsites, 3), dtype=np.float64)
        self._alpha = np.empty((self._n_polsites, 3, 3), dtype=np.float64)
        self._alpha_inverse = np.empty((self._n_polsites, 3, 3), dtype=np.float64)
        self._include_mask = np.zeros(
            (self._n_polsites, self._n_polsites), dtype=np.int8
        )

        for i, pot in enumerate(polsites):
            self._positions[i, 0] = pot.x
            self._positions[i, 1] = pot.y
            self._positions[i, 2] = pot.z

            alpha = _triangle_to_mat(
                np.asarray(pot.polarizability.values, dtype=np.float64)
            )
            self._alpha[i] = alpha
            self._alpha_inverse[i] = np.linalg.inv(alpha)

            for j, pot_other in enumerate(polsites):
                if i == j:
                    continue
                if pot.excludes_site(pot_other.index):
                    continue
                self._include_mask[i, j] = 1

        scheme = self._options.get("summation_induced_fields", "direct")
        if scheme != "direct":
            raise NotImplementedError(
                "Python BMatrix currently supports summation_induced_fields='direct' only."
            )
        if self._options.get("damp_induced", False):
            raise NotImplementedError(
                "Python BMatrix currently does not support damp_induced=True."
            )

    def apply(self, induced_moments):
        vec = np.asarray(induced_moments, dtype=np.float64)
        induced_fields = _induced_fields_direct_kernel(
            self._positions, self._include_mask, vec
        )
        return induced_fields + _apply_diagonal_kernel(self._alpha_inverse, vec)

    def apply_diagonal_inverse(self, in_vector):
        vec = np.asarray(in_vector, dtype=np.float64)
        return _apply_diagonal_inverse_kernel(self._alpha, vec)

    def apply_diagonal(self, in_vector):
        vec = np.asarray(in_vector, dtype=np.float64)
        return _apply_diagonal_kernel(self._alpha_inverse, vec)

    def to_dense_matrix(self):
        n = self._n_polsites
        dense = np.zeros((3 * n, 3 * n), dtype=np.float64)
        for i in range(n):
            i3 = 3 * i
            dense[i3 : i3 + 3, i3 : i3 + 3] = self._alpha_inverse[i]
            xi, yi, zi = self._positions[i]
            for j in range(n):
                if self._include_mask[i, j] == 0:
                    continue
                j3 = 3 * j
                dx = self._positions[j, 0] - xi
                dy = self._positions[j, 1] - yi
                dz = self._positions[j, 2] - zi
                r2 = dx * dx + dy * dy + dz * dz
                inv_r = 1.0 / np.sqrt(r2)
                inv_r5 = (inv_r * inv_r * inv_r) / r2
                txx = (3.0 * dx * dx - r2) * inv_r5
                txy = 3.0 * dx * dy * inv_r5
                txz = 3.0 * dx * dz * inv_r5
                tyy = (3.0 * dy * dy - r2) * inv_r5
                tyz = 3.0 * dy * dz * inv_r5
                tzz = (3.0 * dz * dz - r2) * inv_r5
                dense[i3 + 0, j3 + 0] = -txx
                dense[i3 + 0, j3 + 1] = -txy
                dense[i3 + 0, j3 + 2] = -txz
                dense[i3 + 1, j3 + 0] = -txy
                dense[i3 + 1, j3 + 1] = -tyy
                dense[i3 + 1, j3 + 2] = -tyz
                dense[i3 + 2, j3 + 0] = -txz
                dense[i3 + 2, j3 + 1] = -tyz
                dense[i3 + 2, j3 + 2] = -tzz
        return dense

    def direct_inverse(self):
        return np.linalg.inv(self.to_dense_matrix())


class InducedMoments:
    def __init__(self, potentials, options):
        self._potentials = list(potentials)
        self._polsites = [p for p in self._potentials if p.is_polarizable]
        self._options = dict(options)
        self._printer = lambda s: None

    def set_print_callback(self, printer):
        self._printer = printer

    def compute(self, rhs, guess_or_make_guess, make_guess=None):
        rhs_vec = np.asarray(rhs, dtype=np.float64)
        if make_guess is None:
            make_guess = bool(guess_or_make_guess)
            guess = np.zeros_like(rhs_vec)
        else:
            guess = np.asarray(guess_or_make_guess, dtype=np.float64)

        bmat = BMatrix(self._polsites, self._options)

        if make_guess:
            x = bmat.apply_diagonal_inverse(rhs_vec)
        else:
            x = guess.copy()

        r = rhs_vec - bmat.apply(x)
        z = bmat.apply_diagonal_inverse(r)
        p = z.copy()

        maxiter = int(self._options.get("maxiter", 50))
        thresh = float(self._options.get("induced_thresh", 1e-8))

        for k in range(maxiter):
            ap = bmat.apply(p)
            alpha = float(np.dot(r, z) / np.dot(p, ap))
            x = x + alpha * p
            r_next = r - alpha * ap

            rnorm = float(np.linalg.norm(r_next))
            self._printer(f"{k} --- Norm: {rnorm:.12f}")
            if rnorm < thresh:
                return x

            z_next = bmat.apply_diagonal_inverse(r_next)
            beta = float(np.dot(z_next, r_next) / np.dot(z, r))
            p = z_next + beta * p
            r = r_next
            z = z_next

        raise RuntimeError("Failed to converge induced dipole moments.")
