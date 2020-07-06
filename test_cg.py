import os

import numpy as np
from scipy.sparse import linalg

from cppe import CppeState, BMatrix

from cppe import get_polarizable_sites, InducedMoments, PotfileReader
from cppe.tensors import T_recursive

from scipy.sparse.linalg import LinearOperator


def block(s1, s2):
    return (slice(3 * s1, 3 * s1 + 3), slice(3 * s2, 3 * s2 + 3))


def triangle_to_mat(y):
    mat = np.empty((3, 3))
    mat[0, 0] = y[0]
    mat[0, 1] = mat[1, 0] = y[1]
    mat[0, 2] = mat[2, 0] = y[2]
    mat[1, 1] = y[3]
    mat[1, 2] = mat[2, 1] = y[4]
    mat[2, 2] = y[5]
    return mat


fl = "../mpi_cppe/nilered_in_blg_1000WAT.pot"
# fl = "../mpi_cppe/GFP_cutoff_5.pot"
a = PotfileReader(fl)
potentials = a.read()
polsites = get_polarizable_sites(potentials)
n_polsites = len(polsites)
ind_moms = InducedMoments(potentials, {})
tot_field = 0.001 * np.ones((n_polsites, 3))

dump_bmatrix = True
if dump_bmatrix:
    bmatrix_cpp = BMatrix(polsites, {})
    # A = LinearOperator(2 * (tot_field.size,),
    #                    matvec=bmatrix_cpp.compute_apply)
    # Afull = A @ np.eye(tot_field.size)
    A_cpp = bmatrix_cpp.to_dense_matrix()
    np.save("bmat.npy", A_cpp)
    # np.testing.assert_allclose(A_cpp, Afull, atol=1e-12)

bmat = np.load("bmat.npy")
print(n_polsites, bmat.shape)

# kappa = np.linalg.cond(bmat)
# print("Condition number", kappa)

D_alpha = np.zeros_like(bmat)
for s1, pot1 in enumerate(polsites):
    pol = pot1.polarizability
    alpha = triangle_to_mat(pol.values)
    D_alpha[block(s1, s1)] = alpha

# DB = D_alpha @ bmat
# kappa_d = np.linalg.cond(DB)
# print("Condition number [Jacobi]", kappa_d)

B_local = np.zeros_like(bmat)
for s1, pot1 in enumerate(polsites):
    for s2, pot2 in enumerate(polsites):
        if pot1.excludes_site(pot2.index) or s1 == s2:
            continue
        diff = pot2.position - pot1.position
        dist = np.linalg.norm(diff)
        if dist < 7.5 and s1 > s2:
            T12 = T_recursive(2, diff)
            B_local[block(s1, s2)] = -triangle_to_mat(T12)
            B_local[block(s2, s1)] = -triangle_to_mat(T12)
assert np.max(np.abs(B_local - B_local.T)) < 1e-9

# DB_poly = (D_alpha - D_alpha @ B_local @ D_alpha) @ bmat
# kappa_dpoly = np.linalg.cond(DB_poly)
# print("Condition number [Poly]", kappa_dpoly)

scipy_solvers = [
    # jacobi_diis,
    linalg.cg,
    # linalg.cgs,
    # linalg.gmres,
    # linalg.bicgstab
]

rhs = tot_field.flatten()
x0 = D_alpha @ rhs
n_iter = 0


def callback(xk):
    global n_iter
    n_iter += 1


linear_Pinv = LinearOperator(2 * (tot_field.size,),
                             matvec=lambda xk: D_alpha @ xk)

loc = D_alpha @ B_local @ D_alpha
poly_Pinv = LinearOperator(2 * (tot_field.size,),
                           matvec=lambda xk: (D_alpha - loc) @ xk)
preconds = [None, linear_Pinv, poly_Pinv]
for solver in scipy_solvers:
    for M in preconds:
        solution, status = solver(bmat, rhs, x0, tol=1e-9,
                                  M=M, callback=callback)
        print(n_iter, M)
        n_iter = 0
    # np.testing.assert_allclose(solution, ref, atol=conv_tol)
