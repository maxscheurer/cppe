import os

import numpy as np
from scipy.sparse import linalg

from cppe import CppeState, BMatrix

from cppe import get_polarizable_sites, InducedMoments, PotfileReader, fmm, MultipoleFields
from cppe.tensors import T_recursive

from scipy.sparse.linalg import LinearOperator

from timings import Timer

# fl = "../mpi_cppe/nilered_in_blg_1000WAT.pot"
# fl = "../mpi_cppe/GFP_cutoff_5.pot"
fl = "../mpi_cppe/pft.pot"
potentials = PotfileReader(fl).read()
polsites = get_polarizable_sites(potentials)
n_polsites = len(polsites)

bmatrix_cpp = BMatrix(polsites, {'tree_expansion_order': 7})
print(n_polsites)

timer = Timer()

mf = MultipoleFields(potentials, {})
with timer.record("multipole_fields"):
    multipole_fields = mf.compute()

tot_field = multipole_fields

# TODO: move to C++
positions = np.array([
    p.position for p in polsites
], order="C")

# TODO: move to C++
exclusion_lists = []
for i, p1 in enumerate(polsites):
    lst = []
    for j, p2 in enumerate(polsites):
        if p1.excludes_site(p2.index):
            lst.append(j)
    exclusion_lists.append(lst)

ncrit = 64
order = 7
theta = 0.3
# tree = fmm.build_tree(positions, np.zeros_like(tot_field).flatten(),
#                       n_polsites, ncrit, order, theta,
#                       exclusion_lists)


def field_fmm(indmom):
    sources = indmom.flatten(order="C")
    tree = fmm.build_tree(positions, sources,
                          n_polsites, ncrit, order, theta,
                          exclusion_lists)
    # tree.set_sources(sources)
    field = np.zeros_like(sources)
    tree.compute_field_fmm(field)
    return field.reshape(indmom.shape)


def apply_bmatrix_fmm(indmom):
    induced_fields = field_fmm(indmom)
    return induced_fields + bmatrix_cpp.apply_diagonal(indmom)


# run the solver
rhs = tot_field.flatten()
n_iter = 0

ref_apply = bmatrix_cpp.apply(rhs)
test_apply = apply_bmatrix_fmm(rhs)
test_apply_cpp = bmatrix_cpp.apply_fast_summation(rhs, "fmm")
# np.testing.assert_allclose(ref_apply, test_apply, atol=1e-10)
print("maxdiff", np.max(np.abs(ref_apply - test_apply)))
np.testing.assert_allclose(test_apply, test_apply_cpp)


def callback(xk):
    global n_iter
    n_iter += 1
    print("--- iteration", n_iter)


linear_Pinv = LinearOperator(2 * (tot_field.size,),
                             matvec=bmatrix_cpp.apply_diagonal_inverse)
x0 = linear_Pinv @ rhs
bmat_operator = LinearOperator(2 * (tot_field.size,),
                               matvec=bmatrix_cpp.apply)

with timer.record("direct"):
    solution, status = linalg.cg(bmat_operator, rhs, x0, tol=1e-8,
                                 M=linear_Pinv, callback=callback)
    # indmom = InducedMoments(potentials, {})
    # solution = indmom.compute(rhs, True)
print(n_polsites, n_iter)
n_iter = 0


bmat_fmm = LinearOperator(2 * (tot_field.size,),
                          matvec=apply_bmatrix_fmm)

solver = linalg.cg
with timer.record("fmm"):
    solution_fmm, status = solver(bmat_fmm, rhs, x0, M=linear_Pinv,
                                  tol=1e-8, callback=callback)
print(n_polsites, n_iter)

maxdiff = np.max(np.abs(solution_fmm - solution))
rnorm = np.linalg.norm(solution_fmm - solution)
print(f"maxdiff {maxdiff}, rnorm {rnorm}")

print(timer.describe())
