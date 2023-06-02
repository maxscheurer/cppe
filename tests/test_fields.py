import unittest
import os

import numpy as np
import pytest

from cppe import CppeState, get_polarizable_sites

from .cache import cache
from . import symmetry_helper as symmetry


class TestFields(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = "{}/potfiles/pna_6w.pot".format(dirname)

    def test_nuclear_fields(self):
        mol = cache.molecule["pna"]
        options = {"potfile": self.potfile_path}
        cppe_state = CppeState(options, mol, lambda s: None)
        cppe_state.calculate_static_energies_and_fields()

        nuclear_fields = cppe_state.nuclear_fields

        potentials = cppe_state.potentials
        ref = np.zeros((len(potentials), 3))
        for i, p in enumerate(potentials):
            for a in mol:
                diff = p.position - a.position
                ref[i, :] += a.charge * diff / np.linalg.norm(diff)**3
        np.testing.assert_allclose(ref.flatten(), nuclear_fields, atol=1e-14)

    @pytest.mark.needs_polarizationsolver
    def base_multipole_fields(self, options):
        from polarizationsolver import fields as polfields

        mol = cache.molecule["pna"]
        cppe_state = CppeState(options, mol, lambda s: None)
        cppe_state.calculate_static_energies_and_fields()

        multipole_fields = cppe_state.multipole_fields

        potentials = cppe_state.potentials
        polsites = get_polarizable_sites(potentials)

        ref = np.zeros((len(potentials), 3))
        for i, p1 in enumerate(polsites):
            for j, p2 in enumerate(polsites):
                if not p1.is_polarizable:
                    continue
                if p1.index == p2.index:
                    continue
                if p1.excludes_site(p2.index):
                    continue
                diff = p1.position - p2.position
                for m in p2.multipoles:
                    M_full = symmetry.unfold_tensor(m.values, m.k)
                    if options.get("damp_multipole", False):
                        a_i = p1.polarizability.isotropic_value
                        a_j = p2.polarizability.isotropic_value
                        damp = cppe_state.options["damping_factor_multipole"]
                        a = 1/(a_i*a_j)**(1/6) * damp
                        ref[p1.index, :] += polfields.thole_exp_field(
                            diff, m.k, M_full, 1, a)
                    else:
                        ref[p1.index, :] += polfields.field(
                            diff, m.k, M_full, 1)

        np.testing.assert_allclose(ref.flatten(), multipole_fields, atol=1e-14)

    @pytest.mark.needs_polarizationsolver
    def test_multipole_fields(self):
        options = {"potfile": self.potfile_path}
        self.base_multipole_fields(options)

    @pytest.mark.needs_polarizationsolver
    def test_multipole_fields_damped(self):
        options = {"potfile": self.potfile_path, "damp_multipole": True}
        self.base_multipole_fields(options)
