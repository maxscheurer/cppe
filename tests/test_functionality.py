import unittest
import os

import numpy as np

from .cache import cache

from cppe import PotfileReader
from cppe import PeOptions, CppeState


def print_callback(output):
    print("cb: ", output)


def solve_induced_moments(options, mol):
    cppe_state = CppeState(options, mol, print_callback)
    cppe_state.calculate_static_energies_and_fields()

    static_fields = cppe_state.static_fields
    np.testing.assert_allclose(
        static_fields, cppe_state.nuclear_fields + cppe_state.multipole_fields
    )
    zeros = np.zeros_like(static_fields)
    cppe_state.update_induced_moments(zeros, False)
    return cppe_state


class TestFunctionality(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = "{}/potfiles/pna_6w.pot".format(dirname)
    potfile_iso_path = "{}/potfiles/pna_6w_isopol.pot".format(dirname)

    def test_read_potfile(self):
        p = PotfileReader(self.potfile_path)
        potentials = p.read()
        assert len(potentials) == 18

    def test_compute_nuclear_interaction_energy(self):
        mol = cache.molecule["pna"]
        options = PeOptions()
        options.potfile = self.potfile_path
        cppe_state = CppeState(options, mol, print_callback)
        assert cppe_state.get_polarizable_site_number() == 18
        cppe_state.calculate_static_energies_and_fields()
        en_el_nuc = cppe_state.energies["Electrostatic"]["Nuclear"]
        ref = -0.321349401430  # pelib
        np.testing.assert_almost_equal(en_el_nuc, ref, decimal=9)
        # TODO: split these tests up, use cache
        # test writing to the energy container from Python
        bla = cppe_state.energies["Electrostatic"]
        bla["Nuclear"] = -10.0
        assert bla["Nuclear"] == -10.0
        cppe_state.energies["Electrostatic"]["Nuclear"] = -20.0
        assert cppe_state.energies["Electrostatic"]["Nuclear"] == -20.0

        cppe_state.energies["Electrostatic"]["Electronic"] = 1.0
        cppe_state.energies["Electrostatic"]["Nuclear"] = 2.0
        cppe_state.energies["Electrostatic"]["Multipoles"] = 3.0
        cppe_state.energies["Polarization"]["Electronic"] = 4.0
        cppe_state.energies["Polarization"]["Nuclear"] = 5.0
        cppe_state.energies["Polarization"]["Multipole"] = 6.0
        total_energy = (
            cppe_state.energies["Electrostatic"]["Electronic"]
            + cppe_state.energies["Electrostatic"]["Nuclear"]
            + cppe_state.energies["Electrostatic"]["Multipoles"]
            + cppe_state.energies["Polarization"]["Electronic"]
            + cppe_state.energies["Polarization"]["Nuclear"]
            + cppe_state.energies["Polarization"]["Multipole"]
        )
        assert cppe_state.total_energy == total_energy

    def test_iso_pol(self):
        # use iso_pol option
        mol = cache.molecule["pna"]
        options = PeOptions()
        options.potfile = self.potfile_path
        options.iso_pol = True
        induced_test_state = solve_induced_moments(options, mol)
        induced_test_moments = induced_test_state.get_induced_moments()

        options = PeOptions()
        options.potfile = self.potfile_iso_path
        induced_ref_state = solve_induced_moments(options, mol)
        induced_ref_moments = induced_ref_state.get_induced_moments()

        test_pots = induced_test_state.potentials
        ref_pots = induced_ref_state.potentials

        for pt, pr in zip(test_pots, ref_pots):
            np.testing.assert_almost_equal(pt.polarizability.values,
                                           pr.polarizability.values)
        np.testing.assert_almost_equal(induced_test_moments,
                                       induced_ref_moments, decimal=8)
