import unittest
import os
import h5py

import numpy as np

from cppe import PotfileReader
from cppe import Atom, Molecule
from cppe import PeOptions, CppeState


def solve_induced_moments(options, mol):
    cppe_state = CppeState(options, mol)
    cppe_state.calculate_static_energies_and_fields()

    static_fields = np.array(cppe_state.get_static_fields())
    zeros = np.zeros_like(static_fields)
    cppe_state.update_induced_moments(zeros, False)
    return cppe_state


class TestFunctionality(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = "{}/potfiles/pna_6w.pot".format(dirname)
    potfile_iso_path = "{}/potfiles/pna_6w_isopol.pot".format(dirname)
    pna_path = "{}/pna.hdf5".format(dirname)

    def test_read_potfile(self):
        p = PotfileReader(self.potfile_path)
        potentials = p.read()
        assert len(potentials) == 18

    def test_compute_nuclear_interaction_energy(self):
        f = h5py.File(self.pna_path, 'r')
        mol = Molecule()
        options = PeOptions()
        options.potfile = self.potfile_path
        for z, coord in zip(f['atom_charges'], f['atom_coords']):
            mol.append(Atom(z, *coord))
        cppe_state = CppeState(options, mol)
        assert cppe_state.get_polarizable_site_number() == 18
        cppe_state.calculate_static_energies_and_fields()
        en_el_nuc = cppe_state.energies["Electrostatic"]["Nuclear"]
        ref = -0.321349401430  # pelib
        np.testing.assert_almost_equal(en_el_nuc, ref, decimal=9)
        # test writing to the energy container from Python
        bla = cppe_state.energies["Electrostatic"]
        bla["Nuclear"] = -10.0
        assert bla["Nuclear"] == -10.0
        cppe_state.energies["Electrostatic"]["Nuclear"] = -20.0
        assert cppe_state.energies["Electrostatic"]["Nuclear"] == -20.0

    def test_iso_pol(self):
        # use iso_pol option
        f = h5py.File(self.pna_path, 'r')
        mol = Molecule()
        options = PeOptions()
        options.potfile = self.potfile_path
        options.iso_pol = True
        for z, coord in zip(f['atom_charges'], f['atom_coords']):
            mol.append(Atom(z, *coord))

        induced_test_state = solve_induced_moments(options, mol)
        induced_test_moments = induced_test_state.get_induced_moments()

        options = PeOptions()
        options.potfile = self.potfile_iso_path
        induced_ref_state = solve_induced_moments(options, mol)
        induced_ref_moments = induced_ref_state.get_induced_moments()

        test_pots = induced_test_state.get_potentials()
        ref_pots = induced_ref_state.get_potentials()

        for pt, pr in zip(test_pots, ref_pots):
            for polt, polr in zip(pt.polarizabilities,
                                  pr.polarizabilities):
                np.testing.assert_almost_equal(polt.values,
                                               polr.values)

        np.testing.assert_almost_equal(induced_test_moments,
                                       induced_ref_moments, decimal=8)
