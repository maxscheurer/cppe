import unittest
import os
import h5py

import numpy as np

from cppe import PotfileReader
from cppe import Atom, Molecule
from cppe import PeOptions, CppeState


class TestFunctionality(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = "{}/potfiles/pna_6w.pot".format(dirname)
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
        cppe_state.calculate_static_energies_and_fields()
        en_el_nuc = cppe_state.get_current_energies().get("Electrostatic/Nuclear")
        ref = -0.321349401430  # pelib
        np.testing.assert_almost_equal(en_el_nuc, ref, decimal=9)
