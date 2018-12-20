import unittest
import os

import numpy as np

from cppe import PotfileReader
from cppe import Atom, Molecule
from cppe import PeOptions, CppeState

from pyscf import gto

pna = """
C          8.64800        1.07500       -1.71100
C          9.48200        0.43000       -0.80800
C          9.39600        0.75000        0.53800
C          8.48200        1.71200        0.99500
C          7.65300        2.34500        0.05500
C          7.73200        2.03100       -1.29200
H         10.18300       -0.30900       -1.16400
H         10.04400        0.25200        1.24700
H          6.94200        3.08900        0.38900
H          7.09700        2.51500       -2.01800
N          8.40100        2.02500        2.32500
N          8.73400        0.74100       -3.12900
O          7.98000        1.33100       -3.90100
O          9.55600       -0.11000       -3.46600
H          7.74900        2.71100        2.65200
H          8.99100        1.57500        2.99500
"""


class TestFunctionality(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = "{}/potfiles/pna_6w.pot".format(dirname)

    def test_read_potfile(self):
        p = PotfileReader(self.potfile_path)
        potentials = p.read()
        assert len(potentials) == 18

    def test_compute_nuclear_interaction_energy(self):
        m = gto.Mole(atom=pna, basis="sto3g")
        m.build()
        mol = Molecule()
        options = PeOptions()
        options.potfile = self.potfile_path
        for z, coord in zip(m.atom_charges(), m.atom_coords()):
            mol.append(Atom(z, *coord))
        cppe_state = CppeState(options, mol)
        cppe_state.calculate_static_energies_and_fields()
        en_el_nuc = cppe_state.get_current_energies().get("Electrostatic/Nuclear")
        ref = -0.321349401430  # pelib
        np.testing.assert_almost_equal(en_el_nuc, ref, decimal=9)
