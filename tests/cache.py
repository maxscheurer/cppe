import os

import h5py

from .misc import cached_property

import cppe


def fullfile(fn):
    thisdir = os.path.dirname(__file__)
    if os.path.isfile(os.path.join(thisdir, fn)):
        return os.path.join(thisdir, fn)
    elif os.path.isfile(fn):
        return fn
    else:
        return ""


class TestdataCache():
    @property
    def testcases(self):
        """
        The definition of the test cases: Data generator and reference file
        """
        cases = ["pna", "nilered"]
        return [k for k in cases
                if os.path.isfile(fullfile(k + ".hdf5"))]

    @cached_property
    def molecule(self):
        """
        The HF data a testcase is based upon
        """
        ret = {}
        for k in self.testcases:
            datafile = fullfile(k + ".hdf5")
            f = h5py.File(datafile, 'r')
            mol = cppe.Molecule()
            for z, coord in zip(f['atom_charges'], f['atom_coords']):
                mol.append(cppe.Atom(z, *coord))
            ret[k] = mol
        return ret


# Setup cache object
cache = TestdataCache()
