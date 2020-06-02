import pytest
import unittest
import os

from cppe import CppeState


class TestOptions(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = "{}/potfiles/pna_6w.pot".format(dirname)
    potfile_iso_path = "{}/potfiles/pna_6w_isopol.pot".format(dirname)

    default_options = {
        "potfile": "potential.pot",
        "iso_pol": False,
        "induced_thresh": 1e-8,
        "maxiter": 50,
        "damp_induced": False,
        "damp_multipole": False,
        "damping_factor_induced": 2.1304,
        "damping_factor_multipole": 2.1304,
    }

    def test_defaults(self):
        options_dict = {"potfile": self.potfile_path}
        cppe_state = CppeState(options_dict)

        defaults = self.default_options.copy()
        defaults.pop("potfile")
        for k in defaults:
            assert defaults[k] == cppe_state.options[k]

    def test_set_all(self):
        custom_options = {
            "potfile": self.potfile_path,
            "iso_pol": True,
            "induced_thresh": 1e-4,
            "maxiter": 500,
            "damp_induced": True,
            "damp_multipole": True,
            "damping_factor_induced": 213.04,
            "damping_factor_multipole": 213.04,
        }
        cppe_state = CppeState(custom_options)
        for k in custom_options:
            assert custom_options[k] == cppe_state.options[k]

    def test_invalid_key(self):
        defaults = self.default_options.copy()
        defaults["invalid_key"] = 42
        with pytest.raises(ValueError):
            CppeState(defaults)
