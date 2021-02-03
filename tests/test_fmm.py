import os
import numpy as np
import cppe
import pandas as pd

import unittest


class TestFMM(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = f"{dirname}/potfiles/loprop_solvated_20.pot"
    fs_file = f"{dirname}/potfiles/fmuls_loprop_solvated_20.txt"

    def test_field_computation(self):
        potentials = cppe.PotfileReader(self.potfile_path).read()

        options = {}

        if os.path.isfile(self.fs_file):
            fs = np.loadtxt(self.fs_file)
        else:
            fmuls = cppe.MultipoleFields(potentials, options)
            fs = fmuls.compute_tree()
            np.savetxt(self.fs_file, fs)

        ref = pd.read_csv(os.path.join(self.dirname, "ref_fmm_errors.csv"))

        # data = {}
        test_thetas = [
            # 0.2,
            # 0.3,
            0.5,
            0.7,
            0.99,
        ]
        test_orders = [
            3,
            5,
            # 7,
        ]
        for ii, theta in enumerate(test_thetas):
            for exp_order in test_orders:
                options = {
                    "summation_induced_fields": "fmm",
                    "theta": theta,
                    "tree_expansion_order": exp_order,
                }
                fmuls_tree = cppe.MultipoleFields(potentials, options)
                fs_tree = fmuls_tree.compute()
                fs = fs.reshape(len(potentials), 3)
                fs_tree = fs_tree.reshape(fs.shape)
                err_mu_i = (
                    np.linalg.norm(fs - fs_tree, axis=1)
                    / np.linalg.norm(fs, axis=1)
                )
                # data[f"err_F_{theta}_{exp_order}"] = err_mu_i
                np.testing.assert_allclose(
                    ref[f"err_F_{theta}_{exp_order}"], err_mu_i, atol=1e-8
                )
