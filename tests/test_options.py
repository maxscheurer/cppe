import unittest
from cppe import PeOptions


class TestOptions(unittest.TestCase):
    def test_new_options(self):
        options = PeOptions()
        assert options.options["potfile"] == "potential.pot"
        options.options["potfile"] = "test.pot"
        assert options.options["potfile"] == "test.pot"
