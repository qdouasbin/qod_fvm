# -*- coding: utf-8 -*-

# from .context import qod_fvm

import qod_fvm

import glob
import toml
import unittest


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_absolute_truth_and_meaning(self):
        assert True

    def test_read_toml_example(self):
        input_files = glob.glob("*/*.toml")
        for my_file in input_files:
            try:
                _ = toml.load(my_file)

            # if one input file cannot be read, force wrong assertion
            except:
                print("Issue reading file: %s" % my_file)
                assert (1 == 0)


if __name__ == '__main__':
    unittest.main()
