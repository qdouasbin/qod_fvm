# -*- coding: utf-8 -*-

import unittest
import pickle
import toml
import glob
from qod_fvm import solver
import numpy as np


class TestSteadyStateSolution(unittest.TestCase):
    """Advanced test cases."""

    # def test_example(self):
    #     # do assert here
    #     pass

    def test_check_steady_state_cst_area(self):
        """
        Constant area, 3 cells, final solution should be uniform
        """
        print(" > test steady state with cst area")
        input_file = toml.load('check_steady_solution/steady_state_cst_area.toml')
        params, field_cst = solver.main(input_file)

        _arr = np.ones_like(field_cst.get_T())

        left_bc = input_file['BoundaryConditions']['left']
        right_bc = input_file['BoundaryConditions']['right']

        assert (left_bc['P'] == right_bc['P'])

        _u_target = left_bc['u']
        _T_target = left_bc['T']
        _P_target = right_bc['P']

        assert (np.allclose(field_cst.get_u(), _u_target * np.ones_like(_arr)))
        assert (np.allclose(field_cst.get_T(), _T_target * np.ones_like(_arr)))
        assert (np.allclose(field_cst.get_P(), _P_target * np.ones_like(_arr)))

    def test_check_steady_state_change_area(self):
        """
        Change area and compare to reference solution
        """
        MAKE_REFERENCE = False

        print(" > test steady state with change area")
        input_file = toml.load('check_steady_solution/steady_state_change_area.toml')
        params, field_sol = solver.main(input_file)

        if MAKE_REFERENCE:
            with open(r"check_steady_solution/solution_change_area.pickle", "wb") as _file:
                pickle.dump(field_sol, _file)
        else:
            with open(r"check_steady_solution/solution_change_area.pickle", "rb") as _file:
                field_ref = pickle.load(_file)

            assert(np.allclose(field_sol.get_u(), field_ref.get_u()))
            assert(np.allclose(field_sol.get_T(), field_ref.get_T()))
            assert(np.allclose(field_sol.get_P(), field_ref.get_P()))
            assert(np.allclose(field_sol.get_sos(), field_ref.get_sos()))


if __name__ == '__main__':
    unittest.main()
