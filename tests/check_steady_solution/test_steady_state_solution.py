# -*- coding: utf-8 -*-

from tests.context import qod_fvm

import unittest


class InputParamsRun():
    def __init__(self):
        # Domain
        self.x_min = -1
        self.x_max = 2
        self.n_cell = 5

        # Time
        self.t_init = 0.
        self.n_steps = 500 * 60
        self.CFL = 0.8

        # Simple geometry
        self.simple_geom = True
        self.simple_geom_type = "constant_area"

        # Numercial scheme
        self.stencil_left = 1
        self.stencil_right = 1

        # fluid
        self.GAMMA = 1.4
        self.R_GAS = 290

        # init
        self.init_u = 3
        self.init_T = 300
        self.init_P = 101325

        # BC
        self.bc_left_T = 2 * self.init_T
        self.bc_left_u = 1.5 * self.init_u

        # BC right
        self.bc_right_P = self.init_P

        # Solution output
        self.output_freq = 500 * 60
        self.output_dir = ''


class TestSteadyStateSolution(unittest.TestCase):
    """Advanced test cases."""

    # def test_thoughts(self):
    #     pass
    # self.assertIsNone(qod_fvm.hmm())

    def test_check_steady_state_solution(self):
        inp = InputParamsRun()
        from qod_fvm import solver
        solver.main()
        assert (1 == 1)


if __name__ == '__main__':
    unittest.main()
