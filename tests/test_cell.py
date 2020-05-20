# -*- coding: utf-8 -*-

# from .context import qod_fvm

import qod_fvm
from qod_fvm import cell

import glob
import toml
import unittest

import numpy as np


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_cell_1oad(self):
        _cell = cell.Cell()

    def test_cell_vol(self):
        _cell = cell.Cell()
        _dx = np.random.rand()
        _area = np.random.rand()
        _cell.dx = _dx
        _cell.area = _area
        _cell.compute_volume()
        assert (_cell.vol == _dx * _area)

    def test_cell_vol(self):
        _cell = cell.Cell()
        _cell.set_positions(-1.0, 0.0, 1.0)
        assert (_cell.x_i == 0.0)
        assert (_cell.x_f_0 == -1.)
        assert (_cell.x_f_1 == 1.)
        assert (_cell.dx == 2.)

    def test_set_get_T(self):
        _cell = cell.Cell()
        _temp = np.random.rand()
        _cell.set_T(_temp)
        assert(_cell.get_T() == _temp)

    def test_set_get_u(self):
        _cell = cell.Cell()
        _u = np.random.rand()
        _cell.set_u(_u)
        assert(_cell.get_u() == _u)

    def test_set_get_P(self):
        _cell = cell.Cell()
        _pres = np.random.rand()
        _cell.set_P(_pres)
        assert(_cell.get_P() == _pres)

    def test_set_get_rhou(self):
        _cell = cell.Cell()
        _rho = np.random.rand()
        _u = np.random.rand()
        _cell.rho = _rho
        _cell.u = _u
        _cell.prim_to_cons()
        assert(_cell.get_rhou() == _rho * _u)

    def test_set_get_T_from_TP(self):
        _cell = cell.Cell()
        _cell.r_gas = 290.
        _cell.set_T(400.)
        _cell.set_P(101325.0)
        _cell.set_rho_from_TP()

        assert(_cell.rho == 101325 / (290. * 400.))

        _cell = cell.Cell()
        _cell.r_gas = 290.
        _cell.set_T(400.)
        _cell.set_P(101325.0)
        _cell.set_rho_from_TP(400, 101325.)

        assert(_cell.rho == 101325 / (290. * 400.))

    def test_set_get_T_from_RP(self):
        _cell = cell.Cell()
        _cell.r_gas = 290.
        _cell.rho = 1.
        _cell.set_P(101325.0)
        _cell.set_T_from_RP()

        assert(_cell.temp == 101325 / (290. * 1.))

    def test_cp(self):
        _cell = cell.Cell()
        _cell.temp = 300.
        _cell.r_gas = 290.
        _cell.gamma = 1.4
        # cp = r * gamma / (gamma - 1)
        assert(np.allclose(_cell.get_cp(), (1.4 * 290. / (1.4 - 1.0))))

    def test_cv(self):
        _cell = cell.Cell()
        _cell.r_gas = 290.
        _cell.gamma = 1.4
        # cv =  r / (gamma - 1)
        assert(np.allclose(_cell.get_cv(), (290. / (1.4 - 1.0))))

    def test_enthalpy(self):
        _cell = cell.Cell()
        _cell.r_gas = 290.
        _cell.gamma = 1.4
        _cell.temp = 300.
        _cp = _cell.get_cp()
        assert(np.allclose(_cell.get_enthalpy(), _cp * 300.))

    def test_internal_energy(self):
        _cell = cell.Cell()
        _cell.r_gas = 290.
        _cell.gamma = 1.4
        _cell.set_T(300.)
        _cell.set_P(1e5)
        _cell.set_rho_from_TP()
        _cp = _cell.get_cp()
        _e_int = _cell.get_internal_energy()
        _e_int_check = _cp * 300. - 1e5 / _cell.rho
        assert(np.allclose(_e_int, _e_int_check))

    def test_total_energy(self):
        _cell = cell.Cell()
        _cell.r_gas = 290.
        _cell.gamma = 1.4
        _cell.set_T(300.)
        _cell.set_P(1e5)
        _cell.set_rho_from_TP()
        _cell.u = 0.
        assert(np.allclose(_cell.get_total_energy(),
                           _cell.get_internal_energy()))

        _cell.u = np.sqrt(2)
        assert(np.allclose(_cell.get_total_energy(),
                           _cell.get_internal_energy() + 1.))

    def test_sos(self):
        _cell = cell.Cell()
        _cell.r_gas = 290.
        _cell.gamma = 1.4
        _cell.temp = 300.
        assert(_cell.get_sos() == np.sqrt(1.4 * 290. * 300))

    def test_update_cons_vec(self):
        _cell = cell.Cell()
        _tmp_rho = np.random.rand()
        _tmp_rhou = np.random.rand()
        _tmp_rhoE = np.random.rand()
        _cell.rho = _tmp_rho
        _cell.rho_u = _tmp_rhou
        _cell.rho_E = _tmp_rhoE
        _cell.update_cons_vec()
        assert(_cell.w_cons[0] == _tmp_rho)
        assert(_cell.w_cons[1] == _tmp_rhou)
        assert(_cell.w_cons[2] == _tmp_rhoE)

if __name__ == '__main__':
    unittest.main()
