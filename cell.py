"""

"""
import numpy as np

import utils

utils.plt_style()


class Cell():

    def __init__(self):
        # Cell size [m]
        self.dx = 0.
        # Cell position
        self.x_i = 0.

        # area
        self.diam = 0.
        self.area = 0.

        # Left face position
        self.x_f_0 = 0.

        # right face position
        self.x_f_1 = 0.

        # How does area and volume work?
        self.vol = 0.
        # Do we have trapezoidal elements?

        # variables
        self.pres = 0.
        self.temp = 0.
        self.e_tot = 0.
        self.u = 0.

        # self.cp = 0.
        self.gamma = 0.
        self.r_gas = 0.

        # cons
        self.rho = 0.
        self.rho_u = 0.
        self.rho_e = 0.

        # index of vector
        self._idx_mass = 0
        self._idx_momentum = 1
        self._idx_energy = 2

        # cons_vec
        self.w_cons = np.zeros_like(3)

        # Flux vect (cell-center)
        self.f_cons = np.zeros_like(3)

    def set_positions(self, x_minus_half, x_center, x_plus_half):
        self.x_i = x_center
        self.x_f_0 = x_minus_half
        self.x_f_1 = x_plus_half
        self.dx = self.x_f_1 - self.x_f_0

    def set_T(self, T):
        self.temp = T

    def set_P(self, P):
        self.pres = P

    def set_u(self, u):
        self.u = u

    def set_rho_from_TP(self):
        self.rho = self.pres / (self.r_gas * self.temp)

    def get_cp(self):
        return self.gamma * self.get_cv()

    def get_cv(self):
        return self.r_gas / (self.gamma - 1.)

    def get_enthalpy(self):
        _h = self.get_cp() * self.temp
        return _h

    def get_internal_energy(self):
        _h = self.get_enthalpy()
        _e = _h - self.pres / self.rho
        self.e_tot = _e
        return _e

    def update_cons_vec(self):
        self.w_cons[self._idx_mass] = self.rho
        self.w_cons[self._idx_momentum] = self.rho_u
        self.w_cons[self._idx_energy] = self.rho_e

    def prim_to_cons(self):
        # Mass
        # no need, rho is prim and cons

        # Momentum
        self.rho_u = self.rho * self.u

        # Energy
        _e = self.get_internal_energy()
        self.rho_e = self.rho * _e

        self.update_cons_vec()
