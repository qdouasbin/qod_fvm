"""

"""
import numpy as np

import utils

utils.plt_style()

N_TRANSPORT_EQ = 3


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

        self.normal_l = 1.
        self.normal_r = -1

        # variables
        self.pres = 0.
        self.temp = 0.
        self.e_int = 0.
        self.e_tot = 0.
        self.u = 0.

        # self.cp = 0.
        self.gamma = 0.
        self.r_gas = 0.

        # cons
        self.rho = 0.
        self.rho_u = 0.
        self.rho_E = 0.

        # index of vector
        self.idx_mass = 0
        self.idx_momentum = 1
        self.idx_energy = 2

        # cons_vec
        self.w_cons = np.zeros(N_TRANSPORT_EQ)

        # source term
        self.s_cons = np.zeros(N_TRANSPORT_EQ)

        # Flux vector
        # cell-center
        self.f_cons = np.zeros(N_TRANSPORT_EQ)
        # left
        self.flux_face_l = np.zeros(N_TRANSPORT_EQ)
        # right
        self.flux_face_r = np.zeros(N_TRANSPORT_EQ)

        self.n_transport_eq = N_TRANSPORT_EQ

    def compute_volume(self):
        self.vol = self.area * self.dx

    def set_positions(self, x_minus_half, x_center, x_plus_half):
        self.x_i = x_center
        self.x_f_0 = x_minus_half
        self.x_f_1 = x_plus_half
        self.dx = self.x_f_1 - self.x_f_0

    def get_T(self):
        return self.temp

    def get_P(self):
        return self.pres

    def get_u(self):
        return self.u

    def set_T(self, T):
        self.temp = T

    def set_P(self, P):
        self.pres = P

    def set_u(self, u):
        self.u = u

    def set_rho_from_TP(self):
        self.rho = self.pres / (self.r_gas * self.temp)

    def set_T_from_RP(self):
        self.temp = self.pres / (self.r_gas * self.rho)

    def get_cp(self):
        return self.gamma * self.get_cv()

    def get_cv(self):
        return self.r_gas / (self.gamma - 1.)

    def get_enthalpy(self):
        _h = self.get_cp() * self.temp
        return _h

    def get_internal_energy(self):
        _h = self.get_enthalpy()
        self.e_int = _h - self.pres / self.rho
        return self.e_int

    def get_total_energy(self):
        self.e_tot = self.get_internal_energy() + 0.5 * self.u ** 2
        return self.e_tot

    def get_sos(self):
        sos = np.sqrt(self.gamma * self.r_gas * self.temp)
        assert (sos > 0)
        return sos

    def update_cons_vec(self):
        self.w_cons[self.idx_mass] = self.rho
        self.w_cons[self.idx_momentum] = self.rho_u
        self.w_cons[self.idx_energy] = self.rho_E
        self.w_cons *= self.area

    def update_flux_vec(self):
        self.f_cons[self.idx_mass] = self.rho * self.u
        self.f_cons[self.idx_momentum] = self.rho * self.u ** 2. + self.pres
        self.f_cons[self.idx_energy] = self.u * (self.rho * self.e_tot + self.pres)
        self.f_cons *= self.area

    def prim_to_cons(self):
        # Mass
        # no need, rho is prim and cons

        # Momentum
        self.rho_u = self.rho * self.u

        # Energy
        self.rho_E = self.rho * self.get_total_energy()

    def cons_to_prim(self):
        # mass, ok
        self.rho = self.w_cons[0] / self.area

        # momentum
        self.u = self.w_cons[1] / self.rho / self.area

        # energy
        self.e_tot = self.w_cons[2] / (self.area * self.rho)

    def update_vec_from_var(self):
        self.prim_to_cons()
        self.update_cons_vec()
        self.update_flux_vec()

    def update_var_from_vec(self):
        self.cons_to_prim()

        self.pres = (self.gamma - 1.) * self.rho * \
                    (self.e_tot - 0.5 * self.u ** 2)

        # set T
        self.set_T_from_RP()
