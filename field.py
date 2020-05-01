import numpy as np
import matplotlib.pyplot as plt

import utils

utils.plt_style()

CHECK_TIME_STEP = 0


class Field():

    def __init__(self, list_of_cell):
        self.lst_cell = list_of_cell

    def get_xi(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.x_i)

        return np.array(out)

    def get_diam(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.diam)

        return np.array(out)

    def get_area(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.area)

        return np.array(out)

    def get_T(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.temp)

        return np.array(out)

    def get_P(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.pres)

        return np.array(out)

    def get_rho(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.rho)

        return np.array(out)

    def get_u(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.u)

        return np.array(out)

    def get_rho_u(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.rho_u)

    def get_rho_e(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.rho_e)

        return np.array(out)

    def get_r_gas(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.r_gas)

        return np.array(out)

    def get_gamma(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.gamma)

        return np.array(out)

    def get_dx(self):
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.dx)

        return np.array(out)

    def get_sos(self):
        """
        Get speed of sound
        :return: sos
        """
        _gamma = self.get_gamma()
        _r = self.get_r_gas()
        _T = self.get_T()
        return np.sqrt(_gamma * _r * _T)

    def compute_time_step(self, cfl):
        u_arr = self.get_u()
        sos = self.get_sos()
        delta_x = self.get_dx()

        local_dt = cfl * delta_x / np.maximum(np.abs(u_arr + sos), np.abs(u_arr - sos))

        if CHECK_TIME_STEP:
            fig = plt.figure()
            plt.plot(self.get_xi(), delta_x)
            plt.xlabel("x [m]")
            plt.ylabel(r"$\Delta x$ [m]")
            plt.legend()
            utils.savefig_check('dx')

            fig = plt.figure()
            plt.plot(self.get_xi(), sos)
            plt.xlabel("x [m]")
            plt.ylabel(r"Speed of sound [m/s]")
            plt.legend()
            utils.savefig_check('sos')

            fig = plt.figure()
            plt.plot(self.get_xi(), u_arr)
            plt.xlabel("x [m]")
            plt.ylabel(r"Velocity [m/s]")
            plt.legend()
            utils.savefig_check('velocity')

            fig = plt.figure()
            plt.plot(self.get_xi(), local_dt)
            plt.xlabel("x [m]")
            plt.ylabel(r"local $\Delta t$ [m/s]")
            plt.title('Global dt = %e' % np.amin(local_dt))
            plt.legend()
            utils.savefig_check('local_dt')

            plt.show()

        return np.amin(local_dt)

    def prim_to_cons(self):
        # rho --> no need (both prim and cons)
        # rhoE
        for _cell in self.lst_cell:
            _cell.prim_to_cons()
