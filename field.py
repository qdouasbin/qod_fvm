import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

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
            out.append(_cell.rho_E)

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
        # _gamma = self.get_gamma()
        # _r = self.get_r_gas()
        # _T = self.get_T()
        # return np.sqrt(_gamma * _r * _T)
        out = []
        for _cell in self.lst_cell:
            out.append(_cell.get_sos())

        return np.array(out)

    def get_mach(self):
        return self.get_u() / self.get_sos()

    def add_source_term_p(self):
        x_pos = self.get_xi()
        pres = self.get_P()
        grad_p = np.gradient(pres, x_pos)

        for _idx_cell, _cell in enumerate(self.lst_cell):
            _cell.s_cons[_cell.idx_momentum] = - _cell.area * grad_p[_idx_cell]


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

    def update_vec_from_var(self):
        for _cell in self.lst_cell:
            _cell.update_vec_from_var()

    def update_var_from_vec(self):
        for _cell in self.lst_cell:
            _cell.update_var_from_vec

    def prim_to_cons(self):
        # rho --> no need (both prim and cons)
        # rhoE
        for _cell in self.lst_cell:
            _cell.prim_to_cons()

    def cons_to_prim(self):
        # rho --> no need (both prim and cons)
        # rhoE
        for _cell in self.lst_cell:
            _cell.cons_to_prim()

    def get_cons_matrix(self):
        n_cells = len(self.lst_cell)
        # n_cons = self.lst_cell[0].w_cons.flatten().shape
        n_cons = self.lst_cell[0].n_transport_eq
        w_cons_mat = np.zeros((n_cells, n_cons))
        print(w_cons_mat.shape)

        for _idx, _cell in enumerate(self.lst_cell):
            w_cons_mat[_idx, :] = _cell.w_cons

        return w_cons_mat

    def get_flux_cc_matrix(self):
        """cell centered fluxes"""
        n_cells = len(self.lst_cell)
        # n_cons = self.lst_cell[0].w_cons.flatten().shape
        n_cons = self.lst_cell[0].n_transport_eq
        f_cons_mat = np.zeros((n_cells, n_cons))
        print(f_cons_mat.shape)

        for _idx, _cell in enumerate(self.lst_cell):
            f_cons_mat[_idx, :] = _cell.f_cons

        return f_cons_mat

    def get_source_terms_matrix(self):
        """cell centered fluxes"""
        n_cells = len(self.lst_cell)
        # n_cons = self.lst_cell[0].w_cons.flatten().shape
        n_cons = self.lst_cell[0].n_transport_eq
        s_cons_mat = np.zeros((n_cells, n_cons))
        print(s_cons_mat.shape)

        for _idx, _cell in enumerate(self.lst_cell):
            s_cons_mat[_idx, :] = _cell.s_cons

        return s_cons_mat

    def write_output(self, iteration, time):
        dict_out = {}
        dict_out['x'] = self.get_xi()
        dict_out['diam'] = self.get_diam()
        dict_out['area'] = self.get_area()
        dict_out['P'] = self.get_P()
        dict_out['T'] = self.get_T()
        dict_out['rho'] = self.get_rho()
        dict_out['u'] = self.get_u()
        dict_out['rhou'] = self.get_rho_u()
        dict_out['rhoe'] = self.get_rho_e()
        dict_out['r_gas'] = self.get_r_gas()
        dict_out['gamma'] = self.get_gamma()
        dict_out['mach'] = self.get_mach()
        dict_out['sos'] = self.get_sos()
        dict_out['iteration'] = iteration * np.ones(len(self.lst_cell))
        dict_out['time'] = time * np.ones(len(self.lst_cell))

        df = pd.DataFrame.from_dict(dict_out)
        df.to_csv("solut/solution_%08d.csv" % iteration)

        fig, axes = plt.subplots(2, 3, sharex=True, figsize=(9, 4.5))
        plt.suptitle('Iteration %s\t' % iteration + r' $t$ = %3.3e' % time)
        axes[0, 0].plot(df.x, df.P, '--.')
        axes[0, 0].set_ylabel("Pressure [Pa]")

        axes[0, 1].plot(df.x, df['T'], '--.')
        axes[0, 1].set_ylabel("Temperature [K]")

        axes[0, 2].plot(df.x, df['sos'], '--.')
        axes[0, 2].set_ylabel("Speed of sound [m/s]")

        axes[1, 0].plot(df.x, df.rho, '--.')
        axes[1, 0].set_ylabel(r"Density [kg/m$^3$]")

        axes[1, 1].plot(df.x, df.u, '--.')
        axes[1, 1].set_ylabel("Velocity [m/s]")

        axes[1, 2].plot(df.x, df.mach, '--.')
        axes[1, 2].set_ylabel("Mach [-]")

        for idx in range(3):
            axes[1, idx].set_xlabel("x [m]")

        utils.savefig_solution('solut_%08d.png' % iteration)
        plt.show()

    def plot_cons(self, iteration=None):
        _cell = self.lst_cell[0]
        x_pos = self.get_xi()

        cons = self.get_cons_matrix()
        cons_mass = cons[:, _cell.idx_mass]
        cons_momentum = cons[:, _cell.idx_momentum]
        cons_energy = cons[:, _cell.idx_energy]

        flux = self.get_flux_cc_matrix()
        flux_mass = flux[:, _cell.idx_mass]
        flux_momentum = flux[:, _cell.idx_momentum]
        flux_energy = flux[:, _cell.idx_energy]

        s_term = self.get_source_terms_matrix()
        s_term_mass = s_term[:, _cell.idx_mass]
        s_term_momentum = s_term[:, _cell.idx_momentum]
        s_term_energy = s_term[:, _cell.idx_energy]

        fig, axes = plt.subplots(3, 3, sharex=True, figsize=(9, 7))
        plt.suptitle('Conservatives, cell-centered fluxes and source terms')
        axes[0, 0].plot(x_pos, cons_mass, '--.')
        axes[0, 0].set_ylabel(r"$\rho$")

        axes[0, 1].plot(x_pos, cons_momentum, '--.')
        axes[0, 1].set_ylabel(r"$\rho u$")

        axes[0, 2].plot(x_pos, cons_energy, '--.')
        axes[0, 2].set_ylabel(r"$\rho E$")
        #
        axes[1, 0].plot(x_pos, flux_mass, '--.')
        axes[1, 0].set_ylabel(r"$\rho u$")

        axes[1, 1].plot(x_pos, flux_momentum, '--.')
        axes[1, 1].set_ylabel(r"$\rho u^2$")

        axes[1, 2].plot(x_pos, flux_energy, '--.')
        axes[1, 2].set_ylabel(r"$\rho u (E + p)$")
        #
        axes[2, 0].plot(x_pos, s_term_mass, '--.')
        axes[2, 0].set_ylabel(r"$S_{\rm mass}$")

        axes[2, 1].plot(x_pos, s_term_momentum, '--.')
        axes[2, 1].set_ylabel(r"$S_{\rm momentum}$")

        axes[2, 2].plot(x_pos, s_term_energy, '--.')
        axes[2, 2].set_ylabel(r"$S_{\rm energy}$")

        plt.tight_layout()
        #

        for idx in range(3):
            axes[2, idx].set_xlabel("x [m]")

        plt.show()
