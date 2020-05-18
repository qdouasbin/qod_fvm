import os
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

from qod_fvm.np_to_xmf import NpArray2Xmf
from qod_fvm import utils
from tests.check_steady_solution import input_1d_solver as inp

utils.plt_style()

CHECK_TIME_STEP = 0
TEST_WRITE_FUNCTION = 0
SHOW = 0


class Field:

    def __init__(self, list_of_cell):
        self.lst_cell = list_of_cell
        self.lst_BC = []

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
        return np.array(out)

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
        raise NotImplementedError("This source term should not be used. the source terms are in the fluxes")
        x_pos = self.get_xi()
        pres = self.get_P()
        grad_p = np.gradient(pres, x_pos)

        for _idx_cell, _cell in enumerate(self.lst_cell):
            _cell.s_cons[_cell.idx_momentum] = _cell.area * grad_p[_idx_cell]

    def add_source_term_energy(self):
        raise NotImplementedError("This source term should not be used. the source terms are in the fluxes")
        x_pos = self.get_xi()
        pres = self.get_P()
        area = self.get_area()
        u_vel = self.get_u()

        grad_p_energy = np.gradient(pres * area * u_vel, x_pos)

        for _idx_cell, _cell in enumerate(self.lst_cell):
            _cell.s_cons[_cell.idx_energy] = grad_p_energy[_idx_cell]

    def compute_time_step(self, cfl, step):
        u_arr = self.get_u()
        sos = self.get_sos()
        delta_x = self.get_dx()

        local_dt = cfl * delta_x / np.maximum(np.abs(u_arr + sos), np.abs(u_arr - sos))

        if CHECK_TIME_STEP:
            fig, axes = plt.subplots(2, 2, sharex=True, figsize=(6, 4.5))
            plt.suptitle('Global dt = %e' % np.amin(local_dt))

            axes[0, 0].plot(self.get_xi(), delta_x)
            axes[0, 0].set_ylabel(r"$\Delta x$ [m]")

            axes[0, 1].plot(self.get_xi(), sos)
            axes[0, 1].set_ylabel(r"Speed of sound [m/s]")

            axes[1, 0].plot(self.get_xi(), u_arr)
            axes[1, 0].set_ylabel(r"Velocity [m/s]")

            axes[1, 1].plot(self.get_xi(), local_dt)
            axes[1, 1].set_ylabel(r"local $\Delta t$ [m/s]")

            for idx in range(2):
                axes[1, idx].set_xlabel("x [m]")

            if not step:
                utils.savefig_check('local_dt')
            else:
                utils.savefig_check('local_dt_%6d' % step)

            if SHOW:
                plt.show()

            plt.close()

        dt_min = np.amin(local_dt)
        assert (dt_min == dt_min)

        return dt_min

    def update_vec_from_var(self):
        for _cell in self.lst_cell:
            _cell.update_vec_from_var()

    def update_var_from_vec(self):
        for _cell in self.lst_cell:
            _cell.update_var_from_vec()

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
        # print(w_cons_mat.shape)

        for _idx, _cell in enumerate(self.lst_cell):
            w_cons_mat[_idx, :] = _cell.w_cons

        return w_cons_mat

    def get_flux_cc_matrix(self):
        """cell centered fluxes"""
        n_cells = len(self.lst_cell)
        n_cons = self.lst_cell[0].n_transport_eq
        f_cons_mat = np.zeros((n_cells, n_cons))
        # print(f_cons_mat.shape)

        for _idx, _cell in enumerate(self.lst_cell):
            f_cons_mat[_idx, :] = _cell.f_cons

        return f_cons_mat

    def get_source_terms_matrix(self):
        """cell centered fluxes"""
        n_cells = len(self.lst_cell)
        n_cons = self.lst_cell[0].n_transport_eq
        s_cons_mat = np.zeros((n_cells, n_cons))

        for _idx, _cell in enumerate(self.lst_cell):
            s_cons_mat[_idx, :] = _cell.s_cons

        return s_cons_mat

    def write_output(self, iteration, time, params_IO):
        self.update_var_from_vec()
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

        path = params_IO['directory']
        sol_name = os.path.join(path, "solution_%08d" % iteration)

        print("\t--> Saving %s.csv" % sol_name)
        df = pd.DataFrame.from_dict(dict_out)
        df.to_csv("%s.csv" % sol_name)

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

        utils.savefig_solution('%s.png' % sol_name)

        if SHOW:
            plt.show()

        plt.close()

    def plot_cons(self, params_IO, iteration=None):
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

        fig, axes = plt.subplots(3, 3, sharex=True, figsize=(10, 6))
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

        for idx in range(3):
            axes[2, idx].set_xlabel("x [m]")

        if SHOW:
            plt.show()

        path = params_IO['directory']
        sol_name = os.path.join(path, "cons_solution_%08d" % iteration)
        utils.savefig_solution('%s.png' % sol_name)

        plt.close()

    def write_sol(self, step, time, params_IO):

        if TEST_WRITE_FUNCTION:
            n_pts_x = 100
            n_pts_y = 2

            max_x = 5
            x = np.linspace(0, max_x, n_pts_x)
            y = np.array([0., 0.])

            X, Y = np.meshgrid(x, y)
            f = 4
            Y[0, :] = np.cos(f * x) + 2
            plt.figure()
            plt.plot(x, np.cos(f * x) + 2)
            plt.plot(x, Y[0, :])
            plt.show()
            Y[1, :] = - Y[0, :]
            Z = np.ones((n_pts_y, n_pts_x))

            xmf_out = NpArray2Xmf("./test2D_variableY.h5")
            xmf_out.create_grid(X, Y, Z)

            TEST_U = x * Z

            xmf_out.add_field(TEST_U, "foobar")

            xmf_out.dump()

        x_arr = self.get_xi()
        area_arr = self.get_area()
        radius = 0.5 * np.sqrt(4 * area_arr / np.pi)
        _cell = self.lst_cell[0]

        n_pts_x = len(x_arr)
        n_pts_y = 2

        y = np.array([0., 0.])

        X, Y = np.meshgrid(x_arr, y)
        Z = np.ones((n_pts_y, n_pts_x))

        Y[0, :] = radius
        Y[1, :] = - radius

        path = params_IO['directory']
        sol_name = os.path.join(path, "solut_%08d.h5" % step)

        xmf_out = NpArray2Xmf(sol_name, time=time)
        xmf_out.create_grid(X, Y, Z)

        output_fields = {}
        tmp = self.get_u()
        output_fields['vel-x'] = tmp

        tmp = self.get_mach()
        output_fields['Mach'] = tmp

        tmp = self.get_sos()
        output_fields['SoS'] = tmp

        tmp = self.get_gamma()
        output_fields['gamma'] = tmp

        tmp = self.get_r_gas()
        output_fields['R_GAS'] = tmp

        tmp = self.get_rho_e()
        output_fields['rhoE'] = tmp

        tmp = self.get_rho()
        output_fields['rho'] = tmp

        tmp = self.get_P()
        output_fields['P'] = tmp

        tmp = self.get_T()
        output_fields['T'] = tmp

        tmp = self.get_rho_u()
        output_fields['rhou'] = tmp

        tmp = self.get_dx()
        output_fields['dx'] = tmp

        tmp = np.linspace(0, len(tmp - 1), len(tmp))
        output_fields['idx_cv'] = tmp

        tmp = self.get_area()
        output_fields['A'] = tmp

        fluxes = self.get_flux_cc_matrix()
        flux_mass = fluxes[:, _cell.idx_mass]
        flux_momentum = fluxes[:, _cell.idx_momentum]
        flux_energy = fluxes[:, _cell.idx_energy]
        output_fields['flux_mass'] = flux_mass
        output_fields['flux_momentum'] = flux_momentum
        output_fields['flux_energy'] = flux_energy

        sources = self.get_source_terms_matrix()
        sources_mass = sources[:, _cell.idx_mass]
        sources_momentum = sources[:, _cell.idx_momentum]
        sources_energy = sources[:, _cell.idx_energy]
        output_fields['sources_mass'] = sources_mass
        output_fields['sources_momentum'] = sources_momentum
        output_fields['sources_energy'] = sources_energy

        for key, item in output_fields.items():
            xmf_out.add_field(output_fields[key] * Z, key)

        print("\t--> Writing solution: %s" % sol_name)

        xmf_out.dump()

    def write_sol_adim(self, step, time, params):

        x_arr = self.get_xi()
        area_arr = self.get_area()
        radius = 0.5 * np.sqrt(4 * area_arr / np.pi)
        _cell = self.lst_cell[0]

        n_pts_x = len(x_arr)
        n_pts_y = 2

        y = np.array([0., 0.])

        X, Y = np.meshgrid(x_arr, y)
        Z = np.ones((n_pts_y, n_pts_x))

        Y[0, :] = radius
        Y[1, :] = - radius

        path = params['IO']['directory']
        sol_name = os.path.join(path, "solution_normalized_%08d.h5" % step)

        xmf_out = NpArray2Xmf(sol_name, time=time)
        xmf_out.create_grid(X, Y, Z)

        # normalizing quantities
        _u0 = params["InitialSolution"]["init_u"]
        _T0 = params["InitialSolution"]["init_T"]
        _P0 = params["InitialSolution"]["init_P"]
        _rho0 = (_P0 / _cell.r_gas / _T0)

        output_fields = {}
        tmp = self.get_u()
        output_fields['vel-x'] = tmp / _u0

        tmp = self.get_mach()
        output_fields['Mach'] = tmp
        #
        # tmp = self.get_sos()
        # output_fields['SoS'] = tmp

        # tmp = self.get_gamma()
        # output_fields['gamma'] = tmp
        #
        # tmp = self.get_r_gas()
        # output_fields['R_GAS'] = tmp

        # tmp = self.get_rho_e()
        # output_fields['rhoE'] = tmp
        #

        tmp = self.get_rho()
        output_fields['rho'] = tmp / _rho0

        tmp = self.get_P()
        output_fields['P'] = tmp / _P0

        tmp = self.get_T()
        output_fields['T'] = tmp / _T0

        tmp = self.get_rho_u()
        output_fields['rhou'] = tmp / _T0 / _rho0

        tmp = self.get_dx()
        output_fields['dx'] = tmp

        tmp = np.linspace(0, len(tmp - 1), len(tmp))
        output_fields['idx_cv'] = tmp

        tmp = self.get_area()
        output_fields['A'] = tmp

        fluxes = self.get_flux_cc_matrix()
        flux_mass = fluxes[:, _cell.idx_mass]
        flux_momentum = fluxes[:, _cell.idx_momentum]
        flux_energy = fluxes[:, _cell.idx_energy]
        output_fields['flux_mass'] = flux_mass
        output_fields['flux_momentum'] = flux_momentum
        output_fields['flux_energy'] = flux_energy

        sources = self.get_source_terms_matrix()
        sources_mass = sources[:, _cell.idx_mass]
        sources_momentum = sources[:, _cell.idx_momentum]
        sources_energy = sources[:, _cell.idx_energy]
        output_fields['sources_mass'] = sources_mass
        output_fields['sources_momentum'] = sources_momentum
        output_fields['sources_energy'] = sources_energy

        for key, item in output_fields.items():
            xmf_out.add_field(output_fields[key] * Z, key)

        print("\t--> Writing solution: %s" % sol_name)

        xmf_out.dump()

    def apply_bc(self):
        """
        Apply the boudnary conditions.
        Each item of the list is a boundary condition class an has its own method
        to apply the values at the ghost cells
        """
        for my_bc in self.lst_BC:
            my_bc.apply_bc(self)
