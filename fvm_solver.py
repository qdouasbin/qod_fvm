"""
Solving a toy problem with a 1D finite volume solver:

Assumptions:

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from celluloid import Camera

import input_1d_solver as inp

from cell import Cell, N_TRANSPORT_EQ
from field import Field
import FluxReconstruction as FR
import utils
from glob import glob
from np_to_xmf import create_time_collection_xmf
import os

utils.plt_style()

CHECK_1D_DOMAIN = 0
CHECK_1D_INIT = 0
CHECK_1D_BC = 0


def create_oneD_domain():
    # Domain
    dx = (inp.x_max - inp.x_min) / inp.n_cell

    x_f_0 = inp.x_min - dx * inp.stencil_left
    x_f_N = inp.x_max + dx * inp.stencil_right
    n_faces = inp.n_cell + 1 + inp.stencil_left + inp.stencil_left
    n_cells = n_faces - 1

    x_faces = np.linspace(x_f_0, x_f_N, n_faces)

    x_cell_center = x_faces[0:-1] + 0.5 * dx

    if CHECK_1D_DOMAIN:
        y = np.ones_like(x_faces)

        ycc = np.ones_like(x_cell_center) * 1

        plt.figure()
        plt.plot(x_faces, y, ls='none', marker='+', label='faces')
        plt.plot(x_cell_center, ycc, ls='none', marker='.', label='cells')
        plt.legend()
        plt.xlabel(r"x [m]")
        plt.ylabel(r"Arbitrary [-]")
        plt.savefig("Figures/checks/cells.png")
        plt.show()

    lst_cells = []
    for idx_cell in range(n_cells):
        tmp_cell = Cell()
        tmp_cell.set_positions(x_faces[idx_cell],
                               x_cell_center[idx_cell],
                               x_faces[idx_cell + 1])
        lst_cells.append(tmp_cell)

    return Field(lst_cells)


def initialize_field(dom_1D, case=1):
    def area_x(cell):
        x = cell.x_i
        L = 0.95
        d_in = 1e-1
        d_out = 2 * d_in
        x_center_slope = 0.5
        x_end_slope = x_center_slope + L / 2.
        x_beg_slope = x_center_slope - L / 2.

        slope = (d_out - d_in) / (x_end_slope - x_beg_slope)

        if x < x_center_slope - L / 2.:
            diam = d_in
        elif x > x_center_slope + L / 2:
            diam = d_out
        else:
            diam = d_in
            diam += (x - x_beg_slope) * slope

        return diam, 0.25 * np.pi * diam ** 2.

    def area_constant(cell, diam):
        return diam, 0.25 * np.pi * diam ** 2.

    for _idx, _cell in enumerate(dom_1D.lst_cell):

        _cell.gamma = inp.GAMMA
        _cell.r_gas = inp.R_GAS

        if case == 1:
            _cell.set_T(inp.init_T)
            _cell.set_P(inp.init_P)
            _cell.set_u(inp.init_u)
            _cell.set_rho_from_TP()
            _cell.diam, _cell.area = area_x(_cell)
            # _cell.diam, _cell.area = area_constant(_cell, 1e-1)
            _cell.compute_volume()

    return dom_1D


def plot_prim_var_field(dom_1D, fig_name):
    xi = dom_1D.get_xi()
    area = dom_1D.get_area()
    area_max = area.max()
    temp = dom_1D.get_T()
    temp_max = temp.max()
    pres = dom_1D.get_P()
    pres_max = pres.max()
    rho = dom_1D.get_rho()
    rho_max = rho.max()

    u = dom_1D.get_u()
    u_max = np.amax(u)

    plt.figure()
    plt.plot(xi, area / area_max, ':.', label='A/%2.2e' % area_max)
    plt.plot(xi, temp / temp_max, '--+', label='T/%2.2e' % temp_max)
    plt.plot(xi, pres / pres_max, '--x', label='P/%2.2e' % pres_max)
    plt.plot(xi, rho / rho_max, '--', label=r'$\rho$/%2.2e' % rho_max)

    if (u_max):
        plt.plot(xi, u / u_max, '--', label=r'u/%2.2e' % u_max)
    else:
        plt.plot(xi, u, '--', label=r'u [m/s]')

    plt.legend()
    plt.xlabel(r'x [m]')
    plt.savefig("Figures/checks/%s.png" % fig_name)
    plt.show()


def apply_BC(field):
    """
    Apply boundary conditions

           ghost  BC  domain
        |-- q_0 --|-- q_1 --|
                  |<-- q_target

        q_t = 0.5 * (q_0 + q_1)
        q_0 = 2 q_t - q_1
    :param field:
    :return:
    """
    left_ghost = field.lst_cell[0]
    right_cell = field.lst_cell[1]
    left_ghost.cons_to_prim()
    right_cell.cons_to_prim()
    left_ghost.set_u(inp.bc_left_u)
    left_ghost.set_P(right_cell.get_P())
    left_ghost.set_T(inp.bc_left_T )
    left_ghost.set_rho_from_TP()

    right_ghost = field.lst_cell[-1]
    left_of_right_ghost = field.lst_cell[-2]
    right_ghost.set_P(2. * inp.bc_right_P - left_ghost.get_P())
    right_ghost.set_u(left_of_right_ghost.get_u())
    right_ghost.set_T(left_of_right_ghost.get_T())
    right_ghost.set_rho_from_TP()

    for _ghost in [left_ghost, right_ghost]:
        _ghost.update_vec_from_var()

    return field


def time_marching(field, dt):
    """
    For now 1st order forward euler in time
    """
    # cons2prim
    field.update_var_from_vec()

    # apply BC
    apply_BC(field)

    # prim2cons
    field.update_vec_from_var()

    # Compute source terms

    # ----- Source term -----
    # Using the fluxes instead

    # ----- Advection -----

    # first face to be reconstructed == left bc, last = right bc
    # from (0, 1) cell to (Ng-1, Ng) where Ng = n_cell + n_ghost_cells (here 2)
    stencil_cv_it = enumerate(zip(field.lst_cell[0:-1],
                                  field.lst_cell[1:]))

    # reconstruct inter_cell fluxes
    for _idx, (cell_l, cell_r) in stencil_cv_it:
        inter_cell_flux = FR.get_intercell_flux(cell_l, cell_r, dt)
        cell_l.flux_face_r = inter_cell_flux * cell_l.normal_r
        cell_r.flux_face_l = inter_cell_flux * cell_r.normal_l

    # compute residual and update conservative

    # here we go for cell inside the domain only (1, 2) --> (N-2, N-1)
    stencil_cv_it = enumerate(zip(field.lst_cell[1:-1], field.lst_cell[2:]))
    for _idx, (cell_l, cell_r) in stencil_cv_it:
        # _area = 0.5 * (cell_l.area + cell_r.area)
        _area = min(cell_l.area, cell_r.area)

        # flux_adv_diff = cell_r.flux_face_l * cell_r.area
        flux_adv_diff = cell_r.flux_face_l * _area

        # flux_adv_diff -= cell_l.flux_face_l * cell_l.area
        flux_adv_diff -= cell_l.flux_face_l * _area
        source_terms = cell_l.s_cons

        residual = flux_adv_diff

        # Check for nan
        for idx in range(N_TRANSPORT_EQ):
            assert (residual[idx] == residual[idx])

        cell_l.w_cons -= (dt / cell_l.dx) * residual + dt * source_terms

    # update primitive
    field.update_var_from_vec()

    return field


def advance_time_step(step, time, field):
    dt = field.compute_time_step(inp.CFL, step)

    _ = field.get_sos()

    print(" > Step %06d\ttime = %3.3e\tdt = %3.3e" % (step, time, dt))

    field = time_marching(field, dt)

    time = time + dt

    return time, field


if __name__ == "__main__":
    print("Start")

    utils.clean_directories(inp)

    # Create domain
    domain = create_oneD_domain()

    # Init fields
    field = initialize_field(domain)

    if CHECK_1D_INIT:
        plot_prim_var_field(field, 'init_fields')

    field.update_vec_from_var()

    if CHECK_1D_BC:
        plot_prim_var_field(field, 'apply_BC')

    time = inp.t_init
    field.write_sol(0, time)

    for step in range(inp.n_steps):
        step += 1
        time, field = advance_time_step(step, time, field)

        if not step % inp.output_freq:
            field.write_sol(step, time)
            # field.write_output(step, time)
            # field.plot_cons()

    print("Done")
