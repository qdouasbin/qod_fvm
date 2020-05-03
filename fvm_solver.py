"""
Solving a toy problem with a 1D finite volume solver:

Assumptions:

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from celluloid import Camera

import input_1d_solver as inp

from cell import Cell
from field import Field
import FluxReconstruction as FR
import utils

utils.plt_style()

CHECK_1D_DOMAIN = 0
CHECK_1D_INIT = 1
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
        L = 0.4
        d_in = 1e-2
        d_out = 2e-2
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
            # _cell.diam, _cell.area = area_x(_cell)
            _cell.diam, _cell.area = area_constant(_cell, 1e-2)
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
    # todo: cons2prim for ghost cells
    left_ghost = field.lst_cell[0]
    left_ghost.set_u(inp.bc_left_u)
    left_ghost.set_P(inp.bc_left_P)
    left_ghost.set_T(inp.bc_left_T)
    left_ghost.set_rho_from_TP()

    right_ghost = field.lst_cell[-1]
    right_ghost.set_P(inp.bc_right_P)

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
    field.add_source_term_p()

    # ----- Advection -----


    # first face to be reconstructed == left bc, last = right bc
    # from (0, 1) cell to (Ng-1, Ng) where Ng = n_cell + n_ghost_cells (here 2)
    stencil_cv_it = enumerate(zip(field.lst_cell[0:-1],
                                  field.lst_cell[1:]))

    # reconstruct inter_cell fluxes
    for _idx, (cell_l, cell_r) in stencil_cv_it:
        # print("x_face = %e = %e" % (cell_l.x_f_1, cell_r.x_f_0))
        inter_cell_flux = FR.get_intercell_flux(cell_l, cell_r, dt)
        cell_l.flux_face_r = inter_cell_flux
        cell_r.flux_face_l = inter_cell_flux

    # compute residual and update conservative

    # here we go for cell inside the domain only (1, 2) --> (N-2, N-1)
    stencil_cv_it = enumerate(zip(field.lst_cell[1:-1], field.lst_cell[2:]))
    for _idx, (cell_l, cell_r) in stencil_cv_it:
        # print("cell_l.xi = %e, cell_r.xi = %e" % (cell_l.x_i, cell_r.x_i))

        # print("cell_l.flux_face_l \t:%e\t%e\t%e" % (cell_l.flux_face_l[0], cell_l.flux_face_l[1], cell_l.flux_face_l[2]))
        # print("cell_l.flux_face_r \t:%e\t%e\t%e" % (cell_l.flux_face_r[0], cell_l.flux_face_r[1], cell_l.flux_face_r[2]))
        #
        # print("cell_l.area \t = %e" % cell_l.area)
        # print("cell_r.area \t = %e" % cell_r.area)
        # print("cell_r.area \t = %e" % cell_r.area)
        #
        # print("vol = %e" % cell_l.vol)

        flux_adv_diff = cell_r.flux_face_l * cell_r.area

        flux_adv_diff -= cell_l.flux_face_l * cell_l.area
        source_terms = cell_l.s_cons * cell_l.area
        residual = flux_adv_diff - source_terms
        cell_l.w_cons -= (dt / cell_l.vol) * residual


    # update primitive
    field.update_var_from_vec()

    return field


def advance_time_step(step, time, field):
    dt = field.compute_time_step(inp.CFL)

    print(" > Step %06d\ttime = %3.3e\tdt = %3.3e" % (step, time, dt))

    field = time_marching(field, dt)

    time = time + dt

    return time, field



if __name__ == "__main__":
    print("Start")

    # Create domain
    domain = create_oneD_domain()

    # Init fields
    field = initialize_field(domain)

    if CHECK_1D_INIT:
        plot_prim_var_field(field, 'init_fields')

    # # Apply BC --> no need to do that here
    # field = apply_BC(field)

    field.update_vec_from_var()

    # Check fields after BC
    # field.write_output(0,0)
    # field.plot_cons()

    if CHECK_1D_BC:
        plot_prim_var_field(field, 'apply_BC')

    # field.plot_cons()

    time = inp.t_init

    for step in range(inp.n_steps):
        time, field = advance_time_step(step, time, field)

        if not step % inp.output_freq:
            field.write_output(step, time)
            field.plot_cons()

    print(field.lst_cell[0].rho)
    print(field.lst_cell[0].rho_u)
    print(field.lst_cell[0].rho_E)
    print(field.lst_cell[0].w_cons)

    print(field.lst_cell[0].f_cons)

    print("Done")
