"""
Solving a toy problem with a 1D finite volume solver:

Assumptions:

    - uniform grid spacing (simplified WENO scheme)
    -
"""

import numpy as np
import matplotlib.pyplot as plt

import input_1d_solver as inp

from cell import Cell
from field import Field
import utils

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

    for _idx, _cell in enumerate(dom_1D.lst_cell):

        _cell.gamma = inp.GAMMA
        _cell.r_gas = inp.R_GAS

        if case == 1:
            _cell.set_T(300)
            _cell.set_P(1e5)
            _cell.set_u(0)
            _cell.set_rho_from_TP()
            _cell.diam, _cell.area = area_x(_cell)

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
    field.lst_cell[0].set_u(inp.bc_left_u)
    field.lst_cell[0].set_P(inp.bc_left_P)
    field.lst_cell[0].set_T(inp.bc_left_T)
    field.lst_cell[0].set_rho_from_TP()

    field.lst_cell[-1].set_P(inp.bc_right_P)

    return field


def compute_time_step_size():
    sos = field
    pass


def time_marching(field, dt):
    return field


def advance_time_step(time, field):
    dt = field.compute_time_step(inp.CFL)

    field = time_marching(field, dt)

    time = time + dt

    return time, field


def write_output(step):
    pass


if __name__ == "__main__":
    print("Start")

    # Create domain
    domain = create_oneD_domain()

    # Init fields
    field = initialize_field(domain)

    if CHECK_1D_INIT:
        plot_prim_var_field(field, 'init_fields')

    # Apply BC --> no need to do that here
    field = apply_BC(field)

    if CHECK_1D_BC:
        plot_prim_var_field(field, 'apply_BC')

    field.prim_to_cons()

    time = inp.t_init

    for step in range(inp.n_steps):
        print("\t Step %5d" % step)
        time, field = advance_time_step(time, field)

        if not step % inp.output_freq:
            write_output(step)

    print(field.lst_cell[0].rho)
    print(field.lst_cell[0].rho_u)
    print(field.lst_cell[0].rho_e)
    print(field.lst_cell[0].w_cons)

    print("Done")
