"""
Solving a toy problem with a 1D finite volume solver:

    Details:

        1. Equations: Euler equations
        2. Numerics:
            - 1st order in time (forward Euler)
            - 1st order in space (Rusanov or Lax-Friedrich flux reconstruction)
        3. Boundary conditions:
            - left: inlet_uty
            - right: outlet_p
            - types: Dirichlet or flux formulation
        4. Input file: TOML format
        5. Fluid properties
            - gaseous
            - calorically perfect gas
            - non-viscous

"""

import glob

import numpy as np
import matplotlib.pyplot as plt

from qod_fvm import FluxReconstruction as FR, utils
from qod_fvm import bc
from qod_fvm import mesh

utils.plt_style()

CHECK_1D_INIT = 0


def initialize_field(dom_1D, params_init, params_geom, params_fluid):
    """
    Initialize fields. Both fluid properties and geometries
    :param dom_1D: a Field object to be intialize (the cells inside)
    :param params_init: input params for init (section of TOML file)
    :return: initialize domain (Field onject)
    """

    _implemented_geom = ["constant_area", "change_area"]

    if not params_geom["simple_geom"]:
        raise NotImplementedError("Only 'simple_geom' option is supported")

    def area_x(cell):
        """
        Create an area change
        :param cell:
        :return: diameter array, area array
        """
        x_loc = cell.x_i
        length_area_change = 0.95
        d_in = 1e-1
        d_out = 2 * d_in
        x_center_slope = 0.5
        x_end_slope = x_center_slope + length_area_change / 2.
        x_beg_slope = x_center_slope - length_area_change / 2.

        slope = (d_out - d_in) / (x_end_slope - x_beg_slope)

        if x_loc < x_center_slope - length_area_change / 2.:
            diam = d_in
        elif x_loc > x_center_slope + length_area_change / 2:
            diam = d_out
        else:
            diam = d_in
            diam += (x_loc - x_beg_slope) * slope

        return diam, 0.25 * np.pi * diam ** 2.

    def area_constant(cell, diam):
        return diam, 0.25 * np.pi * diam ** 2.


    if params_geom["simple_geom_type"] == "constant_area":
        for _idx, _cell in enumerate(dom_1D.lst_cell):
            _cell.diam, _cell.area = area_constant(_cell, 1e-1)
            _cell.compute_volume()

    elif params_geom["simple_geom_type"] == "change_area":
        for _idx, _cell in enumerate(dom_1D.lst_cell):
            _cell.diam, _cell.area = area_x(_cell)
            _cell.compute_volume()
    else:
        str_err = 'Only simple geometry options are: '
        for item in _implemented_geom:
            str_err += "\n\t%s" % item
        raise NotImplementedError(str_err)

    if params_init["init_type"] == "Uniform":
        for _idx, _cell in enumerate(dom_1D.lst_cell):
            _cell.gamma = params_fluid["GAMMA"]
            _cell.r_gas = params_fluid["R_GAS"]

            _cell.set_T(params_init['init_T'])
            _cell.set_P(params_init['init_P'])
            _cell.set_u(params_init['init_u'])
            _cell.set_rho_from_TP()
    else:
        raise NotImplementedError("Only a uniform initial field "
                                  "is possible so far.")

    return dom_1D


def time_marching(field, dt):
    """
    For now 1st order forward euler in time
    """
    # cons2prim
    field.update_var_from_vec()

    # apply BC
    field.apply_bc()

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
        flux_adv_diff = cell_r.flux_face_l * cell_r.area
        flux_adv_diff -= cell_l.flux_face_l * cell_l.area

        source_terms = cell_l.s_cons

        # add Euler fluxes
        cell_l.w_cons -= (dt / cell_l.dx) * flux_adv_diff

        # add source terms
        cell_l.w_cons -= dt * source_terms

    # update primitive
    field.update_var_from_vec()

    return field


def advance_time_step(step, time, CFL, field):
    dt = field.compute_time_step(CFL, step)

    _ = field.get_sos()

    print(" > Step %06d\ttime = %3.3e\tdt = %3.3e" % (step, time, dt))

    field = time_marching(field, dt)

    time = time + dt

    return time, field


def main(params):
    print("Start")

    if "IO" in params.keys():
        utils.clean_directories(params['IO'])

    # Create domain
    domain = mesh.create_oneD_domain(params['Mesh'],
                                     params['Numerics'])

    # Init fields
    field = initialize_field(domain,
                             params["InitialSolution"],
                             params["Geometry"],
                             params["Fluid"])

    if CHECK_1D_INIT:
        utils.plot_prim_var_field(field, 'init_fields')

    bc.init_BC(field, params['BoundaryConditions'])

    field.update_vec_from_var()

    time = params["TimeMarching"]['time_init']
    time_max = params["TimeMarching"]['time_end']

    if "IO" in params.keys():
        field.write_sol(0, time, params['IO'])
        field.write_sol_adim(0, time, params)
        # field.write_output(0, time, params['IO'])
        # field.plot_cons(params['IO'], iteration=0)

    for step in range(params["TimeMarching"]["n_steps_max"]):
        if time > time_max:
            break

        step += 1
        time, field = advance_time_step(step, time,
                                        params['TimeMarching']['CFL'],
                                        field)

        # todo: do something better for the output (different frequencies per output)
        if "IO" in params.keys():
            if not step % params["IO"]["frequency"]:
                field.write_sol(step, time, params['IO'])
                field.write_sol_adim(step, time, params)
                # field.write_output(step, time, params['IO'])
                # field.plot_cons(params['IO'], iteration=step)

    print("Done")

    return params, field


if __name__ == "__main__":
    import toml

    path = glob.glob('*/*.toml')
    print(path)
    input_file = toml.load(path[-1])
    print(input_file.keys())
    _ = main(input_file)
