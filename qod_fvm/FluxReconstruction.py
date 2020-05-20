import numpy as np

FR_TYPE = 'RUSANOV'
# FR_TYPE = 'LAX_FRIEDRICH'


def flux_reconstruction_Rusanov(c_l, c_r):
    # like lax friedrich but with different s+
    u_l_max = np.abs(c_l.u) + c_l.get_sos()
    u_r_max = np.abs(c_r.u) + c_r.get_sos()
    s_plus = max(u_l_max, u_r_max)

    flux_l = c_l.f_cons
    flux_r = c_r.f_cons
    state_l = c_l.w_cons
    state_r = c_r.w_cons

    flux = 0.5 * ((flux_l + flux_r) - s_plus * (state_r - state_l))

    return flux


def flux_reconstruction_LaxFriedrich(flux_l, flux_r, state_l, state_r, dx, dt):
    s_plus = dx / dt

    flux = 0.5 * ((flux_l + flux_r) - s_plus * (state_r - state_l))

    return flux


def get_intercell_flux(cell_l, cell_r, dt):
    dx = 0.5 * (cell_l.dx + cell_r.dx)

    flux_l = cell_l.f_cons
    flux_r = cell_r.f_cons

    state_l = cell_l.w_cons
    state_r = cell_r.w_cons

    if FR_TYPE == 'LAX_FRIEDRICH':
        inter_cell_flux = flux_reconstruction_LaxFriedrich(flux_l,
                                                           flux_r,
                                                           state_l,
                                                           state_r,
                                                           dx,
                                                           dt)
    elif FR_TYPE == 'RUSANOV':

        inter_cell_flux = flux_reconstruction_Rusanov(cell_l, cell_r)

    return inter_cell_flux
