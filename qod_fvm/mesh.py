import numpy as np
import matplotlib.pyplot as plt

from qod_fvm.utils import plt_style
from qod_fvm.cell import Cell
from qod_fvm.field import Field

plt_style()

CHECK_1D_DOMAIN = 0


def create_oneD_domain(params_mesh,
                       params_Numerics,
                       check_domain=CHECK_1D_DOMAIN):
    """
    Create a 1D uniform mesh from the input file.
    The

    :param params_mesh: Dictionary, from TOML input file
    :param params_Numerics: Dictionary, from TOML input file
    :return: Field object
    """
    _x_max = params_mesh['x_max']
    _x_min = params_mesh['x_min']
    _n_cells = params_mesh['n_cells']
    _stencil = params_Numerics['stencil']

    # Domain
    dx = (_x_max - _x_min) / _n_cells

    x_f_0 = _x_min - dx * _stencil
    x_f_N = _x_max + dx * _stencil

    n_faces = _n_cells + 1 + 2 * _stencil
    n_cells = n_faces - 1

    x_faces = np.linspace(x_f_0, x_f_N, n_faces)

    x_cell_center = x_faces[0:-1] + 0.5 * dx

    if check_domain:
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
