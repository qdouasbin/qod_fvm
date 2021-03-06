"""
Boundary Condition module

The boundaries are initialized using a function and they are added to the field object

Each boundary is a child class of the meta-class "BoundaryCondition".
"""

import numpy as np


class BoundaryCondition:
    """
    Meta class for boundary conditions.
    Initialize directly from each TOML section:
        [BoundaryConditions.my_boundary_condition]

    Each key/item of the dictionary is set as attribute of the class
    """

    def __init__(self, name, dict_bc, verbose=0):
        """constructor meta-BC"""
        self.name_bc = name
        self.idx_cell_bc = None
        self.method_bc = None
        self.verbose = verbose
        for key in dict_bc:
            setattr(self, key + "_bc", dict_bc[key])

    def apply_bc(self):
        """meta-method"""
        raise NotImplementedError("This method should be implemented for the "
                                  "child classes")


class BC_Inlet_UTY(BoundaryCondition):
    """
    Impose UTY at BC
    """

    def __init__(self, name, dict_bc, verbose=0):
        """constructor of BC_INLET_UTY"""
        self.u_bc = None
        self.T_bc = None
        self.Y_bc = None
        if verbose:
            print("init BC_INLET_UTY")
        BoundaryCondition.__init__(self, name, dict_bc, verbose)

    def apply_bc(self, field):
        """
        Apply boundary conditions

        """

        if self.idx_cell_bc == 0:
            _ghost = field.lst_cell[self.idx_cell_bc]
            first_cell = field.lst_cell[self.idx_cell_bc + 1]

            _ghost.set_u(self.u_bc)
            _ghost.set_P(first_cell.get_P())
            _ghost.set_T(self.T_bc)
            _ghost.set_rho_from_TP()

            _ghost.update_vec_from_var()

        else:
            raise NotImplementedError('only left BC for inlet_uty')


class BC_Inlet_UTYP(BoundaryCondition):
    """
    Impose UTYP at inlet (redundant but reduces the acoustic waves)
    """

    def __init__(self, name, dict_bc, verbose=0):
        """constructor of BC_INLET_UTY"""
        self.u_bc = None
        self.T_bc = None
        self.P_bc = None
        self.Y_bc = None
        if verbose:
            print("init BC_INLET_UTYP")
        BoundaryCondition.__init__(self, name, dict_bc, verbose)

    def apply_bc(self, field):
        """
        Apply boundary conditions
        """

        _ghost = field.lst_cell[self.idx_cell_bc]
        first_cell = field.lst_cell[self.idx_cell_bc + 2]

        if self.idx_cell_bc == 0:
            # Compute target value
            if self.method_bc == "Dirichlet":
                _u_target = self.u_bc
                _T_target = self.T_bc
                _p_target = self.P_bc

            elif self.method_bc == "flux":
                _u_target = 2. * self.u_bc - first_cell.get_u()
                _T_target = 2. * self.T_bc - first_cell.get_T()
                _p_target = 2. * self.P_bc - first_cell.get_P()
                _p_target = self.P_bc
            else:
                raise NotImplementedError("BC methods are 'Dirichlet' or 'flux'")
        else:
            raise NotImplementedError('only left BC for inlet_uty')

        _ghost.set_u(_u_target)
        _ghost.rho = first_cell.rho
        _ghost.set_T(_T_target)
        _ghost.set_P(_p_target)

        if self.verbose > 9:
            print("Boundary condition: %s" % str(type(self)))
            print("u, t, p = %e, %e, %e" % (_u_target, _T_target, _p_target))
            print()

        _ghost.set_rho_from_TP()
        _ghost.update_vec_from_var()


class BC_Outlet_P(BoundaryCondition):

    def __init__(self, name, dict_bc, verbose=0):
        """constructor of BC_Outlet_P"""
        self.P_bc = None
        if verbose:
            print("init BC_Outlet_P")
        BoundaryCondition.__init__(self, name, dict_bc, verbose)

    def apply_bc(self, field):
        """
        Apply boundary conditions

        :param field: field object
        """

        if self.idx_cell_bc == -1:
            _ghost = field.lst_cell[self.idx_cell_bc]

            _ghost.update_vec_from_var()
            first_cell = field.lst_cell[self.idx_cell_bc - 1]
            left_left_cell = field.lst_cell[-2]

            _mach = first_cell.get_u() / first_cell.get_sos()

            # Supersonic treatment
            if _mach >= 1.0:
                _ghost.set_P(first_cell.get_u())
                _ghost.set_u(first_cell.get_u())
                _ghost.set_T(first_cell.get_T())
                _ghost.set_rho_from_TP()

            # Subsonic
            else:
                # Compute target value
                if self.method_bc == "Dirichlet":
                    _p_target = self.P_bc
                elif self.method_bc == "flux":
                    _p_target = 2. * self.P_bc - first_cell.get_P()
                else:
                    raise NotImplementedError("BC methods are 'Dirichlet' or 'flux'")

                if first_cell.get_u() > 0.:
                    _ghost.set_P(_p_target)
                    _ghost.set_u(first_cell.get_u())
                    _ghost.set_T(first_cell.get_T())
                else:
                    _ghost.set_P(first_cell.get_P())
                    _ghost.set_u(np.abs(first_cell.get_u()))
                    _ghost.set_T(first_cell.get_T())

                _ghost.set_rho_from_TP()

            if self.verbose > 9:
                print("Boundary condition: %s" % str(type(self)))
                print("p = %e" % (_p_target))
                print()

        else:
            raise NotImplementedError('only right BC for outlet_P')


implemented_bc = ["inlet_UTY", "inlet_UTYP", "outlet_P"]

dict_BCType = {}
dict_BCType['inlet_UTY'] = BC_Inlet_UTY
dict_BCType['inlet_UTYP'] = BC_Inlet_UTYP
dict_BCType['outlet_P'] = BC_Outlet_P


def init_BC(field, params_BC, verbose):
    for bc_name, bc_item in params_BC.items():
        if bc_item["type"] in implemented_bc:
            bc_type = bc_item["type"]
            bc_class = dict_BCType[bc_type]
            bc_init = bc_class(bc_name, bc_item, verbose)
            field.lst_BC.append(bc_init)
