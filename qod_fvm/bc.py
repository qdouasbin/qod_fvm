"""
Boundary Condition module

The boundaries are initialized using a function and they are added to the field object

Each boundary is a child class of the meta-class "BoundaryCondition".
"""


class BoundaryCondition:
    """
    Meta class for boundary conditions.
    Initialize directly from each TOML section:
        [BoundaryConditions.my_boundary_condition]

    Each key/item of the dictionary is set as attribute of the class
    """

    def __init__(self, name, dict_bc):
        """constructor meta-BC"""
        self.name_bc = name
        self.idx_cell_bc = None
        self.method_bc = None
        for key in dict_bc:
            setattr(self, key + "_bc", dict_bc[key])

    def apply_bc(self):
        """meta-method"""
        raise NotImplementedError("This method should be implemented for the "
                                  "child classes")


class BC_Inlet_UTY(BoundaryCondition):

    def __init__(self, name, dict_bc):
        """constructor of BC_INLET_UTY"""
        self.u_bc = None
        self.T_bc = None
        self.Y_bc = None
        print("init BC_INLET_TPY")
        BoundaryCondition.__init__(self, name, dict_bc)

    def apply_bc(self, field):
        """
        Apply boundary conditions
    
            ```
               ghost  BC  domain
            |-- q_0 --|-- q_1 --|
                      |<-- q_target

            q_t = 0.5 * (q_0 + q_1)
            q_0 = 2 q_t - q_1
            ```
        """

        if self.idx_cell_bc == 0:
            _ghost = field.lst_cell[self.idx_cell_bc]
            first_cell = field.lst_cell[self.idx_cell_bc + 1]

            _ghost.set_u(2. * self.u_bc - first_cell.get_u())
            _ghost.set_P(first_cell.get_P())
            _ghost.set_T(2. * self.T_bc - first_cell.get_T())
            _ghost.set_rho_from_TP()

            _ghost.update_vec_from_var()

        else:
            raise NotImplementedError('only left BC for inlet_uty')


class BC_Outlet_P(BoundaryCondition):

    def __init__(self, name, dict_bc):
        """constructor of BC_Outlet_P"""
        self.P_bc = None
        print("init BC_Outlet_P")
        BoundaryCondition.__init__(self, name, dict_bc)

    def apply_bc(self, field):
        """
        Apply boundary conditions

            ```
               ghost  BC  domain
            |-- q_0 --|-- q_1 --|
                      |<-- q_target

            q_t = 0.5 * (q_0 + q_1)
            q_0 = 2 q_t - q_1
            ```

        :param field: field object
        """

        if self.idx_cell_bc == -1:
            _ghost = field.lst_cell[self.idx_cell_bc]

            _ghost.update_vec_from_var()
            first_cell = field.lst_cell[self.idx_cell_bc - 1]
            left_left_cell = field.lst_cell[-2]
            if self.method_bc == "Dirichlet":
                _ghost.set_P(self.P_bc)
            elif self.method_bc == "flux":
                _ghost.set_P(2. * self.P_bc - first_cell.get_P())
            else:
                raise NotImplementedError("BC methods are 'Dirichlet' or 'flux'")

            _ghost.set_u(first_cell.get_u())
            _ghost.set_T(first_cell.get_T())
            _ghost.set_rho_from_TP()

        else:
            raise NotImplementedError('only right BC for outlet_P')


# def apply_BC(field, inp):
#     """
#     Apply boundary conditions
#
#            ghost  BC  domain
#         |-- q_0 --|-- q_1 --|
#                   |<-- q_target
#
#         q_t = 0.5 * (q_0 + q_1)
#         q_0 = 2 q_t - q_1
#     :param field: Field object
#     :return: modified field object
#     """
#
#     left_ghost = field.lst_cell[0]
#     right_cell = field.lst_cell[1]
#     right_right_cell = field.lst_cell[2]
#     left_ghost.set_u(2. * inp.bc_left_u - right_right_cell.get_u())
#     left_ghost.set_P(right_cell.get_P())
#     left_ghost.set_T(2. * inp.bc_left_T - right_right_cell.get_T())
#     left_ghost.set_rho_from_TP()
#
#     right_ghost = field.lst_cell[-1]
#     left_of_right_ghost = field.lst_cell[-2]
#     left_left_cell = field.lst_cell[-2]
#     right_ghost.set_P(2. * inp.bc_right_P - left_left_cell.get_P())
#     # right_ghost.set_P(inp.bc_right_P)
#     right_ghost.set_u(left_of_right_ghost.get_u())
#     right_ghost.set_T(left_of_right_ghost.get_T())
#     right_ghost.set_rho_from_TP()
#
#     if right_ghost.get_u() < 0.:
#         right_ghost.set_u(0)
#
#     for _ghost in [left_ghost, right_ghost]:
#         _ghost.update_vec_from_var()
#
#     return field

implemented_bc = ["inlet_UTY", "outlet_P"]

dict_BCType = {}
dict_BCType['inlet_UTY'] = BC_Inlet_UTY
dict_BCType['outlet_P'] = BC_Outlet_P


def init_BC(field, params_BC):
    for bc_name, bc_item in params_BC.items():
        if bc_item["type"] in implemented_bc:
            bc_type = bc_item["type"]
            bc_class = dict_BCType[bc_type]
            bc_init = bc_class(bc_name, bc_item)
            field.lst_BC.append(bc_init)
