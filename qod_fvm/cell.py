"""
Base element of the solver. Each Cell object is a control volume and all the elementary operations
are performed at this level.
"""
import numpy as np

from qod_fvm import utils

utils.plt_style()

N_TRANSPORT_EQ = 3


class Cell():
    """
    Control volume object
    """

    def __init__(self):
        """Cell constructor"""

        # Cell size [m]
        self.dx = 0.
        # Cell position
        self.x_i = 0.

        # area
        self.diam = 0.
        self.area = 0.

        # Left face position
        self.x_f_0 = 0.

        # right face position
        self.x_f_1 = 0.

        # How does area and volume work?
        self.vol = 0.
        # Do we have trapezoidal elements?

        self.normal_l = -1.
        self.normal_r = 1.

        # variables
        self.pres = 0.
        self.temp = 0.
        self.e_int = 0.
        self.e_tot = 0.
        self.u = 0.

        # self.cp = 0.
        self.gamma = 0.
        self.r_gas = 0.

        # cons
        self.rho = 0.
        self.rho_u = 0.
        self.rho_E = 0.

        # index of vector
        self.idx_mass = 0
        self.idx_momentum = 1
        self.idx_energy = 2

        # cons_vec
        self.w_cons = np.zeros(N_TRANSPORT_EQ)

        # source term
        self.s_cons = np.zeros(N_TRANSPORT_EQ)

        # Flux vector
        # cell-center
        self.f_cons = np.zeros(N_TRANSPORT_EQ)
        # left
        self.flux_face_l = np.zeros(N_TRANSPORT_EQ)
        # right
        self.flux_face_r = np.zeros(N_TRANSPORT_EQ)

        self.n_transport_eq = N_TRANSPORT_EQ

    def compute_volume(self):
        """Compute and set the volume of the CV"""
        self.vol = self.area * self.dx

    def set_positions(self, x_minus_half, x_center, x_plus_half):
        """
        Initialize the position of the cell center, the faces and the size of the cell

        :param x_minus_half: [m]
        :param x_center: [m]
        :param x_plus_half: [m]
        """
        self.x_i = x_center
        self.x_f_0 = x_minus_half
        self.x_f_1 = x_plus_half
        self.dx = self.x_f_1 - self.x_f_0

    def get_T(self):
        """
        Get temperature

        :return: temperature [K]
        """
        return self.temp

    def get_P(self):
        """
        Get Pressure

        :return: pressure [Pa]
        """
        return self.pres

    def get_u(self):
        """
        Get axial velocity

        :return: velocity [m/s]
        """
        return self.u

    def get_rhou(self):
        """
        Get rhou

        :return: rho [kg.m-2.s-1]
        """
        return self.rho_u

    def set_T(self, temperature):
        """
        Set the temperature

        :param temperature: temperature [K]

        """
        self.temp = temperature

    def set_P(self, pressure):
        """
        Set the pressure
        :param pressure: pressure [Pa]
        """
        self.pres = pressure

    def set_u(self, u):
        """
        Set the axial velocity

        :param u: axial velocity [m/s]
        """
        self.u = u

    def set_rho_from_TP(self, temp=None, pres=None):
        """
        Set the density [kg/m^3] from the temperature [K] and pressure [Pa] from the ideal gas law.
        """
        if not temp:
            temp = self.temp
        else:
            self.temp = temp
        if not pres:
            pres = self.pres
        else:
            self.pres = pres
        self.rho = pres / (self.r_gas * temp)

    def set_T_from_RP(self):
        """
        Set the temperature [K] from the density and the pressure (Ideal gas law)
        """
        self.temp = self.pres / (self.r_gas * self.rho)

    def get_cp(self):
        """
        Get the specific heat capacity [J/K/kg] at constant pressure

        :return: c_p [J/K/kg]
        """
        return self.gamma * self.get_cv()

    def get_cv(self):
        """
        Get the specific heat capacity [J/K/kg] at constant volume

        :return: c_v [J/K/kg]
        """
        return self.r_gas / (self.gamma - 1.)

    def get_enthalpy(self):
        """
        Compute and return the specific sensible enthalpy [J/kg].
        For a calorically perfect gas, the sensible enthalpy is computed as:

        .. math::
            h = c_p T

        :return: enthalpy [J/kg]
        """
        _h = self.get_cp() * self.temp
        return _h

    def get_internal_energy(self):
        """
        Compute and return the internal specific energy [J/kg]. This is computed from the enthalpy:

        .. math::
            e_{int} = h - \\frac{p}{\\rho}

        :return: internal specific energy [J/kg]
        """
        _h = self.get_enthalpy()
        self.e_int = _h - self.pres / self.rho
        return self.e_int

    def get_total_energy(self):
        """
        Compute and returns total energy [J/kg].

        .. math::
            e_{tot} = e_{int} + \\frac{u^2}{2} = h - \\frac{p}{\\rho} + \\frac{u^2}{2}

        :return: total specific energy [J/kg]
        """
        self.e_tot = self.get_internal_energy() + 0.5 * self.u ** 2
        return self.e_tot

    def get_sos(self):
        """
        Computes and returns speed of sound [m/s]. An assertion is made on the positivity of this
        quantity.

        .. math::
            a = \\sqrt{\\gamma r T}

        :return: local speed of sound
        """
        sos = np.sqrt(self.gamma * self.r_gas * self.temp)
        if self.temp < 0.:
            print("Problem, negative speed of sound here:")
            output_fields = self.get_outfield()
            output_fields['x'] = self.x_i
            for key, value in output_fields.items():
                print("\t %s\t:\t%e" % (key, value))
            assert (sos > 0)
        return sos

    def update_cons_vec(self):
        """
        Fill the vector of conservative variables:

        .. math::
            U[0] =  \\rho
            U[1] =  \\rho u
            U[2] =  \\rho E

        """
        self.w_cons[self.idx_mass] = self.rho
        self.w_cons[self.idx_momentum] = self.rho_u
        self.w_cons[self.idx_energy] = self.rho_E

    def update_flux_vec(self):
        """
        Fill the vector of fluxes.

        .. math::
            F[0] = A \\rho u
            F[1] = A \\rho u^2 + p
            F[2] = A u (\\rho E + p)
        """
        self.f_cons[self.idx_mass] = self.rho_u
        self.f_cons[self.idx_momentum] = self.rho * self.u ** 2.# + self.pres
        self.f_cons[self.idx_energy] = self.u * (self.rho_E + self.pres)
        # self.f_cons[self.idx_energy] = self.u * (self.rho_E)
        # self.f_cons *= self.area

    def prim_to_cons(self):
        """
        Convert primitive variables to conservative variables.

            Careful, at this point the conservative vector isn't necessary update yet.
            You should do it manually by calling `update_vec_from_var`.
        """
        # Mass
        # no need, rho is prim and cons

        # Momentum
        self.rho_u = self.rho * self.u

        # Energy
        self.rho_E = self.rho * self.get_total_energy()

    def cons_to_prim(self):
        """
        Retrieve primitive variables from conservative vector


        .. math::
            E = e_i + \\frac{u^2}{2} = h - \\frac{p}{\\rho} - \\frac{u^2}{2} = \\frac{1}{\\gamma - 1} \\frac{p}{\\rho} + \\frac{u^2}{2}

        .. math::
            p = \\rho (\\gamma - 1) \\left( E + \\frac{u^2}{2} \\right)


        """
        # verify this
        # mass, ok
        self.rho = self.w_cons[0] #/ self.area

        # momentum
        self.u = self.w_cons[1] / self.rho #/ self.area

        # energy
        self.e_tot = self.w_cons[2] / self.rho # / self.area

        # pressure
        self.pres = (self.gamma - 1.) * self.rho * (self.e_tot - 0.5 * self.u ** 2.)

        #  Temperature
        self.set_T_from_RP()

        # e_int
        _ = self.get_internal_energy()

    def update_vec_from_var(self):
        """
        Update conservative and flux vectors. First converts primitive variables to conservative
        then fill vectors.
        """
        self.prim_to_cons()
        self.update_cons_vec()
        self.update_flux_vec()

    def update_var_from_vec(self):
        """
        Converts conservative variables to primitives
        """
        self.cons_to_prim()


    def get_outfield(self):
        output_fields = {}
        tmp = self.get_u()
        output_fields['vel-x'] = tmp

        tmp = self.gamma
        output_fields['gamma'] = tmp

        tmp = self.r_gas
        output_fields['R_GAS'] = tmp

        tmp = self.rho_E
        output_fields['rhoE'] = tmp

        tmp = self.rho
        output_fields['rho'] = tmp

        tmp = self.pres
        output_fields['P'] = tmp

        tmp = self.temp
        output_fields['T'] = tmp

        tmp = self.rho_u
        output_fields['rhou'] = tmp

        tmp = self.dx
        output_fields['dx'] = tmp

        output_fields['A'] = tmp
        tmp = self.area

        return output_fields
