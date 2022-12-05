"""
Convenient functions for fitting the energy-volume calculations.
"""

import numpy as np
import scipy.optimize as opt
from typing import Sequence
import warnings


class EqosBirchMurnaghan:

    # __init__ function creates an Object from the class EqosBirchMurnaghan
    # and initializes the object's attributes
    def __init__(self, volumes: np.ndarray, energies: np.ndarray, initial_value: Sequence = None) -> None:
        self.volumes = volumes
        self.energies = energies

        self.initial_value = initial_value

        if self.initial_value is None:
            emin = np.min(self.energies)  # guessing value for energy minimum
            # guessing value for volume minimum
            V0 = self.volumes[np.argmin(self.energies)]
            # guessing value for bulk modulus = 0.01
            self.initial_value = (V0, 0.01, 0.1, emin)
            # guessing value for pressure derivative of bulk modulus = 0.1

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # p_out is the variable, which stores the fitted parameters in the order V0, K0, Kp, and E0
            # pcov indicates the qualtiy of the fit
            p_out, pcov = opt.curve_fit(self.func, self.volumes, self.energies,
                                        jac=self.jac,
                                        p0=self.initial_value,
                                        xtol=1e-10,
                                        maxfev=5000)

        # Extracts the equilibrium properties of the fit from p_out:
        self.V0 = p_out[0]
        self.K0 = p_out[1]
        self.Kp = p_out[2]
        self.E0 = p_out[3]

        self.energy = np.vectorize(self._energy)
        self.pressure = np.vectorize(self._pressure)
        self.bulk_modulus = np.vectorize(self._bulkmodulus)

    # Defines Birch-Murnaghan EOS function
    def func(self, V, V0, K0, Kp, E0):
        eta = (V0 / V) ** (2.0 / 3.0)
        E = E0 + 9.0 * K0 * V0 / 16.0 * (eta - 1.0) ** 2 \
            * (Kp * (eta - 1.0) + 6.0 - 4.0 * eta)
        return E

    def jac(self, V, V0, K0, Kp, E0):
        eta = (V0 / V) ** (2 / 3)

        dfuncdV0 = (3 * K0 * (V0 ** (2 / 3) - V ** (2 / 3)) * ((3 * Kp - 12) * V0 ** (5 / 3) + (
            (6 * Kp - 24) * V0 ** (2 / 3) + (62 - 12 * Kp) * V ** (2 / 3)) * V0 + (3 * Kp - 18) *
            V ** (4 / 3) * V0 ** (1 / 3))) / (
            16 * V ** 2 * V0 ** (1 / 3))

        dfuncdK0 = 9 * V0 / 16.0 * (eta - 1.0) ** 2 \
            * (6.0 + Kp * (eta - 1) - 4 * eta)

        dfuncdKp = 9 * V0 * K0 / 16 * ((eta - 1) ** 3)

        dfuncdE0 = np.ones_like(V)

        return np.stack((dfuncdV0, dfuncdK0, dfuncdKp, dfuncdE0), axis=-1)

    # Calculates the bulk modulus
    def _bulkmodulus(self, volume):
        V = volume
        V0 = self.V0
        K0 = self.K0
        Kp = self.Kp

        return (V) * -(K0 * V0 ** (5 / 3) * ((144 - 27 * Kp) * V ** (14 / 3) + (
            (12 * Kp - 64) * V ** (2 / 3) + (42 * Kp - 196) * V0 ** (2 / 3)) * V ** 4 + (
            108 - 27 * Kp) * V0 ** (4 / 3) * V ** (10 / 3))) / (
            8 * V ** (22 / 3))

    # Calculates energy with self.func for an arbitrary volume
    def _energy(self, volume):
        return self.func(volume, self.V0, self.K0, self.Kp, self.E0)

    # Calculates the final pressure acting on the system
    # If you taka the equilibrium volume V0 as an input the calculated pressure should be exactly 0!
    def _pressure(self, volume):
        V = volume
        V0 = self.V0
        K0 = self.K0
        Kp = self.Kp

        return 3 * K0 / 2 * ((V0 / V) ** (7 / 3) - (V0 / V) ** (5 / 3)) * (
            1 + 3 / 4 * (Kp - 4) * ((V0 / V) ** (2 / 3) - 1))

    # Calculates the pressure derivative of the bulk modulus
    def _deriv_bulkmodulus(self, volume):
        V = volume
        V0 = self.V0
        K0 = self.K0
        Kp = self.Kp

        return -(K0 * V0 ** (5 / 3) * ((48 - 9 * Kp) * V ** (11 / 3) + ((6 * Kp - 32)
                                                                        * V ** (2 / 3) + (12 * Kp - 56) * V0 ** (
            2 / 3)) * V ** 3 + (
            36 - 9 * Kp) * V0 ** (4 / 3) * V ** (7 / 3))) / (4 * V ** (16 / 3))

    # Functions for extracting properties from the eos fit:
    @property
    def V_eq(self):
        return self.V0

    @property
    def E_eq(self):
        return self.E0

    @property
    def K_eq(self):
        return self.K0

    @property
    def Kp_eq(self):
        return self.Kp


class EqosMurnaghan:
    """
    Implements a fit using the Murnaghan equation of state.
    """

    def __init__(self, volumes, energies):
        self.volumes = volumes
        self.energies = energies

        self.v_min = np.amin(volumes)
        self.v_max = np.amax(volumes)

        # Set initial parameters
        i_min = np.argmin(self.energies)
        e0 = self.energies[i_min]
        v0 = self.volumes[i_min]
        K = 1.0
        Kp = 3.5

        p_in = (e0, v0, K, Kp)

        p_out, pcov = opt.curve_fit(self.fun, volumes, energies, p0=p_in)

        self.E0 = p_out[0]
        self.V0 = p_out[1]
        self.K = p_out[2]
        self.Kp = p_out[3]

        # Define convenience methods
        self.energy = np.vectorize(lambda v:
                                   self.fun(v, self.E0, self.V0, self.K, self.Kp))

        self.pressure = np.vectorize(self._pressure)
        self.bulk_modulus = np.vectorize(self._bulk_modulus)

    def fun(self, v, E0, V0, K, Kp):
        nu = v / V0
        return E0 + K * V0 / Kp * (
            (nu**(-Kp + 1) - 1) / (Kp - 1) + nu - 1)

    @property
    def V_eq(self):
        return self.V0

    @property
    def E_eq(self):
        return self.E0

    @property
    def K_eq(self):
        return self.K

    @property
    def Kp_eq(self):
        return self.Kp

    def _bulk_modulus(self, v):
        return self.K * (v / self.V0)**(-self.Kp)

    def _pressure(self, v):
        return self.K / self.Kp * (
            (v / self.V0)**(-self.Kp) - 1)


class EqosPolynomial:
    """
    Implements a fit using the polynomial equation of state.
    """

    def __init__(self, volumes, energies, order=3):
        self.volumes = volumes
        self.energies = energies

        self.v_min = np.amin(volumes)
        self.v_max = np.amax(volumes)

        # Perform the fit
        self.p = np.polyfit(volumes, energies, order)

        self.pd = np.polyder(self.p)
        self.pd2 = np.polyder(self.pd)

        roots = np.roots(self.pd)

        i_min = np.where((self.v_min <= roots) & (roots <= self.v_max))[0]
        if not i_min:
            raise RuntimeError(
                "Equation of state seems not to contain a minimum")

        self.V0 = roots[i_min][0].real
        self.E0 = np.polyval(self.p, self.V0)

        self.energy = np.vectorize(self._energy)
        self.pressure = np.vectorize(self._pressure)
        self.bulk_modulus = np.vectorize(self._bulk_modulus)

    @property
    def V_eq(self):
        return self.V0

    @property
    def E_eq(self):
        return self.E0

    @property
    def K_eq(self):
        return self.bulk_modulus(self.V0)

    def _energy(self, v):
        return np.polyval(self.p, v)

    def _pressure(self, v):
        return -np.polyval(self.pd, v)

    def _bulk_modulus(self, v):
        return np.polyval(self.pd2, v) * v
