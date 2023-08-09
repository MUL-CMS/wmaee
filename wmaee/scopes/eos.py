"""
Convenient functions for fitting the energy-volume calculations.
"""


import numpy as np
from typing import *
from scipy.optimize import curve_fit
from numpy.typing import ArrayLike, NDArray


def murnaghan_fit(volumes: np.ndarray, energies: np.ndarray):
    """
    This function takes two arguments in np.ndarray format, 
    performs a Murnaghan EoS fit, and returns 
    E_eq, V_eq, B_eq, Bp_eq, V_fit, E_fit.

    E_eq ... equilibrium energy.
    V_eq ... equilibrium volume.
    B_eq ... equilibrium bulk modulus.
    Bp_eq .. pressure derivative of the equilibrium bulk modulus.

    V_fit and E_fit are generated volumes and with the fitted
    values E_eq, V_eq, B_eq, Bp_eq calculated corresponding energies.
    V_fit and E_fit can then be used for plotting the fit. 
    """

    # Defining the Murnaghan equation of state function
    def murnaghan_eos(V, E0, V0, B0, Bp):
        return E0 + B0*V0*(1/(Bp*(Bp-1))*(V/V0)**(1-Bp) + V/(Bp*V0) - 1/(Bp-1))

    # Finding V0 and E0
    E0 = np.min(energies)
    min_index = np.argmin(energies)
    V0 = volumes[min_index]

    # Initial guesses
    p_in = [E0, V0, 1, 3.5]

    # Fitting the data to the Murnaghan equation of state
    p_opt, p_cov = curve_fit(murnaghan_eos, volumes, energies, p0=p_in)

    # Extracting fitted parameters
    E_eq = p_opt[0]
    V_eq = p_opt[1]
    B_eq = p_opt[2]
    Bp_eq = p_opt[3]

    # Calculating the standard deviation of the fitted parameters
    p_err = np.sqrt(np.diag(p_cov))

    print(f'Standard deviation of E0, V0, B0, Bp = {p_err}')

    # Preparing the values for the plot
    V_fit = np.linspace(np.min(volumes), np.max(volumes), 100)
    E_fit = murnaghan_eos(V_fit, E_eq, V_eq, B_eq, Bp_eq)

    return E_eq, V_eq, B_eq, Bp_eq, V_fit, E_fit

###############################################################################


def birch_murnaghan_fit_(volumes: ArrayLike, energies: ArrayLike) -> Tuple[float, float, float, float, NDArray, NDArray]:
    """
    This function takes two arguments in np.ndarray format, 
    performs a Birch-Murnaghan EoS fit, and returns 
    E_eq, V_eq, B_eq, Bp_eq, V_fit, E_fit.


    E_eq ... equilibrium energy.
    V_eq ... equilibrium volume.
    B_eq ... equilibrium bulk modulus.
    Bp_eq .. pressure derivative of the equilibrium bulk modulus.

    V_fit and E_fit are generated volumes and with the fitted
    values E_eq, V_eq, B_eq, Bp_eq calculated corresponding energies.
    V_fit and E_fit can then be used for plotting the fit. 
    """

    # Defining the Birch-Murnaghan equation of state function
    def birch_murnaghan_eos(V, E0, V0, B0, Bp):
        return E0 + (9*B0*V0/16)*(((V0/V)**(2/3)-1)**3*Bp + ((V0/V)**(2/3)-1)**2*(6-4*(V0/V)**(2/3)))

    # Finding V0 and E0
    E0 = np.min(energies)
    min_index = np.argmin(energies)
    V0 = volumes[min_index]

    # Initial guesses
    p_in = [E0, V0, 1, 3.5]

    # Fitting the data to the Murnaghan equation of state
    p_opt, p_cov = curve_fit(birch_murnaghan_eos, volumes, energies, p0=p_in)

    # Extracting fitted parameters
    E_eq, V_eq, B_eq, Bp_eq = p_opt

    # Calculating the standard deviation of the fitted parameters
    p_err = np.sqrt(np.diag(p_cov))

    print(f'Standard deviation of E0, V0, B0, Bp = {p_err}')

    # Preparing the values for the plot
    V_fit = np.linspace(np.min(volumes), np.max(volumes), 100)
    E_fit = birch_murnaghan_eos(V_fit, E_eq, V_eq, B_eq, Bp_eq)

    return E_eq, V_eq, B_eq, Bp_eq, V_fit, E_fit


def polynomial_fit(volumes: ArrayLike, energies: ArrayLike, order: int = 3):
    """
    This function takes two arguments in np.ndarray format, 
    performs a polynomial fit, and returns 
    E_eq, V_eq, B_eq, V_fit, E_fit.
    The third parameter is the order of the polynomial and is optional.
    The order should be higher than 2, because the E-V curves 
    typically don't have a quadratic shape! 

    E_eq ... equilibrium energy.
    V_eq ... equilibrium volume.
    B_eq ... equilibrium bulk modulus.

    V_fit and E_fit are generated volumes and with the fitted
    values E_eq, V_eq, B_eq, Bp_eq calculated corresponding energies.
    V_fit and E_fit can then be used for plotting the fit. 
    """

    # Extract the smallest and the biggest volume
    V_min = np.amin(volumes)
    V_max = np.amax(volumes)

    # Polynomial fit of the data
    coefficients = np.polyfit(volumes, energies, order)

    # Calculate the coefficients of the first derivative
    first_derivative = np.polyder(coefficients, 1)

    # Calculate at which volume the first derivative is zero
    roots = np.roots(first_derivative)

    # Now we need to extract indices of the correct
    # minimum volume from roots
    v_min_indices = np.where((V_min <= roots) & (roots <= V_max))

    # With v_min_indices we can extract the correct minimum
    # volume from roots
    V_eq = float(roots[v_min_indices])

    # Further, we can now evaluate the polynomial at V_eq
    # to get the equilubrium energy
    E_eq = np.polyval(coefficients, V_eq)

    # Now, we will calculate the bulk modulus, B0
    # In general, B = - V*(dp/dV)
    # p(V) = first derivative, hence dp/dV is the second
    # derivative of the energy polynomial
    second_derivative = np.polyder(coefficients, 2)

    # Calculating the bulk modulus
    B_eq = -V_eq * (-np.polyval(second_derivative, V_eq))

    # Preparing the values for the plot
    V_fit = np.linspace(np.min(volumes), np.max(volumes), 100)
    E_fit = np.polyval(coefficients, V_fit)

    return E_eq, V_eq, B_eq, V_fit, E_fit


def birch_murnaghan(volumes: ArrayLike, e0: float, v0: float, b0: float, bp: float):
    two_thirds = 2/3
    return e0 + (9 * v0 * b0 / 16) * (bp * ((volumes / v0) ** two_thirds - 1) ** 3 + ((volumes / v0) ** two_thirds - 1) ** 2 * (6 - 4 * (volumes / v0) ** two_thirds))


def birch_murnaghan_fit(volumes: ArrayLike, energies: ArrayLike, p0: Optional[Iterable[float]] = None) -> Tuple[Iterable[float], Callable[[NDArray], NDArray]]:
    """
    Performs a Birch-Murnaghan-fit on a given set of volumes and energies
    :param volumes: volumes
    :type volumes: ArrayLike
    :param energies: energies
    :type energies: ArrayLike
    :param p0: (iterable) starting parameters in the order E0_guess, V0_guess, B0_guess, Bp_guess (default: None)
    :param return_analytic: (bool) wether also to return a SymPy analytic expression of the fit
    :returns: (optimum:np.array, E(V):callable, optional analytic:Sympy.Expression)
    """

    if p0 is None:
        V0_guess = np.amin(volumes)
        E0_guess = np.amin(energies)
        B0_guess = 100
        Bp_guess = 1
        p0 = (E0_guess, V0_guess, B0_guess, Bp_guess)

    optimum, _ = curve_fit(birch_murnaghan, np.array(volumes), np.array(energies), p0=p0)
    return optimum, lambda volumes: birch_murnaghan(volumes, *optimum)
