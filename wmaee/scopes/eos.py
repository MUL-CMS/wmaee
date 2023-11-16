"""
Convenient functions for fitting the energy-volume data.
"""

import numpy as np
from scipy.optimize import curve_fit
from typing import Tuple, Optional, Iterable, Callable
from numpy.typing import ArrayLike, NDArray


def eos_birch_murnaghan(V: ArrayLike, E0: float, V0: float, B0: float, Bp: float) -> NDArray:
    """
    Birch-Murnaghan equation of state.

    Parameters
    ----------
    V : ArrayLike
        Volume value(s).
    E0 : float
        Equilibrium energy.
    V0 : float
        Equilibrium volume.
    B0 : float
        Equilibrium bulk modulus.
    Bp : float
        Pressure derivative of the equilibrium bulk modulus.

    Returns
    -------
    NDArray
        Energies corresponding to the given volumes.

    Notes
    -----
    The Birch-Murnaghan equation of state is given by:
    E(V) = E0 + (9*B0*V0/16) * (((V0/V)**(2/3) - 1)**3 * Bp + ((V0/V)**(2/3) - 1)**2 * (6 - 4*(V0/V)**(2/3)))
    where E(V) is the energy as a function of volume (V).
    """
    return E0 + (9 * B0 * V0 / 16) * (((V0 / V)**(2/3) - 1)**3 * Bp + ((V0 / V)**(2/3) - 1)**2 * (6 - 4 * (V0 / V)**(2/3)))


def birch_murnaghan_fit(volumes: ArrayLike,
                        energies: ArrayLike,
                        p0: Optional[Iterable[float]] = None,
                        show_errors: bool = False,
                        eos_function: bool = False) -> Tuple[Iterable[float], Callable[[NDArray], NDArray]]:
    """
    Perform a Birch-Murnaghan fit on a given set of volumes and energies.

    Parameters
    ----------
    volumes : ArrayLike
        Volumes.
    energies : ArrayLike
        Energies.
    p0 : Optional[Iterable[float]], optional
        Starting parameters in the order E0_guess, V0_guess, B0_guess, Bp_guess, by default None.
        If not provided, default values are calculated based on the input data.
    show_errors : bool, optional
        Whether to print the standard deviation of the fitted parameters, by default False.
    eos_function : bool, optional
        Whether to return the EOS function as a callable, by default False.

    Returns
    -------
    Tuple[Iterable[float], Callable[[NDArray], NDArray]]
        Tuple containing the optimum parameters and, optionally, the EOS function.

    optimum : Iterable[float]
        Optimal parameters E0, V0, B0, Bp.
    eos_function : Callable[[NDArray], NDArray], optional
        The EOS function as a callable.

    Notes
    -----
    The Birch-Murnaghan equation of state is used for fitting.

    If `eos_function` is True, the function returns the EOS function as a callable.

    If `show_errors` is True, the function prints the standard deviation of the fitted parameters.
    
    If `p0` is not provided, default values are calculated as follows:
    - V0_guess: (np.amin(volumes) + np.amax(volumes))/2
    - E0_guess: np.amin(energies)
    - B0_guess: 1  # Presumably in eV/Ang^3; 160.2 in GPa
    - Bp_guess: 1
    """ 
    if p0 is None:
        V0_guess = (np.amin(volumes) + np.amax(volumes))/2
        E0_guess = np.amin(energies)
        B0_guess = 1  # Presumably in eV/Ang^3; 160.2 in GPa
        Bp_guess = 1
        p0 = (E0_guess, V0_guess, B0_guess, Bp_guess)

    optimum, p_cov = curve_fit(eos_birch_murnaghan, np.array(volumes), np.array(energies), p0=p0)
    
    if show_errors:
        p_err = np.sqrt(np.diag(p_cov))
        print(f'Standard deviation of E0, V0, B0, Bp = {p_err}')

    if eos_function:
        return optimum, lambda V: eos_birch_murnaghan(V, *optimum)
    
    return optimum



def eos_murnaghan(V: ArrayLike, E0: float, V0: float, B0: float, Bp: float) -> NDArray:
    """
    Murnaghan equation of state.

    Parameters
    ----------
    V : ArrayLike
        Volume value(s).
    E0 : float
        Equilibrium energy.
    V0 : float
        Equilibrium volume.
    B0 : float
        Equilibrium bulk modulus.
    Bp : float
        Pressure derivative of the equilibrium bulk modulus.

    Returns
    -------
    NDArray
        Energies corresponding to the given volumes.

    Notes
    -----
    The Murnaghan equation of state is given by:
    E(V) = E0 + B0 * V0 * (1/(Bp*(Bp-1)) * (V/V0)**(1-Bp) + V/(Bp*V0) - 1/(Bp-1))
    where E(V) is the energy as a function of volume (V).
    """
    return E0 + B0 * V0 * (1/(Bp * (Bp - 1)) * (V / V0)**(1 - Bp) + V / (Bp * V0) - 1 / (Bp - 1))



def murnaghan_fit(volumes: ArrayLike,
                  energies: ArrayLike,
                  p0: Optional[Iterable[float]] = None,
                  show_errors: bool = False,
                  eos_function: bool = False) -> Tuple[Iterable[float], Callable[[NDArray], NDArray]]:
    """
    Perform a Murnaghan fit on a given set of volumes and energies.

    Parameters
    ----------
    volumes : ArrayLike
        Volumes.
    energies : ArrayLike
        Energies.
    p0 : Optional[Iterable[float]], optional
        Starting parameters in the order E0_guess, V0_guess, B0_guess, Bp_guess, by default None.
        If not provided, default values are calculated based on the input data.
    show_errors : bool, optional
        Whether to print the standard deviation of the fitted parameters, by default False.
    eos_function : bool, optional
        Whether to return the EOS function as a callable, by default False.

    Returns
    -------
    Tuple[Iterable[float], Callable[[NDArray], NDArray]]
        Tuple containing the optimum parameters and, optionally, the EOS function.

    optimum : Iterable[float]
        Optimal parameters E0, V0, B0, Bp.
    eos_function : Callable[[NDArray], NDArray], optional
        The EOS function as a callable.

    Notes
    -----
    The Murnaghan equation of state is used for fitting.

    If `eos_function` is True, the function returns the EOS function as a callable.

    If `show_errors` is True, the function prints the standard deviation of the fitted parameters.

    If `p0` is not provided, default values are calculated as follows:
    - V0_guess: (np.amin(volumes) + np.amax(volumes))/2
    - E0_guess: np.amin(energies)
    - B0_guess: 1  # Presumably in eV/Ang^3; 160.2 in GPa
    - Bp_guess: 1
    """
    if p0 is None:
        V0_guess = (np.amin(volumes) + np.amax(volumes))/2
        E0_guess = np.amin(energies)
        B0_guess = 1  # Presumably in eV/Ang^3; 160.2 in GPa
        Bp_guess = 1
        p0 = (E0_guess, V0_guess, B0_guess, Bp_guess)

    optimum, p_cov = curve_fit(eos_murnaghan, np.array(volumes), np.array(energies), p0=p0)
    
    if show_errors:
        p_err = np.sqrt(np.diag(p_cov))
        print(f'Standard deviation of E0, V0, B0, Bp = {p_err}')

    if eos_function:
        return optimum, lambda V: eos_murnaghan(V, *optimum)
    
    return optimum


def polynomial_fit(volumes: ArrayLike,
                   energies: ArrayLike,
                   order: int = 3,
                   eos_function: bool = False) -> Tuple[float, float, float, Optional[Callable[[NDArray], NDArray]]]:
    """
    Perform a polynomial fit on a given set of volumes and energies.

    Parameters
    ----------
    volumes : ArrayLike
        Volumes.
    energies : ArrayLike
        Energies.
    order : int, optional
        Order of the polynomial fit, by default 3.
    eos_function : bool, optional
        Whether to return the EOS function as a callable, by default False.

    Returns
    -------
    Tuple[float, float, float, Optional[Callable[[NDArray], NDArray]]]
        Tuple containing the equilibrium energy, equilibrium volume, bulk modulus,
        and, optionally, the EOS function.

    E_eq : float
        Equilibrium energy.
    V_eq : float
        Equilibrium volume.
    B_eq : float
        Bulk modulus.
    eos_function : Callable[[NDArray], NDArray], optional
        The EOS function as a callable.

    Notes
    -----
    The function performs a polynomial fit of the data and calculates equilibrium properties.

    If `eos_function` is True, the function returns the EOS function as a callable.
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
    # to get the equilibrium energy
    E_eq = np.polyval(coefficients, V_eq)

    # Now, we will calculate the bulk modulus, B_eq
    # In general, B_eq = - V_eq * (dp/dV)
    # p(V) = first derivative, hence dp/dV is the second
    # derivative of the energy polynomial
    second_derivative = np.polyder(coefficients, 2)

    # Calculating the bulk modulus
    B_eq = -V_eq * (-np.polyval(second_derivative, V_eq))

    fitted_params = (E_eq, V_eq, B_eq)
    
    if eos_function:
        return fitted_params, lambda V: np.polyval(coefficients, V)
    
    return fitted_params
