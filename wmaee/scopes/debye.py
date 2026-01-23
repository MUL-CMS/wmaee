"""
Collection of routines for thermodynamics with the Debye model.

David Holec
david.holec@unileoben.ac.at
"""

from typing import Optional, Union, Any
from copy import deepcopy
from math import exp
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from ase import Atoms


def get_rho(V: float, material: Union[str, Atoms]) -> float:
    """
    Compute mass density rho [g/cm^3] from per-atom volume V [Å^3/atom]
    and either an ASE Atoms object or a chemical formula.

    Parameters
    ----------
    V : float
        Per-atom volume in Å^3/atom.
    material : str or ase.Atoms
        Chemical formula string (e.g., 'Fe2O3', 'Si', 'AlN') or an ASE Atoms object.

    Returns
    -------
    float
        Mass density in g/cm^3.

    Notes
    -----
    - Conversions: 1 amu = 1.66054e-24 g and 1 Å^3 = 1e-24 cm^3,
      so rho = 1.66054 * (average_mass_in_amu) / V.
    """
    if isinstance(material, str):
        atoms = Atoms(material)
    elif isinstance(material, Atoms):
        atoms = material
    else:
        raise TypeError("material must be a chemical formula string or an ase.Atoms object")

    M_amu = atoms.get_masses().mean()  # average mass per atom (amu)
    return 1.66054 * M_amu / V


    
    
class Debye:
    """
    Implementation of the harmonic Debye approximation.

    This class stores and provides accessors for thermodynamic and elastic
    descriptors required to evaluate properties in the harmonic Debye model,
    such as the Helmholtz free energy, heat capacity, and entropy.

    Parameters
    ----------
    at_per_fu : int, optional
        Number of atoms in the primitive cell, $N_\\mathrm{at}$, by default 1.

    Attributes
    ----------
    data : dict
        Internal storage for material fields:
        - 'V': float or None
            Volume per atom, $V$ [Å^3/atom].
        - 'rho': float or None
            Mass density, $\\rho$ [g/cm^3].
        - 'Nat': int
            Number of atoms in the primitive cell, $N_\\mathrm{at}$ [1].
        - 'B': float or None
            Bulk modulus, $B$ [GPa].
        - 'nu': float or None
            Poisson's ratio, $\\nu$ [1].
        - 'v': float or None
            Mean (Debye) sound velocity, $v$ [m/s].
        - 'ThetaD': float or None
            Debye temperature, $\\Theta_D$ [K].
        - 'E0': float or None
            Equilibrium total energy per atom, $E_0$ [eV/atom].

    Notes
    -----
    In the harmonic Debye approximation, vibrational contributions depend
    primarily on the Debye temperature $\\Theta_D$ (or equivalently the
    Debye frequency), which can be related to the elastic constants (here
    represented by $B$ and $\\nu$) and the density $\\rho$ through the mean
    sound velocity $v$. The model is commonly used to approximate the
    vibrational free energy
    $$
    F_\\mathrm{vib}(T) = k_B T
    \\left[
      3 N_\\mathrm{at} \\ln\\!\\left(1 - e^{-\\Theta_D/T}\\right)
      - 9 N_\\mathrm{at} \\frac{T}{\\Theta_D} \\int_0^{\\Theta_D/T} \\frac{x^3}{e^x - 1}\\,dx
    \\right],
    $$
    and the total Helmholtz free energy can be expressed as
    $F(T) = E_0 + F_\\mathrm{vib}(T)$ at fixed volume in the harmonic limit.

    Unit conventions are per-atom (where applicable) to avoid ambiguity when
    comparing different structures or formula units.
    """

    def __init__(self, at_per_fu: int = 1):
        self.data = {
            'V': None,        # float: Å^3/atom
            'rho': None,      # float: g/cm^3
            'Nat': at_per_fu, # int: atoms in primitive cell
            'B': None,        # float: GPa
            'nu': None,       # float: dimensionless
            'v': None,        # float: m/s
            'ThetaD': None,   # float: K
            'E0': None        # float: eV/atom
        }

    # -------------------------
    # Volume V [Å^3/atom]
    # -------------------------
    @property
    def V(self) -> Optional[float]:
        """
        Volume per atom, $V$.

        Returns
        -------
        float or None
            Volume per atom $V$ [Å^3/atom], or ``None`` if unset.

        Notes
        -----
        - Use per-atom volume to avoid dependence on the choice of the
          crystallographic cell or formula unit size.
        """
        return self.data['V']

    @V.setter
    def V(self, V: float) -> None:
        """
        Set the volume per atom.

        Parameters
        ----------
        V : float
            Volume per atom $V$ [Å^3/atom].

        Notes
        -----
        - No validation is performed here; consider enforcing $V > 0$ upstream.
        """
        self.data['V'] = V

    # -------------------------
    # Density rho [g/cm^3]
    # -------------------------
    @property
    def rho(self) -> Optional[float]:
        """
        Mass density, $\\rho$.

        Returns
        -------
        float or None
            Mass density $\\rho$ [g/cm^3], or ``None`` if unset.
        """
        return self.data['rho']

    @rho.setter
    def rho(self, rho: float) -> None:
        """
        Set the mass density.

        Parameters
        ----------
        rho : float
            Mass density $\\rho$ [g/cm^3].

        Notes
        -----
        - In Debye-model workflows, $\\rho$ and elastic properties determine
          the mean sound velocity $v$ (and thus $\\Theta_D$).
        """
        self.data['rho'] = rho

    # -------------------------
    # Number of atoms Nat [1]
    # -------------------------
    @property
    def Nat(self) -> int:
        """
        Number of atoms in the primitive cell, $N_\\mathrm{at}$.

        Returns
        -------
        int
            $N_\\mathrm{at}$ [1].
        """
        return self.data['Nat']

    @Nat.setter
    def Nat(self, Nat: int) -> None:
        """
        Set the number of atoms in the primitive cell.

        Parameters
        ----------
        Nat : int
            $N_\\mathrm{at}$ [1].

        Notes
        -----
        - Should be a positive integer.
        """
        self.data['Nat'] = Nat

    # -------------------------
    # Bulk modulus B [GPa]
    # -------------------------
    @property
    def B(self) -> Optional[float]:
        """
        Bulk modulus, $B$.

        Returns
        -------
        float or None
            Bulk modulus $B$ [GPa], or ``None`` if unset.
        """
        return self.data['B']

    @B.setter
    def B(self, B: float) -> None:
        """
        Set the bulk modulus.

        Parameters
        ----------
        B : float
            Bulk modulus $B$ [GPa].
        """
        self.data['B'] = B

    # -------------------------
    # Poisson's ratio nu [1]
    # -------------------------
    @property
    def nu(self) -> Optional[float]:
        """
        Poisson's ratio, $\\nu$.

        Returns
        -------
        float or None
            Poisson's ratio $\\nu$ [dimensionless], or ``None`` if unset.

        Notes
        -----
        - For mechanically stable, isotropic media, $-1 < \\nu < 0.5$.
        """
        return self.data['nu']

    @nu.setter
    def nu(self, nu: float) -> None:
        """
        Set Poisson's ratio.

        Parameters
        ----------
        nu : float
            Poisson's ratio $\\nu$ [dimensionless].
        """
        self.data['nu'] = nu

    # -------------------------
    # Mean sound velocity v [m/s]
    # -------------------------
    @property
    def v(self) -> Optional[float]:
        """
        Mean (Debye) sound velocity, $v$.

        Returns
        -------
        float or None
            Mean sound velocity $v$ [m/s], or ``None`` if unset.

        Notes
        -----
        - Often derived from elastic constants and density
          (e.g., via averages of longitudinal and transverse speeds).
        """
        return self.data['v']

    @v.setter
    def v(self, v: float) -> None:
        """
        Set the mean (Debye) sound velocity.

        Parameters
        ----------
        v : float
            Mean sound velocity $v$ [m/s].
        """
        self.data['v'] = v

    # -------------------------
    # Debye temperature ThetaD [K]
    # -------------------------
    @property
    def ThetaD(self) -> Optional[float]:
        """
        Debye temperature, $\\Theta_D$.

        Returns
        -------
        float or None
            Debye temperature $\\Theta_D$ [K], or ``None`` if unset.

        Notes
        -----
        - $\\Theta_D$ can be computed from $v$ and number density.
        - In the harmonic Debye model, $\\Theta_D$ sets the vibrational
          energy scale that governs $F_\\mathrm{vib}(T)$ and $C_V(T)$.
        """
        return self.data['ThetaD']

    @ThetaD.setter
    def ThetaD(self, ThetaD: float) -> None:
        """
        Set the Debye temperature.

        Parameters
        ----------
        ThetaD : float
            Debye temperature $\\Theta_D$ [K].
        """
        self.data['ThetaD'] = ThetaD

    # -------------------------
    # Equilibrium energy E0 [eV/atom]
    # -------------------------
    @property
    def E0(self) -> Optional[float]:
        """
        Equilibrium total energy per atom, $E_0$.

        Returns
        -------
        float or None
            $E_0$ [eV/atom], or ``None`` if unset.

        Notes
        -----
        - In the harmonic limit at fixed volume, the Helmholtz free energy is
          $F(T) = E_0 + F_\\mathrm{vib}(T)$.
        """
        return self.data['E0']

    @E0.setter
    def E0(self, E0: float) -> None:
        """
        Set the equilibrium total energy per atom.

        Parameters
        ----------
        E0 : float
            $E_0$ [eV/atom].
        """
        self.data['E0'] = E0

    
    def calculate_mean_sound_velocity(self) -> None:
        """
        Compute and store the mean (Debye) sound velocity, v.
    
        The mean sound velocity is estimated from Poisson's ratio and the
        bulk modulus and density via:
        $$
        v = f(\nu)\,\sqrt{B/\rho},
        $$
        where
        $$
        f(\nu) =
        \left[
            \frac{3}{
                2\left(\tfrac{2}{3}\tfrac{1+\nu}{1-2\nu}\right)^{3/2}
              + \left(\tfrac{1}{3}\tfrac{1+\nu}{1-\nu}\right)^{3/2}
            }
        \right]^{1/3}.
        $$
    
        Parameters
        ----------
        None
    
        Returns
        -------
        None
            The computed mean sound velocity is stored in ``self.data['v']`` [m/s].
    
        Notes
        -----
        - Units:
          - Input: ``B`` in [GPa], ``rho`` in [g/cm^3], ``nu`` dimensionless.
          - Output: ``v`` in [m/s].
          - The factor ``1e3`` accounts for the conversion
            $\\sqrt{\\text{GPa}/(\\text{g cm}^{-3})} \\to \\text{m s}^{-1}$,
            since $1\\,\\text{GPa}=10^9\\,\\text{Pa}$ and
            $1\\,\\text{g cm}^{-3}=10^3\\,\\text{kg m}^{-3}$.
        - This method updates ``self.data['v']`` in place.
        - Requires that ``nu``, ``B``, and ``rho`` are set beforehand.
        """
        v = self.data['nu']
        fnu = (3/(2*(2/3*(1+v)/(1-2*v))**(3/2) + (1/3*(1+v)/(1-v))**(3/2)))**(1/3)
        self.data['v'] = 1e3 * fnu * (self.data['B']/self.data['rho'])**0.5
    
    
    def calculate_Debye_temperature(self) -> None:
        """
        Compute and store the Debye temperature, Θ_D.
    
        The Debye temperature is evaluated as
        $$
        \\Theta_D \\,=\\, \\frac{h}{k_B}\\, v \\,\\left(\\frac{3}{4\\pi\\,\\Omega}\\right)^{1/3},
        $$
        where $\\Omega = N_\\mathrm{at}\\,V$ is the primitive-cell volume and
        $v$ is the mean sound velocity.
    
        Parameters
        ----------
        None
    
        Returns
        -------
        None
            The computed Debye temperature is stored in ``self.data['ThetaD']`` [K].
    
        Notes
        -----
        - Units:
          - Input: ``V`` in [Å^3/atom], ``Nat`` dimensionless (atoms/cell),
            so $\\Omega$ is in [Å^3/cell]; ``v`` in [m/s].
          - Constants: ``h`` in [J·s], ``kB`` in [J·K^-1].
          - Output: ``ThetaD`` in [K].
          - The factor ``1e10`` converts [Å^-1] to [m^-1].
        - This method calls ``calculate_mean_sound_velocity()`` to ensure ``v`` is up to date.
        - Requires that ``V``, ``Nat``, and the inputs needed for ``v`` are set.
        """
        kB = 1.38064852e-23  # J K^-1  (m2 kg s-2 K-1)
        h = 6.62607015e-34   # J s     (m2 kg s^-1)
        from math import pi
    
        self.calculate_mean_sound_velocity()
        Omega = self.data['Nat'] * self.data['V']  # Å^3 per cell
        self.data['ThetaD'] = 1e10 * h/kB * (3/(4*pi*Omega))**(1/3) * self.data['v']
    
    
    def DebyeIntegral(self, TD: float | list[float], order: int = 3) -> float | list[float]:
        """
        Evaluate the Debye integral of order n.
    
        This computes the dimensionless Debye function
        $$
        D_n(x) \\,=\\, \\frac{n}{x^n} \\int_0^x \\frac{t^n}{e^t - 1}\\,dt,
        $$
        where typically $n=3$ and $x=\\Theta_D/T$.
    
        Parameters
        ----------
        TD : float or array-like
            Upper bound $x = \\Theta_D/T$ (dimensionless).
        order : int, optional
            Debye integral order $n$ (default is 3).
    
        Returns
        -------
        float or list of float
            Value(s) of $D_n(x)$ at the provided bound(s). Returns a float if
            ``TD`` is scalar, otherwise a list of floats.
    
        Notes
        -----
        - Numerical evaluation is performed via ``scipy.integrate.quad``.
        - This helper is used in the vibrational free-energy expression of
          the harmonic Debye model.
        """
    
        def __f(x):
            return x**order/(exp(x)-1)
    
        if hasattr(TD, '__len__'):
            return [order/(TDD**order) * integrate.quad(__f, 0, TDD)[0] for TDD in TD]
        else:
            return order/(TD**order) * integrate.quad(__f, 0, TD)[0]
    
    
    def Fvib_harm(self, T: float) -> float:
        """
        Vibrational Helmholtz free energy in the harmonic Debye model, per atom.
    
        The vibrational contribution is computed as
        $$
        F_\\mathrm{vib}(T)
          \\,=\\, \\frac{9}{8} k_B \\Theta_D
          \\, - \\, k_B T
          \\left[
            D_3\\!\\left(\\frac{\\Theta_D}{T}\\right)
            \\, -\\, 3\\ln\\!\\left(1 - e^{-\\Theta_D/T}\\right)
          \\right],
        $$
        where $D_3(x)$ is the Debye function of order 3.
    
        Parameters
        ----------
        T : float
            Temperature $T$ [K].
    
        Returns
        -------
        float
            Vibrational free energy per atom, $F_\\mathrm{vib}(T)$ [eV/atom].
    
        Notes
        -----
        - Uses ``kB = 8.617333262145e-5`` eV/K.
        - Calls ``calculate_Debye_temperature()`` internally; ensure the material
          parameters needed for $\\Theta_D$ are set.
        - The result is per atom, consistent with the class unit conventions.
        """
        kB = 8.617333262145e-5  # eV K^-1
        from numpy import log, exp
    
        self.calculate_Debye_temperature()
        ThetaD = self.data['ThetaD']
        debInt = self.DebyeIntegral(ThetaD / T)
    
        return 9/8 * kB * ThetaD - kB * T * (debInt - 3 * log(1 - exp(-ThetaD / T)))
    
    
    def F_harm(self, T: float) -> float:
        """
        Total Helmholtz free energy in the harmonic Debye model, per atom.
    
        The total free energy at fixed volume in the harmonic limit is
        $$
        F(T) \\,=\\, E_0 \\, + \\, F_\\mathrm{vib}(T).
        $$
    
        Parameters
        ----------
        T : float
            Temperature $T$ [K].
    
        Returns
        -------
        float
            Total Helmholtz free energy per atom, $F(T)$ [eV/atom].
    
        Notes
        -----
        - Requires that ``E0`` is set (equilibrium energy per atom), and that the
          parameters needed to compute $\\Theta_D$ are available.
        - Internally calls ``Fvib_harm(T)``.
        """
        return self.data['E0'] + self.Fvib_harm(T)




class DebyeQHA:
    """
    Quasi-harmonic Debye model container for equation-of-state parameters.

    This class stores per-atom reference quantities and material information
    needed for quasi-harmonic (volume-dependent) thermodynamic calculations,
    typically combining the Debye model for vibrational contributions with an
    equation of state (e.g., Birch–Murnaghan) for the static lattice energy.

    Parameters
    ----------
    E0 : float
        Reference (equilibrium) energy per atom, $E_0$ [eV/atom].
    V0 : float
        Reference (equilibrium) volume per atom, $V_0$ [Å^3/atom].
    B0 : float
        Bulk modulus at $V_0$, $B_0$ [GPa].
    Bp : float
        Pressure derivative of the bulk modulus at $V_0$, $B'_0$ [dimensionless].
    material : str or ase.Atoms
        Chemical formula string (e.g., 'Al', 'Fe2O3') or an ASE ``Atoms`` object
        describing the material for which thermodynamic properties are computed.
    at_per_fu : int, optional
        Number of atoms in the primitive cell, $N_\\mathrm{at}$ [1], by default 1.

    Attributes
    ----------
    data : dict
        Internal storage for model parameters:
        - 'E0': float
            $E_0$ [eV/atom].
        - 'V0': float
            $V_0$ [Å^3/atom].
        - 'B0': float
            $B_0$ [GPa].
        - 'Bp': float
            $B'_0$ [1].
        - 'material': str or ase.Atoms
            Material descriptor (formula or ``Atoms``).
        - 'Nat': int
            $N_\\mathrm{at}$ [1].
        - 'volumes': sequence of float or None
            Mesh of volumes per atom, $\\{V\\}$ [Å^3/atom], for QHA scans.

    Notes
    -----
    - Units are per atom throughout to keep a consistent convention.
    - The pair $(B_0, B'_0)$ is commonly used with the third-order Birch–Murnaghan
      equation of state to describe $E(V)$ around $V_0$.
    """

    def __init__(
        self,
        E0: float,
        V0: float,
        B0: float,
        Bp: float,
        material: Union[str, Atoms],
        at_per_fu: int = 1
    ) -> None:
        """
        Initialize the DebyeQHA container with equation-of-state parameters.

        See class docstring for a detailed description of parameters and units.
        """
        self.data = {
            'E0': E0,           # float: eV/atom
            'V0': V0,           # float: Å^3/atom
            'B0': B0,           # float: GPa
            'Bp': Bp,           # float: 1 (dimensionless)
            'nu': None,         # float: 1 (dimensionless)
            'material': material,
            'Nat': at_per_fu,   # int: atoms in primitive cell
            'volumes': None     # Sequence[float]: Å^3/atom
        }
        self._deb = []          # array of Debye (HA) objects

    # Nat [atoms in primitive cell]
    @property
    def Nat(self) -> int:
        """
        Number of atoms in the primitive cell, $N_\\mathrm{at}$.
    
        Returns
        -------
        int
            $N_\\mathrm{at}$ [1].
        """
        return self.data['Nat']
    
    @Nat.setter
    def Nat(self, Nat: int) -> None:
        """
        Set the number of atoms in the primitive cell.
    
        Parameters
        ----------
        Nat : int
            $N_\\mathrm{at}$ [1].
    
        Returns
        -------
        None
    
        Raises
        ------
        ValueError
            If ``Nat`` is not a positive integer.
        """
        if not isinstance(Nat, int) or Nat <= 0:
            raise ValueError("Nat must be a positive integer")
        self.data['Nat'] = Nat

    # E0 [eV/atom]
    @property
    def E0(self) -> float:
        """
        Equilibrium energy per atom, $E_0$.
    
        Returns
        -------
        float
            $E_0$ [eV/atom].
        """
        return self.data['E0']
    
    @E0.setter
    def E0(self, E0: float) -> None:
        """
        Set the equilibrium energy per atom.
    
        Parameters
        ----------
        E0 : float
            $E_0$ [eV/atom].
    
        Returns
        -------
        None
        """
        self.data['E0'] = float(E0)


    # V0 [Å^3/atom]
    @property
    def V0(self) -> float:
        """
        Equilibrium volume per atom, $V_0$.
    
        Returns
        -------
        float
            $V_0$ [Å^3/atom].
        """
        return self.data['V0']
    
    @V0.setter
    def V0(self, V0: float) -> None:
        """
        Set the equilibrium volume per atom.
    
        Parameters
        ----------
        V0 : float
            $V_0$ [Å^3/atom].
    
        Returns
        -------
        None
        """
        self.data['V0'] = float(V0)


    # B0 [GPa]
    @property
    def B0(self) -> float:
        """
        Bulk modulus at $V_0$, $B_0$.
    
        Returns
        -------
        float
            $B_0$ [GPa].
        """
        return self.data['B0']
    
    @B0.setter
    def B0(self, B0: float) -> None:
        """
        Set the bulk modulus at $V_0$.
    
        Parameters
        ----------
        B0 : float
            $B_0$ [GPa].
    
        Returns
        -------
        None
        """
        self.data['B0'] = float(B0)


    # Bp [dimensionless]
    @property
    def Bp(self) -> float:
        """
        Pressure derivative of the bulk modulus at $V_0$, $B'_0$.
    
        Returns
        -------
        float
            $B'_0$ [dimensionless].
        """
        return self.data['Bp']
    
    @Bp.setter
    def Bp(self, Bp: float) -> None:
        """
        Set the pressure derivative of the bulk modulus at $V_0$.
    
        Parameters
        ----------
        Bp : float
            $B'_0$ [dimensionless].
    
        Returns
        -------
        None
        """
        self.data['Bp'] = float(Bp)


    # material [formula or ASE Atoms]
    @property
    def material(self) -> Union[str, Atoms]:
        """
        Material descriptor.
    
        Returns
        -------
        str or ase.Atoms
            Chemical formula string or an ASE ``Atoms`` object.
        """
        return self.data['material']
    
    @material.setter
    def material(self, material: Union[str, Atoms]) -> None:
        """
        Set the material descriptor.
    
        Parameters
        ----------
        material : str or ase.Atoms
            Chemical formula string (e.g., 'Al', 'Fe2O3') or an ASE ``Atoms`` object.
    
        Returns
        -------
        None
    
        Raises
        ------
        TypeError
            If ``material`` is not a string or ``ase.Atoms`` instance.
        """
        if not isinstance(material, (str, Atoms)):
            raise TypeError("material must be a chemical formula string or an ase.Atoms object")
        self.data['material'] = material

    # Poisson's ratio
    @property
    def nu(self) -> Optional[float]:
        """
        Poisson's ratio, $\\nu$.

        Returns
        -------
        float or None
            Poisson's ratio $\\nu$ [dimensionless], or ``None`` if unset.

        Notes
        -----
        - For mechanically stable, isotropic media, $-1 < \\nu < 0.5$.
        """
        return self.data['nu']

    @nu.setter
    def nu(self, nu: float) -> None:
        """
        Set Poisson's ratio.

        Parameters
        ----------
        nu : float
            Poisson's ratio $\\nu$ [dimensionless].
        """
        self.data['nu'] = nu

    # volumes [Å^3/atom] (array-like)
    @property
    def volumes(self) -> np.ndarray:
        """
        Volume mesh per atom, {V}.
    
        Returns
        -------
        numpy.ndarray or None
            1D array of volumes per atom [Å^3/atom], or ``None`` if unset.
    
        Notes
        -----
        - Used to sample the quasi-harmonic free energy across different volumes.
        - Stored internally in ``self.data['volumes']``.
        """
        return self.data['volumes']
    
    
    @volumes.setter
    def volumes(self, volumes: Any) -> None:
        """
        Set the volume mesh per atom.
    
        Parameters
        ----------
        volumes : array-like
            Any array-like of volumes per atom [Å^3/atom]. Examples include
            list, tuple, numpy.ndarray, pandas Series.
    
        Returns
        -------
        None
    
        Raises
        ------
        ValueError
            If the provided mesh cannot be converted to a 1D float array or is empty.
    
        Notes
        -----
        - The input is converted to a 1D numpy array of dtype float and stored in
          ``self.data['volumes']``.
        - Ensure values are per atom and, typically, positive.
        """
        arr = np.asarray(volumes, dtype=float)
        if arr.ndim != 1 or arr.size == 0:
            raise ValueError("volumes must be a non-empty 1D array-like of floats.")
        self.data['volumes'] = arr
        self._deb = []
        for V in self.data['volumes']:
            d = Debye()
            d.nu = self.nu
            d.V = V
            d.E0 = self.get_E0(V)
            d.B = self.get_bulk_modulus(V)
            d.rho = get_rho(V, self.material)
            self._deb.append(deepcopy(d))

    
    def get_E0(self, V: float) -> float:
        """
        Energy at volume from stored Birch–Murnaghan parameters.
    
        Evaluates the energy per atom at a given per-atom volume using the
        Birch–Murnaghan equation of state with the parameters stored in the instance.
    
        Parameters
        ----------
        V : float
            Volume per atom, $V$ [Å^3/atom].
    
        Returns
        -------
        float
            Energy per atom, $E(V)$ [eV/atom].
    
        Notes
        -----
        - Uses the instance's ``E0``, ``V0``, ``B0`` (in GPa), and ``Bp`` values.
        - Delegates to ``BM_EoS(V, E0, V0, B0, Bp)``.
        """
        return self.BM_EoS(V, self.E0, self.V0, self.B0, self.Bp)
    
    
    def BM_EoS(self, V: float, E0: float, V0: float, B0: float, B0p: float) -> float:
        """
        Third-order Birch–Murnaghan equation of state, energy vs. volume per atom.
    
        Computes the per-atom energy as a function of per-atom volume using:
        $$
        E(V) \;=\; E_0 \;+\; \frac{9\,V_0\,B_0}{16}
        \\left\\{
          \\left[\\left(\\frac{V_0}{V}\\right)^{2/3}-1\\right]^3 B'_0
          +
          \\left[\\left(\\frac{V_0}{V}\\right)^{2/3}-1\\right]^2
          \\left[6 - 4 \\left(\\frac{V_0}{V}\\right)^{2/3}\\right]
        \\right\\}.
        $$
    
        Parameters
        ----------
        V : float
            Volume per atom, $V$ [Å^3/atom].
        E0 : float
            Reference energy per atom, $E_0$ [eV/atom].
        V0 : float
            Reference volume per atom, $V_0$ [Å^3/atom].
        B0 : float
            Bulk modulus at $V_0$, $B_0$ [GPa].
        B0p : float
            Pressure derivative of the bulk modulus at $V_0$, $B'_0$ [dimensionless].
    
        Returns
        -------
        float
            Energy per atom, $E(V)$ [eV/atom].
    
        Notes
        -----
        - If a consistent energy unit is required, convert $B_0$ to [eV/Å^3] via
          $1\\,\\text{eV}/\\text{Å}^3 = 160.21766208\\,\\text{GPa}$ before use.
        """
        return E0 + 9*V0*B0/16 * (((V0/V)**(2/3)-1)**3 * B0p + ((V0/V)**(2/3)-1)**2* (6-4*(V0/V)**(2/3)))
    
    
    def fit_BM_EoS(self, vols, enes):
        """
        Fit Birch–Murnaghan (third-order) EoS parameters to energy–volume data.
    
        Performs a nonlinear least-squares fit of the BM equation of state to the
        provided per-atom volume and energy data, returning fitted parameters.
    
        Parameters
        ----------
        vols : array-like
            Volumes per atom, $\\{V\\}$ [Å^3/atom].
        enes : array-like
            Energies per atom, $\\{E\\}$ [eV/atom].
    
        Returns
        -------
        numpy.ndarray
            Array of fitted parameters with shape (4,), ordered as:
            [E0, V0, B0, Bp], where
            - E0: float, E0 [eV/atom]
            - V0: float, V0 [Å^3/atom]
            - B0: float, B0 [GPa]
            - Bp: float, B0' [dimensionless]
    
        Notes
        -----
        - Initial guesses used for the fit are:
          $E_0 = \\min(E)$,
          $V_0 = \\tfrac{1}{2}(\\min(V)+\\max(V))$,
          $B_0 \\approx 1.5\\,\\text{eV}/\\text{Å}^3$,
          $B'_0 = 4$.
        - Converts the fitted bulk modulus from [eV/Å^3] to [GPa] using
          $1\\,\\text{eV}/\\text{Å}^3 \\approx 160.2\\,\\text{GPa}$.
        - Assumes a callable ``EoS(V, E0, V0, B, Bp)`` compatible with ``curve_fit``.
        """
        eosfit, _ = curve_fit(self.BM_EoS, vols, enes, [min(enes), 0.5*(min(vols)+max(vols)), 1.5, 4])
        eosfit[2] *= 160.2
        return eosfit

        
    import numpy as np
    from numpy.typing import ArrayLike, NDArray
    
    def calc_equilibrium(self, T: ArrayLike) -> NDArray[np.floating]:
        """
        Fit BM3 EoS at each temperature and return the fitted parameters.
    
        For each temperature in T, computes the harmonic free energies across the
        stored volume mesh (via self._deb[i].F_harm(T_i)), fits the third-order
        Birch–Murnaghan EoS to E(V), and returns the fitted parameters.
    
        Parameters
        ----------
        T : array-like
            Temperatures, {T} [K]. Can be a scalar or array-like.
    
        Returns
        -------
        numpy.ndarray
            - If T is scalar: shape (4,) ordered as [E0, V0, B0, Bp]
            - If T is array-like: shape (N, 4) where N = len(T)
              ordered as [E0, V0, B0, Bp] for each temperature.
    
        Notes
        -----
        - Assumes self.volumes is the per-atom volume mesh [Å^3/atom].
        - Assumes self._deb is an iterable of Debye objects aligned with self.volumes,
          each providing F_harm(T) -> energy [eV/atom] at temperature T.
        """
        is_scalar = np.isscalar(T) or (np.ndim(T) == 0)
        T_arr = np.atleast_1d(T).astype(float)
    
        fits = []
        for t in T_arr:
            enes = [d.F_harm(t) for d in self._deb]
            eosfit = self.fit_BM_EoS(self.volumes, enes)  # returns [E0, V0, B0, Bp]
            fits.append(eosfit)
    
        fits = np.asarray(fits, dtype=float)
        return fits[0] if is_scalar else fits
    
    
    def calc_F(self, T: ArrayLike) -> NDArray[np.floating] | float:
        """
        Equilibrium energy E0(T) from BM3 fit at each temperature.
    
        Parameters
        ----------
        T : array-like
            Temperatures, {T} [K]. Can be a scalar or array-like.
    
        Returns
        -------
        float or numpy.ndarray
            - If T is scalar: scalar float E0(T) [eV/atom]
            - If T is array-like: array of E0(T) [eV/atom] with shape (N,)
        """
        eosfits = self.calc_equilibrium(T)
        if eosfits.ndim == 1:
            return float(eosfits[0])
        return eosfits[:, 0]
    
    
    def calc_V(self, T: ArrayLike) -> NDArray[np.floating] | float:
        """
        Equilibrium volume V0(T) from BM3 fit at each temperature.
    
        Parameters
        ----------
        T : array-like
            Temperatures, {T} [K]. Can be a scalar or array-like.
    
        Returns
        -------
        float or numpy.ndarray
            - If T is scalar: scalar float V0(T) [Å^3/atom]
            - If T is array-like: array of V0(T) [Å^3/atom] with shape (N,)
        """
        eosfits = self.calc_equilibrium(T)
        if eosfits.ndim == 1:
            return float(eosfits[1])
        return eosfits[:, 1]
    
    
    def calc_B(self, T: ArrayLike) -> NDArray[np.floating] | float:
        """
        Equilibrium bulk modulus B0(T) from BM3 fit at each temperature.
    
        Parameters
        ----------
        T : array-like
            Temperatures, {T} [K]. Can be a scalar or array-like.
    
        Returns
        -------
        float or numpy.ndarray
            - If T is scalar: scalar float B0(T) [GPa]
            - If T is array-like: array of B0(T) [GPa] with shape (N,)
        """
        eosfits = self.calc_equilibrium(T)
        if eosfits.ndim == 1:
            return float(eosfits[2])
        return eosfits[:, 2]

    
    def get_pressure(self, V: float) -> float:
        """
        Third-order Birch–Murnaghan pressure, P(V), as a function of per-atom volume.
    
        Computes the pressure using the BM equation of state:
        $$
        P(V) \;=\; \frac{3}{2} B_0
          \left[
            \left(\frac{V_0}{V}\right)^{7/3} - \left(\frac{V_0}{V}\right)^{5/3}
          \right]
          \left\{
            1 + \frac{3}{4}\,(B'_0 - 4)\left[\left(\frac{V_0}{V}\right)^{2/3} - 1\right]
          \right\}.
        $$
    
        Parameters
        ----------
        V : float
            Volume per atom, $V$ [Å^3/atom].
    
        Returns
        -------
        float
            Pressure, $P(V)$ [GPa].
    
        Notes
        -----
        - Units:
          - ``V0`` in [Å^3/atom], ``B0`` in [GPa], ``Bp`` dimensionless.
          - The returned pressure is in [GPa], consistent with ``B0``.
        """
        V0 = self.V0
        B0 = self.B0
        Bp = self.Bp
    
        eta = V0 / V
        return (3.0 * B0 / 2.0) * (eta**(7.0/3.0) - eta**(5.0/3.0)) * (1.0 + 0.75 * (Bp - 4.0) * (eta**(2.0/3.0) - 1.0))


    def get_bulk_modulus(self, V: float) -> float:
        """
        Volume-dependent (isothermal) bulk modulus, B(V), in the BM3 framework.
    
        Uses the linearized relation based on the definition of the pressure derivative
        at zero pressure:
        $$
        B(V) \;\approx\; B_0 \;+\; B'_0\,P(V),
        $$
        where $P(V)$ is the Birch–Murnaghan (third-order) pressure and $B'_0 = (dB/dP)_{P=0}$.
    
        Parameters
        ----------
        V : float
            Volume per atom, $V$ [Å^3/atom].
    
        Returns
        -------
        float
            Bulk modulus, $B(V)$ [GPa].
    
        Notes
        -----
        - This is a linear approximation in pressure around $P=0$; for large compressions
          or expansions, a full BM3 expression for $B(V)$ may be preferred.
        - Units: ``B0`` in [GPa], ``Bp`` dimensionless, returned $B(V)$ in [GPa].
        """
        p = self.get_pressure(V)
        return self.B0 + self.Bp * p