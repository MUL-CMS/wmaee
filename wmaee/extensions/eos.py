from scipy.optimize import curve_fit
from numpy import min, max

def BM_EoS(V, E0, V0, B0, B0p):
    """Definition of Birch-Murnaghan EoS

    :param V: volume at which we want energy
    :type V: float or np.array
    :param E0: equilibrium energy
    :type E0: float
    :param V0: equilibrium volume
    :type V0: float
    :param B0: bulk modulus
    :type B0: float
    :param B0p: pressure derivative of bulk modulus
    :type B0p: float
    :return: energy at volume(s) V
    :rtype: float or np.array
    """
    return E0 + 9*V0*B0/16 * (((V0/V)**(2/3)-1)**3 * B0p + ((V0/V)**(2/3)-1)**2 * (6-4*(V0/V)**(2/3)))

def fit_EV_curve(volumes, energies, eos='BM', p0=None):
    """Perform energy-volume curve fitting.

    :param volumes: list of volumes
    :type volumes: list or np.array
    :param energies: energies
    :type energies: list or np.array
    :param eos: equation of state. Allowed values are:
        'BM' (Birch-Murnaghan EoS). 
        Defaults to 'BM'
    :type eos: str, optional
    :param p0: set of initial values for EoS fitting, defaults to None
    :type p0: tuple, optional
    :raises Exception: unknown type of EoS
    :return: fitted EoS parameters in the form of list and dictionary
    :rtype: (list, dict)
    """
    if eos == 'BM':
        if p0 == None:
            # some usually good initial values for the fit            
            p0 = (min(energies), # guess for E0
                  0.5*(min(volumes)+max(volumes)), #mean value of volume range
                  1, # B0~160GPa
                  4)
        # perform Birch-Murnaghan fit              
        optimum, _ = curve_fit(BM_EoS, volumes, energies, p0=p0)
        opt = {'E0': optimum[0],
               'V0': optimum[1],
               'B0': optimum[2]*160.2,
               'B0p': optimum[3]
               }
    else:
        raise Exception('Unknown EoS.')
    return optimum, opt