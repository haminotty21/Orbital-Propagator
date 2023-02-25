from datetime import datetime
from math import copysign, pow, pi

from numpy import arccos, dot, array, cross
from numpy.linalg import norm


def get_julian_datetime(date):
    """
    Convert a datetime object into julian float.
    Args:
        date: datetime-object of date in question

    Returns: float - Julian calculated datetime.
    Raises:
        TypeError : Incorrect parameter type
        ValueError: Date out of range of equation
    """

    # Ensure correct format
    if not isinstance(date, datetime):
        raise TypeError('Invalid type for parameter "date" - expecting datetime')
    elif date.year < 1801 or date.year > 2099:
        raise ValueError('Datetime must be between year 1801 and 2099')

    # Perform the calculation
    julian_datetime = 367 * date.year - int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) + int(
        (275 * date.month) / 9.0) + date.day + 1721013.5 + (
                          date.hour + date.minute / 60.0 + date.second / pow(60,
                                                                                  2)) / 24.0 - 0.5 * copysign(
        1, 100 * date.year + date.month - 190002.5) + 0.5

    return julian_datetime




def SI_to_AU(r=None, v=None, mu=None):
    if r is not None:
        r = r/149597871
    if v is not None:
        v = v/149597871 * 3600 * 24
    if mu is not None:
        mu = mu * (3600*24)**2/(149597871)**3
    return r, v, mu

def AU_to_SI(r=None, v=None, mu=None):
    if r is not None:
        r = r*149597871

    if v is not None:
        v = v * 149597871 / (3600 * 24)
    if mu is not None:
        mu = mu / (3600*24)**2 * 149597871**3
    return r, v, mu

def orbital_elements(r, v, mu):
    r_m = norm(r)
    v_m = norm(v)

    v_r = dot(r, v/r_m)

    h = cross(r, v)          # Specific Angular momentum (1st Element)
    h_m = norm(h)
    i = arccos(h[2]/h_m)    # Inclination (2nd Element)

    if i > pi or i < 0:
        print(f'Inclination is {i} and should be within range 0 to pi')

    N = cross(array([0, 0, 1]), h)
    N_m = norm(N)
    raan = arccos(N[0] / N_m)             # Right Ascension of the ascending node (3rd Element)

    if N[1] < 0:
        2*pi - raan

    e = 1/mu*((v_m**2 - mu/r_m)*r - r_m*v_r*v)      # Eccentricity (4th Element)
    e_m = norm(e)
    aop = arccos(dot(N/N_m, e/e_m))     # Argument of Perigee (5th)

    if e[2] < 0:
        2*pi - aop

    theta = arccos(dot(e/e_m, r/r_m))               # True Anomaly (6th)

    if v_r < 0:
        2*pi - theta

    return h, i, raan, e, aop, theta
